# Parametrizing new chemical entities for AMBER

## Ligands and modified residues with AmberTools

---

## About this tutorial

Sooner or later every molecular dynamics project runs into a molecule that the
force field does not know: a ligand that is not in the standard libraries, a
cofactor, a post-translational modification, a non-natural amino acid. This
tutorial is about what to do then.

By the end you will have produced, for a molecule of your choice, the two files
that are the actual deliverable of any parametrization job:

- **a library file** — `.lib` for a ligand, `.prepi` for a residue — which tells
  LEaP what atoms the molecule has, what types they are, what charges they
  carry, and (for a residue) how it connects to its neighbours;
- **a `.frcmod` file** — which supplies the bonded and non-bonded parameters
  that the parent force field is missing.

Everything else you will produce along the way is scaffolding.

**Two worked examples run through the whole document.** For ligands we use
**ceramide (CER)**. For modified residues we use **methionine sulfoxide (SME)**,
the oxidised form of methionine. The two tracks share most of their steps; where
they diverge the section headings say so explicitly.

> **This is not a black box.** Every command below will produce *some* output
> even when the result is wrong. The checkpoints scattered through the text are
> there because the difference between a parametrization that works and one that
> quietly poisons your simulation is almost always a line of output that nobody
> read. Please read them.

### What you need

- AmberTools (any recent version; the output shown here is from AmberTools 24).
  `antechamber`, `parmchk2`, `prepgen` and `tleap` all ship with it and are free.
- A molecular visualizer. We use PyMOL for structure editing and VMD for
  checking results, but anything equivalent works.
- A correct starting structure of your molecule, and two numbers you must know
  before you start: its **net charge** and its **spin multiplicity** (almost
  always 1 — but if your molecule is a radical or has an open-shell metal
  centre, stop and reconsider whether this workflow applies at all).

### Conventions used here

Commands to type are shown in boxes. Output you should expect is shown
indented, usually abbreviated. Whenever we paste a long log, the sentence
immediately after it tells you **which single line matters** — ignore the rest.

Anything written in `ANGLE BRACKETS` is a placeholder for you to replace.
Numbers that are specific to our examples — net charges, ion counts, atom names
— are flagged as such. **Do not copy them into your own system without
thinking.**

We write file-editing steps with `vi` because it exists on every machine. If you
prefer `nano`, `emacs`, VS Code or anything else, use it. For `vi` the two
commands you need are `i` to start typing and `:wq` to save and quit
(`:q!` quits without saving). A short cheat sheet:
<https://www.atmos.albany.edu/daes/atmclasses/atm350/vi_cheat_sheet.pdf>

If you get stuck anywhere in this procedure, contact us.

---

# Part 0 — What a force field actually needs

Before typing anything, it is worth being precise about what we are producing,
because the rest of the tutorial is just four shopping-list items.

A classical force field computes the energy of a molecule from a fixed functional
form. To evaluate it for *your* molecule, it needs:

1. **Atom types** — a label on every atom (`c3`, `oh`, `CX`, …) that says which
   row of the parameter tables applies to it. Types encode element *and* chemical
   environment: an sp³ carbon and an aromatic carbon are different types.
2. **Partial charges** — one number per atom, for the Coulomb term. These must
   sum to the formal charge of the molecule, exactly.
3. **Bonded parameters** — equilibrium bond lengths and force constants, angles,
   proper and improper torsions, for every combination of types that occurs in
   your molecule.
4. **Non-bonded parameters** — Lennard-Jones radius and well depth per type.

Items 1 and 2 are produced by `antechamber`. Items 3 and 4 are looked up, or
guessed by analogy, by `parmchk2`. That is the whole ligand workflow, and it is
why the ligand branch of this tutorial is short.

**A residue needs a fifth thing: connectivity to its neighbours.** A ligand is a
closed molecule. A residue is a unit of a polymer — it has a head that bonds to
the previous residue and a tail that bonds to the next one, and LEaP has to know
which atoms those are and what types the neighbouring atoms will have. Nothing in
a `.mol2` file can express that. This is the single reason the residue branch is
longer, and it is worth holding on to, because every extra step in that branch
exists to supply this one missing piece of information.

---

# Part 1 — Step 0: check whether someone has already done it

Parametrization is work and it is error-prone. Before doing any of it, look.

```bash
ls $AMBERHOME/dat/leap/lib          # standard residue and ion libraries
ls $AMBERHOME/dat/leap/prep         # prepi templates, including amino12.in
ls $AMBERHOME/dat/leap/parm         # parameter files: parm10.dat, gaff2.dat, frcmod.*
ls $AMBERHOME/dat/contrib           # contributed parameter sets
```

Also worth checking before you start:

- **The Bryce Group AMBER Parameter Database** (University of Manchester), which
  collects published parameter sets for cofactors, modified amino acids,
  nucleotides and common ligands.
- **The literature for your specific molecule.** Phosphorylated and methylated
  amino acids, haem, nucleotide cofactors and the common lipids all have
  published, validated parameter sets. A published set that has been tested
  against experiment beats anything you will produce in an afternoon.

If you find one, use it, and cite it. If you do not, continue.

---

# Part 2 — Step 1: the starting structure

This is the step people rush, and it is the step that determines whether
everything downstream is meaningful. The charges you are about to compute
describe *the structure you hand in* — its protonation, its stereochemistry, its
conformation. Nothing later in the pipeline will fix a mistake made here.

## 2.1 Requirements common to both branches

- **All hydrogens present and explicit.**
- **Correct protonation state** for the pH you intend to simulate. If in doubt,
  use a tool: **Dimorphite-DL** for ligands on their own; **PDB2PQR**, **ProToss**
  or **H++** for a ligand or residue in its protein context. The answer can depend
  on the environment, which is why the protein-context tools exist.
- **Correct stereochemistry**, including at any centre your modification creates.
  Methionine sulfoxide, our example, has a stereogenic sulfur — the R and S
  sulfoxides are different molecules and you must decide which one you are
  parametrizing.
- **A sensible conformation.** See the box in Part 3 on why this matters more
  than it looks.
- **Open the structure in a visualizer and look at it.** Every time. Most
  parametrization problems are visible in ten seconds at this stage.

## 2.2 Ligand branch — CER

For a ligand, the structure you prepare *is* the molecule you are parametrizing.
Save it as PDB (other formats work; `antechamber` reads many of them). We start
from `CER.pdb`.

Two practical points:

- The file should contain **only** your ligand.
- **Note the atom names now.** They must match, exactly, the atom names used for
  this ligand in whatever PDB you later load into LEaP.
This is the single most common cause of LEaP failures later.

## 2.3 Modified-residue branch — SME

Here the structure you prepare is **not** the thing you want. You want an
internal residue of a protein chain, but what you can download is an isolated
molecule with a free amine and a free carboxylic acid. Those termini are
chemically wrong for a chain unit, and if you compute charges on them, the
charges on the nearby backbone atoms will be wrong too.

So we get the fragment and then we **cap it** so it looks like a piece of
a chain.

### Getting the fragment

If your protein already contains the modified residue, you can cut it out of your
own PDB. Otherwise, download it from the PDB chemical component dictionary:

1. Find the residue's three-character code in your PDB file (it is alphanumeric;
   ours is `SME`).
2. Go to <https://www.rcsb.org/>, click **Chemical Search** below the search bar.
3. Choose the **Chemical Component** tab, enter the code, hit search.
4. Download the component file. We save `SME.cif`.

> **Use the CCD (Chemical Component Dictionary) as the reference on atom names.** Whatever you build, your final
> `.prepi` must use the same atom names as the PDB files you will later load. For
> SME the CCD names are:
>
> `N CA CB CG S OE CE C O OXT H H2 HA HB2 HB3 HG2 HG3 HE1 HE2 HE3 HXT`
>
> Note that **the sulfur is called `S`, not `SD`.** In standard methionine it is
> `SD`. This difference will matter in Part 6, Route 2 — remember it.

### Capping: ACE–SME–NME

We delete the two atoms that make the fragment a free amino acid, and replace
them with an acetyl group (ACE) at the N-terminus and an N-methylamide group
(NME) at the C-terminus. The result mimics the chemical environment of an
internal residue.

In PyMOL:

1. Open `SME.cif`.
2. Delete **`H2`** (the second hydrogen on N) and **`OXT` + `HXT`** (the hydroxyl
   of the carboxylic acid). Select atoms with `Ctrl` + middle mouse button, then
   right-click → *atom* → *delete*.
3. Select the N-terminal nitrogen, then **Build → Residue → Acetyl**.
4. Select the C-terminal carbon, then **Build → Residue → N-Methyl**.
5. **File → Export Molecule**, choose PDB, save as `SME.pdb`.

Any tool that can add ACE and NME caps will do; there is nothing special about
PyMOL here.

> **Checkpoint.** Open `SME.pdb` and confirm: the caps are there, the sulfoxide
> oxygen is there, the stereochemistry is what you intended, and the geometry is
> not distorted. Count the atoms — ours has 30.

> **The net charge you will need next is the charge of the *residue*, not of the
> capped molecule.** ACE and NME are neutral, so the two are the same number.
> For SME that number is 0. For phosphoserine it would be −2. Getting this
> wrong is silent and fatal: the charges will be computed for a different
> chemical species than the one you meant.

---

# Part 3 — Step 2: charges and atom types

Both branches now run the same program, which does both jobs at once.

## Ligand

```bash
antechamber -i CER.pdb -fi pdb -o CER.mol2 -fo mol2 -at gaff2 -nc 0 -rn CER -c bcc -s 2
```

## Modified residue

```bash
antechamber -i SME.pdb -fi pdb -o SME.mol2 -fo mol2 -at gaff2 -nc 0 -rn SME -c bcc -s 2
```

(The `-nc 0` in both lines is **our** net charge. Use yours.)

### What the flags mean

| flag | meaning |
|---|---|
| `-i` / `-fi` | input file and its format |
| `-o` / `-fo` | output file and its format |
| `-at` | which **atom type** set to assign. Default `gaff`; we ask for `gaff2` |
| `-c` | charge method. `bcc` = AM1-BCC |
| `-nc` | net charge of the input molecule |
| `-rn` | residue name to write into the output. **Maximum 3 characters** |
| `-s` | how verbose the printout is (`0`, `1`, `2`) |
| `-pf y` | optional: delete the `ANTECHAMBER_*` intermediate files on the way out |

`antechamber -help` lists everything.

> **Warning about `-s`.** In `antechamber`, `-s 2` means "verbose". In
> `parmchk2`, which we run in Part 7, `-s 2` means "use GAFF2". Same letter,
> completely different meaning.

### What is happening under the hood

`antechamber` runs several programs in sequence; you will see them scroll past:

```
Running: .../bondtype -j full -i ANTECHAMBER_BOND_TYPE.AC0 ...
Running: .../atomtype -i ANTECHAMBER_AC.AC0 -o ANTECHAMBER_AC.AC -p gaff2
Info: Total number of electrons: 318; net charge: 0
Running: .../sqm -O -i sqm.in -o sqm.out
Running: .../am1bcc -i ANTECHAMBER_AM1BCC_PRE.AC ...
```

`bondtype` perceives bond orders, `atomtype` assigns GAFF2 types, `sqm` does the
semi-empirical quantum calculation, and `am1bcc` applies the bond-charge
corrections to the resulting AM1 charges.

> **Line to check:** `Total number of electrons: 318; net charge: 0`. If that net
> charge is not the one you intended, stop — everything downstream is wrong.

### Checkpoint: did the quantum calculation converge?

```bash
tail sqm.out
```

```
 QMMM: 109      109      H       -33.6246  170.6083    0.8407
 QMMM: 110      110      H       -33.3839  172.2030    0.0149
 QMMM: 111      111      H       -34.3504  172.1117    1.5442
          --------- Calculation Completed ----------
```

> **Line to check:** `Calculation Completed`. If it is absent, `sqm` failed and
> your charges are garbage or missing. See Appendix A.

Once you have seen that line, the `ANTECHAMBER_*` and `sqm.*` files are
intermediates and can be deleted.

### Reading the `.mol2`

```bash
vi CER.mol2
```

```
@<TRIPOS>ATOM
      1 O1         -33.0660  172.4520   20.1460 oh     1542 CER       -0.610800
      2 HO1        -33.5210  172.3850   20.9890 ho     1542 CER        0.412000
      3 C1S        -33.6690  171.5530   19.1750 c3     1542 CER        0.110400
```

Columns, left to right: atom number, **atom name**, x, y, z, **atom type**,
residue number, residue name, **partial charge**.

Those three bolded columns are the entire product of this step.

> **Checkpoint: do the charges sum to your formal charge?**
> ```bash
> python3 check_charge.py CER.mol2 0
> ```
> (`check_charge.py` is provided with this tutorial. The expected charge is the
> second argument.) They should sum to your net charge to within rounding.
> Checking now costs five seconds; discovering it after a week of simulation
> does not.

> **Checkpoint: do the atom names match your target PDB?** If this ligand or
> residue appears in a PDB you will later load into LEaP, the names in column 2
> must match that file exactly. Fix it now, in whichever file is wrong.

---

> ## Box: why AM1-BCC and not RESP?
>
> The charges in the standard AMBER protein force fields (ff94 through ff19SB)
> were derived by **RESP** fitting to an electrostatic potential computed at
> HF/6-31G\*. AM1-BCC is a different and cheaper model: it computes AM1
> charges semi-empirically and then applies empirical bond-charge corrections
> that were fitted to reproduce HF/6-31G\* ESP charges. It is good, it is fast,
> it is free, and it is the standard choice for GAFF ligands.
>
> **We use it here for a practical reason: RESP as originally formulated needs
> Gaussian, and we do not have a Gaussian licence.**
>
> What you give up,
>
> - AM1-BCC approximates HF/6-31G\* ESP charges; it does not reproduce them. A
>   residue you parametrize this way speaks a slightly different electrostatic
>   language from the ff14SB residues around it.
> - The charges depend on the **conformation** you submit. `sqm` performs an AM1
>   geometry optimisation before the charges are derived, so they correspond to
>   the AM1-optimised geometry reached from your input. For a long flexible
>   molecule like our ceramide this is not a detail.
>   **Start from an extended conformation.** If you need charges for one
>   specific geometry, freeze it with `-ek "maxcyc=0"`.
> - Single-conformation charges in general are a known weakness. Published RESP
>   parameter sets usually fit several conformations at once.
>
> **If you want RESP and have no Gaussian, you are not actually blocked.** The
> RESP fitting program itself (`resp`) ships with AmberTools; what you need
> externally is only the electrostatic potential. Free or academically free QM
> codes that can produce it include **GAMESS-US**, **Psi4** and **ORCA**. The
> **R.E.D. Tools / R.E.D. Server** suite automates the whole multi-conformation
> RESP procedure and can drive GAMESS-US or Firefly rather than Gaussian. This
> is more work than `-c bcc`, and out of scope for a single session, but it is
> the right answer if you are producing parameters you intend to publish.
>
> A newer alternative also worth knowing: `-c abcg2` (ABCG2) is a re-fitted
> bond-charge-correction model available in recent AmberTools, intended as a
> successor to AM1-BCC.

---

# Part 4 — Step 3 (modified residues only): fix the backbone atom types

**Ligand track: skip to Part 5.**

We asked `antechamber` for GAFF2 types, and it obediently assigned GAFF2 types to
every atom, including the backbone. That is wrong, and it is wrong for two
separate reasons.

**The immediate reason: LEaP will refuse to build your chain.** When LEaP joins
your residue to its neighbours, it creates a bond between the previous residue's
carbonyl carbon — an ff14SB atom of type `C` — and your residue's nitrogen. If
your nitrogen carries the GAFF2 type `n`, LEaP needs a `C–n` bond parameter.
There is no such parameter, in ff14SB or in GAFF2, and you get:

```
Could not find bond parameter for: C - n
```

**The deeper reason: φ and ψ.** The ff14SB backbone torsions — the parameters
that make a protein backbone behave like a protein backbone — are defined for
the type sequence `N–CX–C–N`. Those parameters are the product of years of
refitting against QM and NMR data, and they are the main thing that distinguishes
one AMBER protein force field generation from the next. If your residue's
backbone carries GAFF types, it gets GAFF torsions instead: generic ones, fitted
for small organic molecules, with no knowledge of protein conformational
preferences. Your modified residue would sample a different Ramachandran
landscape from every other residue in the chain. Even if LEaP built it without
complaining, this would be wrong.

**GAFF is for the part of your molecule that is genuinely unusual.** GAFF and
GAFF2 were developed to cover a broad space of drug-like organic chemistry, which
is exactly what you need for a sulfoxide, a phosphate ester or a non-natural
sidechain. They were not developed to represent protein backbones, and ff14SB
was. Use each where it belongs.

### What to do

Open the file and replace the backbone atom types by hand:

| atom name | GAFF2 assigns | **change it to** |
|---|---|---|
| `N` | `n` | **`N`** |
| `H` | `hn` | **`H`** |
| `CA` | `c3` | **`CX`** |
| `HA` | `h1` | **`HA`** |
| `C` | `c` | **`C`** |
| `O` | `o` | **`O`** |

Everything else — the whole sidechain, and the cap atoms — keeps its GAFF2 type.

Note the case convention, which is not decorative: **AMBER protein types are
upper case, GAFF types are lower case**. `C` and `c` are different types, and so
are `N`/`n`, `O`/`o`, `S`/`s`. This makes mixed files easy to read at a glance.

> **An alternative worth knowing.** You can ask `antechamber` for AMBER types
> throughout with `-at amber`, which is what the official AMBER tutorial on
> modified residues does. Then there is nothing to edit. The trade-off:
> AMBER's type set covers protein-like chemistry well and unusual chemistry
> badly. The decision rule we suggest:
>
> - modification is *protein-like* (methylated lysine, a non-natural amino acid
>   built from standard functional groups) → `-at amber`, nothing to edit;
> - modification is *genuinely unusual chemistry* (sulfoxide, phosphate,
>   halogen, organometallic) → `-at gaff2` plus the backbone fix above.
>
> SME is a sulfoxide, so we take the second route.

---

# Part 5 — Branch A: the ligand library

**Modified-residue track: skip to Part 6.**

For a ligand there is no connectivity problem to solve, so we go straight to
building the library. 

(You will need `CER.frcmod` first — that is Part 7. If you are following
linearly, do Part 7 now and come back. We have kept the library step here
because conceptually it belongs with the ligand branch.)

Create `lib.leap`:

```bash
vi lib.leap
```

```
source leaprc.gaff2

CER = loadmol2 CER.mol2
loadamberparams CER.frcmod

charge CER
check CER

saveoff CER CER.lib
saveAmberParmNetCDF CER CER.prmtop CER.ncrst
savepdb CER CER_check.pdb
quit
```

Run it:

```bash
tleap -f lib.leap
```

The log is long. Three lines matter:

```
Total unperturbed charge:   0.000000
Checking parameters for unit 'CER'.
Unit is OK.
```

> **Line to check 1:** `Total unperturbed charge` must equal your formal charge,
> as an integer.
> **Line to check 2:** `Unit is OK`. If instead you see `Could not find bond
> parameter` or similar, your `.frcmod` is incomplete.

You will also see, harmlessly:

```
 (Residues lacking connect0/connect1 -
  these don't have chain types marked:
     res     total affected
     CER     1
 )
```

That is LEaP noticing your ligand has no head and no tail. For a ligand this is
correct and expected.

> **Checkpoint.** `vmd CER.prmtop CER.ncrst`. Does the molecule look right? Any
> absurd bond lengths, flattened centres, atoms sitting on top of each other?

**You now have `CER.lib` and `CER.frcmod`. Those two files are the deliverable.**
Skip to Part 8.

---

# Part 6 — Branch B: the modified-residue library

## 6.1 Why `.lib` is not enough

Recall Part 0: a residue has to tell LEaP where its head and tail are, and what
atom types it expects its neighbours to present. The file format that carries
this is **prepi** (sometimes `.prepin`), an older AMBER format that differs from
`.mol2` in two relevant ways:

- it stores **internal coordinates** (bond length, angle, dihedral relative to
  three previously-defined atoms) rather than Cartesian coordinates, which is
  what lets LEaP place the residue into an existing chain;
- it carries the head/tail definitions and the improper torsion list.

We will produce it in two different ways, and then compare them — the comparison
is the part of this section you should actually remember.

## 6.2 The prepi format

Here is the ff14SB template for ordinary methionine, from
`$AMBERHOME/dat/leap/prep/amino12.in`:

```
METHIONINE
 MET  INT     1
 CORR OMIT DU   BEG
   0.00000
   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000    0.00000
   2  DUMM  DU    M    1   0  -1     1.449     0.000     0.000    0.00000
   3  DUMM  DU    M    2   1   0     1.522   111.100     0.000    0.00000
   4  N     N     M    3   2   1     1.335   116.600   180.000   -0.41570
   5  H     H     E    4   3   2     1.010   119.800     0.000    0.27190
   6  CA    CX    M    4   3   2     1.449   121.900   180.000   -0.02370
   7  HA    H1    E    6   4   3     1.090   109.500   300.000    0.08800
   8  CB    2C    3    6   4   3     1.525   111.100    60.000    0.03420
   ...
  14  SD    S     S   11   8   6     1.810   110.000   180.000   -0.27370
  15  CE    CT    3   14  11   8     1.780   100.000   180.000   -0.05360
   ...
  19  C     C     M    6   4   3     1.522   111.100   180.000    0.59730
  20  O     O     E   19   6   4     1.229   120.500     0.000   -0.56790

IMPROPER
 -M   CA   N    H
 CA   +M   C    O
DONE
```

Columns: index, **atom name**, **atom type**, **topology letter**, three
reference atom indices, bond length, bond angle, dihedral, **charge**.

### The three dummy atoms

Every prepi starts with three `DUMM` atoms of type `DU`. They exist only to
define a coordinate frame for the first real atom, since internal coordinates
need three prior atoms to refer to. They are required; they are not part of your
molecule; they carry no charge.

### The topology column

This is the column everyone finds mysterious. It describes the **build tree** —
the order in which LEaP reconstructs the molecule from internal coordinates — and
the letter counts **how many children an atom has in that tree, hydrogens
included**:

| letter | meaning | children |
|---|---|---|
| `M` | main chain: on the path from HEAD to TAIL | — |
| `E` | end of a branch | 0 |
| `S` | side chain, linear continuation | 1 |
| `B` | branch | 2 |
| `3` | three-way branch | 3 |

`M` takes precedence: a backbone atom is `M` regardless of how many children it
has.

Two worked examples from the file above, which make the rule concrete:

- `CB` is **`3`**, and its three children are `HB2`, `HB3`, `CG`. **Two of those
  are hydrogens.** The letter counts branches in the tree, *not* heavy-atom
  coordination.
- `SD` is **`S`**: in methionine the sulfur has exactly one child, `CE`.

Now the point that matters for us: **in methionine sulfoxide the same sulfur
gains a second child, `OE`, and therefore must become `B`.** Adding an atom
changes its parent's topology letter. Forgetting this is one of the easy ways to
produce a prepi that builds a mangled molecule.

The official reference is <https://ambermd.org/doc/prep.html>, which also states
the tree-ordering rule:

> The tree structure begins at the first dummy atom, and traverses the main chain
> until a branch point (node) is reached. That branch is traversed until its end
> or until the next node is reached. When you come to a node with more than one
> branch (topological type "B" or "3"), it doesn't matter which branch is
> traversed first as long as you return to the next higher node when an end is
> reached.

### The IMPROPER block

```
 -M   CA   N    H
 CA   +M   C    O
```

`-M` means "the main-chain atom of the **previous** residue" and `+M` "of the
**next** residue". These two entries keep the two peptide groups planar. Unless
your modification is *on* the backbone, copy them as they are.

---

## 6.3 Route 1 — `prepgen` (the recommended working method)

`prepgen` takes your capped molecule, removes the caps, and writes the prepi,
redistributing the caps' charge so the residue comes out with the right total.
It needs the molecule in `.ac` format and a **mainchain file** telling it what to
remove and where the head and tail are.

### Convert to `.ac`

```bash
antechamber -i SME.mol2 -fi mol2 -o SME.ac -fo ac -at gaff2
```

`prepgen` only reads `.ac`, which is why we take this detour.

```bash
vi SME.ac
```

```
CHARGE     -0.00 ( 0 )
Formula: H16 C8 N2 O3 S1
ATOM      1  CZ  SME    -1  -2.967   0.010  -1.116  0.664100        c
ATOM      7  N   SME     0  -1.797  -0.400  -1.472 -0.570900        N
ATOM      8  CA  SME     0  -0.566   0.400  -1.524  0.044700       CX
...
BOND     7    7    8  1      N     CA
```

Same information as the `.mol2`, plus an explicit connectivity block at the
bottom. Residue number `-1` is the ACE cap, `0` is SME, `1` is the NME cap.

> This is the file in which you apply the backbone atom-type fix from Part 4, if
> you have not already done it in the `.mol2`. Check it before continuing:
> `N H CA HA C O` must read `N H CX HA C O` in the last column.

### Rename the cap atoms

`prepgen` identifies atoms **by name**, so every name in the file must be unique
and unambiguous. Straight out of PyMOL they are not: the NME nitrogen is called
`N`, exactly like your residue's backbone nitrogen, and both methyl groups are
called `CH3` with hydrogens `HH31/HH32/HH33`.

Two rules:

1. **No name may be shared between an atom you keep and an atom you omit.** If
   you write `OMIT_NAME N`, `prepgen` cannot tell which `N` you mean.
2. **No name may be shared between two atoms you omit, either.** It is tempting
   to think it does not matter since both disappear, but `prepgen` will not
   accept the ambiguity.

**Rename the caps, not the residue.**
Your residue's backbone nitrogen must stay `N`, because that is what every
PDB in the world calls it and LEaP matches by name. Rename the caps instead:

| cap | original | rename to |
|---|---|---|
| ACE | `C`, `O`, `CH3`, `HH31-33` | `CY`, `OY`, `CAY`, `HY1`, `HY2`, `HY3` |
| NME | `N`, `H`, `CH3`, `HH31-33` | `NZ`, `HZ`, `CAT`, `HT1`, `HT2`, `HT3` |

Nothing belonging to SME changes, so nothing has to be renamed back.

> The trailing name columns in the `BOND` block are informational —
> connectivity is defined by the atom *indices* — so they may still show the old
> names after you edit the `ATOM` block. That is harmless, though you may tidy
> them if you want the file to read cleanly.

### Write the mainchain file

```bash
vi SME.mc
```

```
HEAD_NAME N
TAIL_NAME C
MAIN_CHAIN CA
OMIT_NAME CY
OMIT_NAME OY
OMIT_NAME CAY
OMIT_NAME HY1
OMIT_NAME HY2
OMIT_NAME HY3
OMIT_NAME NZ
OMIT_NAME HZ
OMIT_NAME CAT
OMIT_NAME HT1
OMIT_NAME HT2
OMIT_NAME HT3
PRE_HEAD_TYPE C
POST_TAIL_TYPE N
CHARGE 0.0
```

Reading it line by line:

- `HEAD_NAME` / `TAIL_NAME` — the atoms that bond to the previous and next
  residue.
- `MAIN_CHAIN` — the atoms on the path between them, *excluding* head and tail.
  For an α-amino acid that is just `CA`.
- `OMIT_NAME` — every cap atom. Twelve of them; `30 − 12 = 18` atoms will remain.
- `PRE_HEAD_TYPE` / `POST_TAIL_TYPE` — the **atom types** your neighbours will
  present: an ff14SB carbonyl `C` before you, an ff14SB amide `N` after you.
- `CHARGE` — the formal charge of the residue. **0.0 for SME; use your own.**

### Run it

```bash
prepgen -i SME.ac -o SME.prepi -m SME.mc -rn SME
```

Result:

```
   1  DUMM  DU    M    0  -1  -2     0.000      .0      .0    0.000000
   2  DUMM  DU    M    1   0  -1     1.449      .0      .0    0.000000
   3  DUMM  DU    M    2   1   0     1.523  111.21      .0    0.000000
   4  N     N     M    3   2   1     1.540  111.208 -180.000  -0.563045
   5  H     H     E    4   3   2     1.010   94.429   15.609   0.328965
   6  CA    CX    M    4   3   2     1.469  146.926  176.281   0.045315
   7  CB    c3    3    6   4   3     1.529  109.477  122.599  -0.060555
   ...
   9  S     s4    B    8   7   6     1.814  109.484 -179.976   0.273309
  10  OE    o     E    9   8   7     1.420  109.462  -65.041  -0.498249
   ...
LOOP
IMPROPER
 -M   CA   N    H
 CA   +M   C    O
DONE
STOP
```

### What `prepgen` did with the charges

The 18 surviving atoms carried, in the `.ac` file, a total charge of
**−0.054001**, not zero. AM1-BCC distributed charge over the
whole capped molecule, and the ACE and NME groups did not happen to come out
exactly neutral. Cut the caps off naively and you are left with a residue that
has a fractional charge.

`prepgen` fixes this. It takes the residual and spreads it over the remaining
atoms **in proportion to each atom's |q|**. You can verify this on the numbers
above: the residual is +0.054001, the sum of |q| over the 18 atoms is 3.924999,
and every single atom was scaled by

```
Δq / |q| = 0.054001 / 3.924999 = 0.013758
```

to within 2 × 10⁻⁶. The result sums to **+0.000001** — zero, to the precision of
the file.

> **Checkpoint.**
> ```bash
> python3 check_charge.py SME.prepi 0
> ```
> → `18 atoms, sum of charges = +0.000001` → `OK`

Remember this mechanism. In Route 2 nothing does it for you.

> **Also check the topology column** against the rules in 6.2, particularly at
> any atom your modification touches. Here `S` is `B` because it now has two
> children. Fix anything that looks wrong.

---

## 6.4 Route 2 — merging into a standard-residue template

This route builds the prepi by hand, starting from the template for the
*unmodified* parent residue. We present it for two reasons: 
it makes you read the format rather than trust a program, and
it is useful when the automatic route gives you something odd and you need to
understand why.

### Get the template

```bash
cp $AMBERHOME/dat/leap/prep/amino12.in .
vi amino12.in
```

Inside `vi`, type `/MET` and press Enter to jump to methionine (use whatever
three-letter code is the standard parent of *your* residue). Copy the block from
the title line down to `DONE`, and paste it into a new file:

```bash
vi SME.prepi
```

That block is reproduced in section 6.2 above.

### Get the computed charges

```bash
antechamber -i SME.ac -fi ac -o SME_ac.prepi -fo prepi -at gaff2 -rn SME
```

Note `-at gaff2`: without it `antechamber` falls back to GAFF, and the banner
will tell you so (`Info: The atom type is set to gaff`). We call the output
`SME_ac.prepi` to keep it distinct from the template-based `SME.prepi`.

This file describes the full capped molecule, ACE–SME–NME, with AM1-BCC charges.

### The merge

You now have two files, and you are going to combine them:

| from the **template** (`SME.prepi`) | from the **antechamber** file (`SME_ac.prepi`) |
|---|---|
| atom ordering and tree structure | partial charges |
| idealised internal coordinates | the extra atoms of your modification |
| ff14SB backbone atom types | GAFF2 sidechain atom types |
| the `IMPROPER` block | — |

Four operations, in order:

NOTE: the explanation below is for you to understand what shoud me modified in the file.
But if you decide to do it this way the advice is that you give both files to an AI, 
and then check the resulting file.

**1. Replace the charges.** For every atom present in both files, overwrite the
template's ff14SB charge with the AM1-BCC charge from `SME_ac.prepi`.

**2. Insert the new atoms.** SME has one atom that methionine does not: the
sulfoxide oxygen `OE`. Add a line for it, with:
- a sensible bond length, angle and dihedral (we use 1.485 Å, 106.5°, 60.0° —
  typical sulfoxide geometry);
- its charge, from `SME_ac.prepi`;
- **reference atom indices renumbered for this file.** The indices in
  `SME_ac.prepi` refer to the capped molecule and are meaningless here. `OE`
  refers to `S`, `CG`, `CB`, which in the template are atoms 14, 11, 8.

**3. Fix the topology letters.** `S` was `S` in methionine (one child, `CE`). It
now has two children and must become `B`.

**4. Fix the atom types and the names.** Sidechain atoms take their GAFF2 types
from `SME_ac.prepi`; backbone atoms keep the template's ff14SB types
(`N H CX HA C O`). And **rename the sulfur from the template's `SD` to the CCD's
`S`** — this is the trap flagged back in Part 2. A prepi that calls it `SD` will
not match any SME-containing PDB you ever load.

### Close the charge balance — the step that is easy to miss

`prepgen` absorbed the cap residual for you. Here, nothing has. If you simply
copy the 18 AM1-BCC charges across, your residue sums to **−0.054001**, and LEaP
will cheerfully build a system with a non-integer total charge. Under PME that
introduces a spurious uniform neutralising background, and your ion counts will
be wrong on top of it. Nothing will crash. You will not notice.

So do explicitly what `prepgen` does implicitly:

```bash
python3 check_charge.py SME.prepi 0
```

If the sum is off, the script prints the corrected charge for every atom, using
the same proportional-to-|q| redistribution `prepgen` uses. Paste those values in
and run it again until it says `OK`.

The finished file is provided as `SME_template_merge_balanced.prepi`:

```
    4 N    N    M    3   2   1     1.335   116.600   180.000   -0.563045
    5 H    H    E    4   3   2     1.010   119.800     0.000    0.328965
    6 CA   CX   M    4   3   2     1.449   121.900   180.000    0.045315
    7 HA   H1   E    6   4   3     1.090   109.500   300.000    0.088907
    8 CB   c3   3    6   4   3     1.525   111.100    60.000   -0.060555
   ...
   11 S    s4   B   11   8   6     1.810   110.000   180.000    0.273309
   12 OE   o    E   14  11   8     1.485   106.500    60.000   -0.498249
   ...
```

---

## 6.5 Comparing the two routes

| | Route 1 — `prepgen` | Route 2 — template merge |
|---|---|---|
| Effort | one small `.mc` file | line-by-line manual edit |
| Charge balance | **automatic** | **manual — must be checked** |
| Geometry | from your capped structure (AM1-optimised) | idealised template values |
| Atom ordering | determined by `prepgen` | inherited from the template |
| Easy to get wrong | atom naming in `OMIT_NAME` | everything |
| What it teaches you | what `prepgen` does for you | what a prepi *is* |

---

# Part 7 — Step 5: generating the missing parameters

`parmchk2` reads your library file, works out every bond, angle, torsion and
Lennard-Jones term the molecule will need, checks which ones the parameter
database already has, and writes the rest into a `.frcmod` by analogy.

## 7.1 Ligands — one pass

```bash
parmchk2 -i CER.mol2 -f mol2 -o CER.frcmod -s 2
```

> **`-s 2` is not optional.** `parmchk2`'s `-s` flag selects which parameter
> database to match against, and it **defaults to GAFF (`-s 1`)**. We assigned
> GAFF2 atom types in Part 3, so without `-s 2` the program draws its analogies
> and penalty scores from the wrong database. Pair them: `-at gaff` with `-s 1`,
> `-at gaff2` with `-s 2`. (Use the numeric form; the AMBER developers recommend
> `-s 2` over the spelled-out `-s gaff2`.)

## 7.2 Modified residues — two passes

Your residue has ff14SB types on the backbone and GAFF2 types on the sidechain.
That means it also has **cross terms that exist in neither database** — a `CX–c3`
bond, an `N–CX–c3` angle, and so on. A single `parmchk2` run against GAFF2 would
also override the perfectly good backbone parameters that `parm10.dat` already
contains.

The procedure, following the official AMBER tutorial on modified residues:

```bash
# Pass 1 - everything parm10.dat can supply, printed in full (-a Y)
parmchk2 -i SME.prepi -f prepi -o frcmod.SME -a Y \
         -p $AMBERHOME/dat/leap/parm/parm10.dat

# Keep only what parm10 genuinely provided
grep -v "ATTN" frcmod.SME > frcmod1.SME

# Pass 2 - let GAFF2 fill in the rest
parmchk2 -i SME.prepi -f prepi -o frcmod2.SME -s 2
```

Why `-a Y`: it prints *all* parameters, including perfect matches, not just the
missing ones. We need them written out explicitly so they can override the GAFF2
versions later.

Why `-p parm10.dat`: it is the main parameter database behind ff14SB. You can
confirm which file your force field uses by reading its `leaprc`. Omit the flag
and `parmchk2` matches against GAFF, which is not what you want for atoms you
deliberately gave AMBER types.

**The load order in LEaP then matters** — see Part 9. `frcmod2.SME` (GAFF2) is
loaded **first** and `frcmod1.SME` (parm10) **second**, so that wherever both
offer a parameter, the ff14SB value wins.

## 7.3 Reading the output

```
MASS
ns   14.010   0.530    same as n

BOND
hn-ns   403.20   1.013   same as hn- n, penalty score=  0.0
c3-ns   328.70   1.462   same as c3- n, penalty score=  0.0

IMPROPER
c3-ns-c -o   10.5  180.0  2.0   Using general improper torsional angle X- X- c- o, penalty score=  6.0
```

Every line tells you where the value came from and how far `parmchk2` had to
reach to get it. The **penalty score** is the distance of the analogy:

| penalty | interpretation |
|---|---|
| 0 | exact match in the database |
| up to ~10 | close analogue; normally fine |
| 10–50 | a reach; look at the line and decide |
| > 50 | poor analogy; verify before using |
| `ATTN: NEEDS REVISION` | no acceptable analogy found at all |

> **About `ATTN: NEEDS REVISION`.** `parmchk2` always writes *a* number — often
> zeros, sometimes a generic guess. The warning does not mean the file is
> incomplete; it means the value present is a placeholder with no justification.
> **Do not simulate with it.** Either supply the parameter from the literature,
> or compute it (a relaxed QM scan for a torsion, a Hessian-based fit for a bond
> or angle), or find a published parameter set for the fragment.
>
> ```bash
> grep -c ATTN *.frcmod     # should be 0 in the files you finally load
> ```

For our ceramide, the only line worth a second look is the improper
`c3-ns-c-o` with penalty 6 — a modest reach on a term that controls amide
planarity. We accept it here and verify by simulation in Part 8.

---

# Part 8 — Step 6: validating what you built

**1. The net charge is an integer.**
```bash
python3 check_charge.py SME.prepi 0
```
and later, in the LEaP log, `Total unperturbed charge: 0.000000`.

**2. No `ATTN` lines survive into the files you load.**
```bash
grep -c ATTN frcmod1.SME frcmod2.SME
```

**3. LEaP reports `Unit is OK`,** and nowhere says `Could not find ... parameter`.

**4. The geometry looks right.**
```bash
vmd SME.prmtop SME.ncrst
```
Bond lengths, planarity of conjugated groups, no atoms superimposed, correct
stereochemistry at every centre.

**5. For a residue: it builds into a peptide.** This is the test that a ligand
cannot fail and a residue often does, because it exercises the head/tail
definitions.

```bash
vi peptide.leap
```

```
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberprep SME.prepi
loadamberparams frcmod2.SME
loadamberparams frcmod1.SME

pep = sequence { ALA SME ALA }
charge pep
check pep
saveAmberParmNetCDF pep SMEpeptide.prmtop SMEpeptide.ncrst
savepdb pep SMEpeptide.pdb
quit
```

```bash
tleap -f peptide.leap
vmd SMEpeptide.prmtop SMEpeptide.ncrst
```

**6. A short minimisation runs without exploding.** Visual inspection does not
catch a wrong force constant or a nonsensical equilibrium value; minimisation
does, in under a minute.

```
Minimisation
 &cntrl
  imin=1, maxcyc=500, ncyc=200, ntb=0, igb=0, cut=999.0,
 /
```

```bash
sander -O -i min.in -p SMEpeptide.prmtop -c SMEpeptide.ncrst \
          -o min.out -r min.ncrst
```

Energies should decrease smoothly to a finite value. `NaN`, an energy of 10¹⁵,
or a structure that flies apart means a parameter is wrong — read `min.out` for
which term is large.

---

# Part 9 — Step 7: building your production system

Both branches converge here. You have your library file and your `.frcmod`; now
build the system you actually want to simulate.

### Prepare the PDB

- The residue or ligand must appear with the **same three-letter name** and the
  **same atom names** as in your library file.
- For a modified residue, simply put it in the sequence at the right position;
  LEaP will connect it using the head/tail you defined.
- Remove anything you do not want (crystallographic waters you have decided
  against, alternate conformations, hydrogens you want LEaP to rebuild).

### The LEaP script

```bash
vi system.leap
```

**Ligand version:**

```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberparams CER.frcmod
loadoff CER.lib

sys = loadPDB system.pdb
charge sys
check sys
```

**Modified-residue version:**

```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberprep SME.prepi
loadamberparams frcmod2.SME     # GAFF2 first ...
loadamberparams frcmod1.SME     # ... so ff14SB values override them

sys = loadPDB system.pdb
charge sys
check sys
```

> **The order of the two `loadamberparams` lines is deliberate.** LEaP lets later
> definitions override earlier ones. Loading the GAFF2 set first and the parm10
> set second means that wherever both supply a parameter, the ff14SB value is the
> one used — which is what we want for the backbone.

Then, common to both, solvation and ions:

```
solvateOct sys TIP3PBOX 10.0

# see the note below before filling these in
addIons sys Na+ 0
addIonsRand sys Na+ <N> Cl- <N>

charge sys
saveAmberParmNetCDF sys system.prmtop system.ncrst
savepdb sys system_built.pdb
quit
```

```bash
tleap -f system.leap
```

### About the ions — this needs two passes

You cannot write the ion counts before you know how many waters you have, and
you do not know that until after `solvateOct` has run. So:

**Pass 1.** Run the script with the ion lines commented out. Read the number of
water molecules from the LEaP output (or count `WAT` residues in
`system_built.pdb`).

**Pass 2.** Compute the counts and re-run:

```
number of ion pairs  =  n_waters × C / 55.5 # C = desired salt molarity, mol/L
cations to add       =  n_pairs − (system charge)/2
anions to add        =  n_pairs + (system charge)/2
```

55.5 mol/L is the molarity of pure water. The second and third lines only make
sense when `n_pairs` is larger than half the system charge; for a highly charged
system, neutralise first and add salt on top, which is what the two-line form
above does:

- `addIons sys Na+ 0` — the `0` means "add exactly enough to neutralise".
- `addIonsRand sys Na+ <N> Cl- <N>` — then add `<N>` pairs of bulk salt.
  `addIonsRand` places ions randomly and is much faster than `addIons` on large
  systems; `addIons` places them by electrostatic potential, which sounds better
  but tends to pile counterions onto charged surface patches.

> **Line to check:** the final `charge sys` must print `0.000000`. If it prints
> anything else, either your ion arithmetic is wrong or — more likely — your
> library file does not sum to an integer, and you should go back to Part 6.

### Final check

```bash
vmd system.prmtop system.ncrst
```

Your molecule should be present, correctly bonded, correctly placed, and
surrounded by water and ions.

**You are now ready to run.** And from this point onwards, the only two files you
need to carry around for this molecule are:

- **ligand:** `CER.lib` and `CER.frcmod`
- **modified residue:** `SME.prepi`, `frcmod1.SME` and `frcmod2.SME`

---

# Part 10 — Where this procedure stops being good enough

Everything above produces parameters that are reasonable for exploratory work.

**Charges.** AM1-BCC approximates HF/6-31G\* ESP charges rather than reproducing
them, is derived from a single conformation, and is not the model the surrounding
protein force field was built with. For published work on a system where the
modified group's electrostatics are the point — a phosphorylation site, a charged
ligand in a binding pocket — derive RESP charges from multiple conformations,
with the caps restrained to zero total charge so the residue stays compatible
with its neighbours. See the Part 3 box for how to do that without Gaussian.

**Torsions.** Any dihedral `parmchk2` filled in by analogy is a guess about a
conformational preference. If a torsion in your modification controls something
you care about — a sidechain rotamer, a ligand binding pose — scan it with QM and
fit it. This is the single most common reason published parameter sets exist at
all.

**Force-field generation.** This tutorial uses ff14SB with TIP3P because that
combination is ubiquitous and well documented. The current recommendation for new
protein work is **ff19SB with OPC water**. If you move to it, note a consequence
that is specific to what we have just done: ff19SB replaces the backbone torsion
corrections with **CMAP** terms, and CMAPs are defined only for the standard
amino acids. A residue you build by this procedure gets no CMAP, so its backbone
will effectively behave as it would under ff14SB while every residue around it
behaves as ff19SB. Whether that is acceptable depends on your question, but you
should know it is happening.

**Validation.** Parameters are not right because they were generated; they are
right because they reproduce something. If experimental or high-level QM data
exist for your molecule — a geometry, a vibrational frequency, a hydration free
energy, a rotamer population, an NMR observable — compare.

**Reproducibility.** The charges you obtained depend on the AmberTools version,
the input conformation and the charge method. Record all three. Keep the commands
in a shell script rather than typing them, and keep the input structure alongside
the outputs.
---

# Appendix A — Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `sqm` ends without `Calculation Completed` | Bad starting geometry, or wrong `-nc`. Minimise the structure first. If it is a convergence problem, loosen the criteria: `-ek "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-9, ndiis_attempts=700"` |
| `Fatal Error! Weird atomic valence` | `acdoctor` found an impossible valence: usually a missing or extra hydrogen, or a wrong protonation state. **Fix the structure.** Disabling the check with `-dr no` should be a last resort, not a first move |
| `Cannot properly run bondtype` | Bond perception failed, often from distorted geometry or a format that lost connectivity. Try a format that carries bond orders (`mol2`, `sdf`) |
| Charges do not sum to an integer | Expected for a capped fragment — see 6.3. In a ligand `.mol2` it means something is wrong with `-nc` |
| `Could not find bond parameter for: C - n` | GAFF type on a backbone atom. Part 4 |
| `Could not find angle parameter: CX - c3 - ...` | A genuine mixed-type cross term. Part 7.2 |
| LEaP: `Unknown residue: SME` | Library file not loaded, or the residue name in the PDB does not match the one in the library |
| LEaP: `Atom ... does not have a type` | Atom-name mismatch between your library file and your PDB. Compare the two name lists directly |
| LEaP: `FATAL: Atom .R<SME>.A<...> does not have a position` | Your PDB is missing an atom the library expects |
| Non-integer `Total unperturbed charge` | Part 6.4 |
| Residue at a chain terminus | The prepi you built is the **internal** form. An N- or C-terminal variant is a separate library entry, built the same way but capped on one side only |
| Minimisation gives `NaN` | A zero or absurd force constant, usually from an `ATTN` line you kept, or an atom overlap in the input |
| `parmchk2` output looks unrelated to your types | You forgot `-s 2` |

---

# Appendix B — Command summary

### Ligand

```bash
antechamber -i CER.pdb -fi pdb -o CER.mol2 -fo mol2 -at gaff2 -nc 0 -rn CER -c bcc -s 2 -pf y
tail sqm.out                                    # -> Calculation Completed
python3 check_charge.py CER.mol2 0              # -> OK
parmchk2  -i CER.mol2 -f mol2 -o CER.frcmod -s 2
grep -c ATTN CER.frcmod                         # -> 0
tleap -f lib.leap                               # -> CER.lib
```

### Modified residue

```bash
antechamber -i SME.pdb  -fi pdb  -o SME.mol2 -fo mol2 -at gaff2 -nc 0 -rn SME -c bcc -s 2
tail sqm.out                                    # -> Calculation Completed
antechamber -i SME.mol2 -fi mol2 -o SME.ac   -fo ac   -at gaff2
# edit SME.ac: backbone types -> N H CX HA C O ; rename cap atoms
prepgen   -i SME.ac -o SME.prepi -m SME.mc -rn SME
python3 check_charge.py SME.prepi 0             # -> OK
parmchk2  -i SME.prepi -f prepi -o frcmod.SME -a Y -p $AMBERHOME/dat/leap/parm/parm10.dat
grep -v ATTN frcmod.SME > frcmod1.SME
parmchk2  -i SME.prepi -f prepi -o frcmod2.SME -s 2
tleap -f peptide.leap                           # test build
tleap -f system.leap                            # production system
```

### The four checkpoints, in one place

```bash
tail sqm.out | grep "Calculation Completed"     # charges were actually computed
python3 check_charge.py <libfile> <formal charge>   # charge is an integer
grep -c ATTN *.frcmod                           # no unjustified parameters
grep -E "Unit is OK|Total unperturbed charge" leap.log
```

---

# Appendix C — Files accompanying this tutorial

| file | what it is |
|---|---|
| `check_charge.py` | sums the charges in a `.prepi`, `.mol2` or `.ac` and, if they are off, prints the corrected values |
| `SME_template_merge_balanced.prepi` | the finished Route 2 file, charge-balanced — the answer key for the 6.5 exercise |

---

## Sources and further reading

- AMBER Basic Tutorial 4b — parametrizing a ligand with antechamber.
- AMBER Basic Tutorial A26 (Tutorial 5) — simulating a modified residue; the
  source of the two-pass `parmchk2` procedure in Part 7.2.
- <https://ambermd.org/doc/prep.html> — the prep file format reference.
- Carlos Ramos-Guzmán, *AMBER custom residue parameterization*,
  <https://carlosramosg.com/amber-custom-residue-parameterization>.
- Jakalian, Jack & Bayly, *Fast, efficient generation of high-quality atomic
  charges. AM1-BCC model: II.* J. Comput. Chem. 2002 —
  <https://pubmed.ncbi.nlm.nih.gov/12395429/>.
- Maier et al., *ff14SB: Improving the accuracy of protein side chain and
  backbone parameters from ff99SB.* J. Chem. Theory Comput. 2015.
- Tian et al., *ff19SB: Amino-acid-specific protein backbone parameters trained
  against quantum mechanics energy surfaces in solution.* JCTC 2020.
- The Bryce Group AMBER Parameter Database, University of Manchester.
