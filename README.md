# solerita

Parametrize a **ligand** or a **modified amino acid** for AMBER, one step at a
time, following the class tutorial. At every step solerita:

1. tells you **why** the step exists,
2. shows the **exact command** it runs,
3. runs it, and
4. **checks** the result and says, in plain words, what was checked.

Nothing is hidden. At the end, `commands.sh` contains every command in
order: the same procedure as the tutorial, as a script you can read and run.

---

## Option A: Google Colab (nothing to install on your computer)

1. Open `notebooks/solerita_class.ipynb` in Colab
   (File → Upload notebook, or open it from Google Drive).
2. Run the first cell. It installs AmberTools (1–2 minutes) and asks you to
   upload `solerita.zip`.
3. Run the cells from top to bottom, reading each output.

## Option B: your own computer (Linux, macOS, or Windows with WSL)

You need conda or mamba. Then, from the `solerita` folder:

```bash
conda env create -f environment.yml     # AmberTools + solerita, once
conda activate solerita                 # every time you open a terminal
solerita --check                        # should end with "Everything is in place."
```

---

## Using it

### From a terminal

```bash
solerita CER --charge 0      # ceramide, a ligand
solerita SME --charge 0      # methionine sulfoxide, a modified residue
```

It pauses before each step: read, then press Enter. A folder named after the
molecule (`CER/`, `SME/`) holds everything it makes.

Other options:

| | |
|---|---|
| `solerita my_ligand.pdb --charge -1 --resname LIG` | use your own structure file |
| `--offline` | use the copies of CER and SME shipped with solerita, if there is no internet |
| `--no-pause` | run all steps without stopping |
| `--kind residue` | say it is a residue if solerita cannot tell |

### From Python or Jupyter

```python
from solerita import Parametrization

p = Parametrization("SME", charge=0)
p.inspect()       # Part 2   look at the starting structure
p.cap()           # Part 2.3 add ACE and NME caps          (residues only)
p.charges()       # Part 3   AM1-BCC charges, GAFF2 types
p.backbone()      # Part 4   ff14SB types on the backbone   (residues only)
p.library()       # Part 6   the .prepi (Part 5 .lib for a ligand)
p.parameters()    # Part 7   missing parameters (frcmod)
p.test()          # Part 8   build it and minimise it
p.finish()        # summary, and commands.sh
```

For a **ligand** the order is `inspect, charges, parameters, library, test`:
LEaP needs the parameters before it can write the `.lib`. If you call a step
too early, solerita tells you which one comes first.

---

## The molecule code

Three-letter codes come from the **PDB Chemical Component Dictionary (CCD)**,
the PDB's dictionary of every small molecule and modified residue it
contains. Search a code at <https://www.rcsb.org/> to see its page. The CCD
entry gives the official atom names, and for a modified residue its parent
amino acid. That is how solerita knows whether it is a ligand or a residue.

## The net charge is yours to give

solerita will not guess it. The charges describe exactly the structure you
give it, with exactly its hydrogens. Saying `--charge -1` does **not** remove
a proton: if you want the deprotonated form, remove the hydrogen from the
structure yourself and give solerita that file. Step 2 prints the hydrogens it
found, so check them against the charge you gave.

## What the checks mean

| mark | meaning |
|---|---|
| `[  ok  ]` | checked, and fine |
| `[ look ]` | not an error, but look at it before using the parameters; it says why |
| `[ FAIL ]` | something is wrong; it says what to do. The run stops there. |
| `[ note ]` | information only |

After fixing the cause of a `FAIL`, run that step again, then the steps after it.

## The files you end up with

| | ligand | modified residue |
|---|---|---|
| the molecule for LEaP | `CER.lib` | `SME.prepi` |
| extra parameters | `CER.frcmod` | `frcmod2.SME` then `frcmod1.SME`, in this order |
| everything that was run | `commands.sh` | `commands.sh` |
| test topology | `CER.prmtop` | `SMEpeptide.prmtop` (ALA-SME-ALA) |

Every program's full output is kept in a `.log` file named after the step.

## Where the choices are written down

`solerita/protocol.py` lists every scientific choice (charge method, atom
types, backbone types, penalty thresholds, minimisation settings), each with
a short explanation. It is the only file that decides anything; the rest
only carries out those choices and checks the results.

## Two things solerita does differently from typing the tutorial by hand

* **Caps are built by solerita** from standard peptide geometry, instead of
  in PyMOL. You can still give it a structure you capped yourself.
* **The alpha hydrogen is typed `H1`.** That is the ff14SB type for it in
  every standard residue. The tutorial table that says `HA` is being
  corrected (`HA` is the aromatic-hydrogen type).
