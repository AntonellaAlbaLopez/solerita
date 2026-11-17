"""The tutorial, one step at a time.

    from solerita import Parametrization
    p = Parametrization("SME", charge=0)     # a PDB chemical-component code
    p.inspect()          # Part 2   look at the starting structure
    p.cap()              # Part 2.3 (residues) add ACE and NME caps
    p.charges()          # Part 3   AM1-BCC charges, GAFF2 atom types
    p.backbone()         # Part 4   (residues) ff14SB types on the backbone
    p.library()          # Part 5/6 the file LEaP loads (.lib or .prepi)
    p.parameters()       # Part 7   whatever the force field is missing
    p.test()             # Part 8   build it and minimise it
    p.finish()           # writes commands.sh with everything that was run

Every step prints why it exists, the exact command it runs, and what was
checked afterwards. Nothing is hidden: commands.sh at the end is a script you
can read, edit and run yourself.
"""
from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import textwrap

from . import formats as F
from . import protocol as P
from .cap import cap_residue, looks_precapped, rename_existing_caps
from .checks import CheckLog
from .inventory import take_inventory
from .structure import (adopt_reference_flags, ccd_reference, fetch_ccd,
                        read_ccd_cif, read_structure)

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
WIDTH = 76

LIGAND_STEPS = ["inspect", "charges", "parameters", "library", "test"]
RESIDUE_STEPS = ["inspect", "cap", "charges", "backbone", "library",
                 "parameters", "test"]


class StepFailed(RuntimeError):
    """A check failed badly enough that the next steps would be meaningless."""


def _wrap(text, indent="  "):
    return "\n".join(textwrap.fill(par, WIDTH, initial_indent=indent,
                                   subsequent_indent=indent)
                     for par in text.strip().split("\n\n"))


def amberhome():
    exe = shutil.which("antechamber")
    if not exe:
        raise RuntimeError(
            "AmberTools was not found: 'antechamber' is not on your PATH.\n"
            "  On Colab, run the installation cell at the top of the notebook."
            "\n  On your own computer, activate the environment first:\n"
            "      conda activate solerita")
    return os.environ.get("AMBERHOME") or os.path.dirname(
        os.path.dirname(os.path.realpath(exe)))


class Parametrization:
    """One molecule, taken through the tutorial step by step.

    source  a three-letter PDB chemical-component code ("SME", "CER") or the
            path to a structure file (.pdb, .mol2, .sdf, .cif)
    charge  the net charge of the molecule you want to parametrize. You must
            state it: it decides which chemical species the charges describe.
    resname residue name to use in the AMBER files (default: the code)
    kind    "ligand" or "residue"; detected from the PDB entry when possible
    offline use the copy of CER or SME shipped with solerita, not the PDB
    """

    def __init__(self, source, charge, resname=None, kind=None,
                 workdir=None, offline=False, pause=False):
        if charge is None:
            raise ValueError("Please state the net charge, e.g. charge=0. "
                             "solerita will not guess it for you.")
        self.source = str(source)
        self.charge = int(charge)
        self.is_file = os.path.exists(self.source)
        default_name = (os.path.splitext(os.path.basename(self.source))[0]
                        if self.is_file else self.source).upper()[:3]
        self.resname = (resname or default_name).upper()
        if not re.fullmatch(r"[A-Z0-9]{1,3}", self.resname):
            raise ValueError(f"residue name {self.resname!r}: use 1-3 "
                             f"letters or digits, e.g. resname='LIG'")
        self.kind = kind
        self.offline = offline
        self.pause = pause
        self.workdir = os.path.abspath(workdir or self.resname)
        os.makedirs(self.workdir, exist_ok=True)
        self.checks = CheckLog()
        self.commands = []           # (comment, shell text) for commands.sh
        self.done = []
        self.structure = None
        self.omit = []
        self.n_step = 0
        self._cur = "intake"
        self.frcmods = []

    # ------------------------------------------------------------ plumbing
    @property
    def steps(self):
        return RESIDUE_STEPS if self.kind == "residue" else LIGAND_STEPS

    def path(self, name):
        return os.path.join(self.workdir, name)

    def _start(self, step, part, title, why):
        if step != "inspect":
            if self.kind is None:
                raise StepFailed("Run p.inspect() first.")
            if step not in self.steps:
                raise StepFailed(f"'{step}' is not a step for a {self.kind}. "
                                 f"The steps are: {', '.join(self.steps)}")
            need = self.steps[:self.steps.index(step)]
            missing = [s for s in need if s not in self.done]
            if missing:
                raise StepFailed(f"Run p.{missing[0]}() first. The order "
                                 f"is: {' -> '.join(self.steps)}")
        self._cur = "intake" if step == "inspect" else step
        self.checks.checks = [c for c in self.checks.checks
                              if c.step != self._cur]
        self.commands = [c for c in self.commands if c[0] != self._cur]
        if step in self.done:
            # re-running a step invalidates everything after it
            self.done = self.done[:self.done.index(step)]
        self.n_step = (self.steps.index(step) + 1 if self.kind else 1)
        total = len(self.steps) if self.kind else "?"
        print("\n" + "=" * WIDTH)
        print(f"STEP {self.n_step}/{total}   {part}   {title}")
        print("=" * WIDTH)
        print(_wrap(why))
        print()
        if self.pause:
            input("  Press Enter to run this step... ")

    def _finish(self, step):
        if step not in self.done:
            self.done.append(step)
        nxt = (self.steps[self.steps.index(step) + 1]
               if self.kind and step in self.steps
               and self.steps.index(step) + 1 < len(self.steps) else None)
        print(f"\n  Step complete." + (f" Next: p.{nxt}()" if nxt
                                       else " Last step: p.finish()"))

    def _check(self, name, question, ok, message, guidance="",
               warn_only=False, evidence=""):
        c = self.checks.record(name, question, ok, message, evidence=evidence,
                               guidance=guidance, warn_only=warn_only,
                               step=self._cur)
        tag = {"PASS": "  ok  ", "WARN": " look ", "FAIL": " FAIL "}[c.result]
        print(f"  [{tag}] {question}")
        print(f"           {message}")
        if c.result != "PASS" and guidance:
            print(_wrap(("why it matters: " if c.result == "WARN"
                         else "what to do: ") + guidance, indent=" " * 11))
        if evidence and c.result != "PASS":
            for ln in evidence.strip().splitlines()[:4]:
                print(f"           | {ln}")
        return c.result

    def _note_command(self, comment, text):
        self.commands.append((self._cur, comment, text))

    def _run(self, argv, comment, log):
        """Run one AmberTools command, show it, and keep it for commands.sh."""
        pretty = " ".join(shlex.quote(a) for a in argv)
        print(f"  $ {pretty}")
        self._note_command(comment, pretty)
        env = dict(os.environ)
        env["AMBERHOME"] = amberhome()
        proc = subprocess.run(argv, cwd=self.workdir, env=env,
                              capture_output=True, text=True)
        with open(self.path(log), "w") as fh:
            fh.write(f"$ {pretty}\n\n{proc.stdout}\n{proc.stderr}")
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout).strip().splitlines()[-6:]
            self._check(f"{self._cur}:exit", f"did {argv[0]} finish normally?",
                        False, f"it stopped with exit code {proc.returncode}",
                        guidance=f"read {log} in {self.workdir} for the full "
                                 f"output", evidence="\n".join(tail))
            raise StepFailed(f"{argv[0]} failed; see {self.path(log)}")
        return proc.stdout + proc.stderr

    def _stop_if_failed(self, step):
        bad = [c for c in self.checks.checks
               if c.result == "FAIL" and c.step == self._cur]
        if bad:
            raise StepFailed(f"{len(bad)} check(s) failed in this step. Fix "
                             f"the problem above, then run the step again.")

    def _write(self, name, text, comment):
        with open(self.path(name), "w") as fh:
            fh.write(text)
        body = text if text.endswith("\n") else text + "\n"
        self._note_command(comment, f"cat > {name} <<'EOF'\n{body}EOF")

    # ---------------------------------------------------------------- Part 2
    def inspect(self):
        """Part 2 -- read the starting structure and say what is in it."""
        self._start(
            "inspect", "Part 2", "The starting structure",
            "Everything that follows describes exactly the structure you hand "
            "in: its atoms, its protons, its charge. Nothing downstream can "
            "repair a mistake made here, so the first job is to look.\n\n"
            "For a three-letter code, solerita downloads the entry from the "
            "PDB Chemical Component Dictionary (CCD), the reference "
            "dictionary of every small molecule and modified residue in the "
            "PDB. It gives the official atom names, the formal charge, and "
            "for a modified residue its parent amino acid.")
        st = None
        if self.is_file:
            print(f"  Reading {self.source}")
            st = read_structure(self.source)
        else:
            bundled = os.path.join(DATA, f"{self.source.upper()}.cif")
            url = f"https://files.rcsb.org/ligands/download/{self.source}.cif"
            if not self.offline:
                print(f"  Downloading {url}")
                try:
                    st = fetch_ccd(self.source,
                                   cache_dir=self.path(".ccd-cache"))
                except Exception as exc:
                    if not os.path.exists(bundled):
                        raise StepFailed(
                            f"could not download {self.source} ({exc}). "
                            f"Check the code and your internet connection.")
                    print(f"  Download failed ({type(exc).__name__}); using "
                          f"the copy shipped with solerita instead.")
            if st is None:
                if not os.path.exists(bundled):
                    raise StepFailed(f"no offline copy of {self.source}; "
                                     f"only CER and SME are bundled")
                print(f"  Using the bundled copy of {self.source}")
                st = read_ccd_cif(bundled)
            self._note_command("Part 2: the CCD entry this was built from",
                               f"wget {url}")
        reference = (ccd_reference(st, cache_dir=self.path(".ccd-cache"))
                     if self.is_file and not self.offline else None)
        adopt_reference_flags(st, reference)

        meta = dict(st.meta or {})
        if reference is not None:
            for k, v in (reference.meta or {}).items():
                meta.setdefault(k, v)
        ccd_type = str(meta.get("ccd_type") or "")
        if self.kind is None:
            self.kind = ("residue" if "PEPTIDE LINKING" in ccd_type.upper()
                         else "ligand")
            how = (f"the CCD calls it '{ccd_type}'" if ccd_type
                   else "no CCD information, so assumed; pass kind='residue' "
                        "if it is part of a protein chain")
            print(f"  Treating it as a {self.kind.upper()} ({how}).")
        print()
        inv = take_inventory(st, expect_charge=self.charge,
                             expect_kind=self.kind, reference=reference)
        print("\n".join("    " + ln for ln in inv.describe().splitlines()))
        parent = meta.get("parent_comp_id") or meta.get("mon_nstd_parent_comp_id")
        if parent and str(parent) not in ("?", "."):
            print(f"    parent residue  {parent}  (from the CCD)")
        leaving = [a.name for a in st.atoms if getattr(a, "leaving", False)]
        if leaving and self.kind == "residue":
            print(f"    leaving atoms   {' '.join(leaving)}  -- the CCD marks "
                  f"these as removed when the residue joins a chain")
        print()
        for f in inv.findings:
            if f.severity == "note":
                print(f"  [ note ] {f.message}")
                continue
            self._check(f"intake:{f.topic}",
                        f"{f.topic}: is the input usable as it is?",
                        f.severity != "error", f.message,
                        guidance=f.detail, warn_only=(f.severity == "warn"))
        self._stop_if_failed("intake")
        self.structure = st
        if self.kind == "ligand":
            out = f"{self.resname}.pdb"
            st.write_pdb(self.path(out), resname=self.resname)
            print(f"\n  Wrote {out} ({len(st.atoms)} atoms).")
        self._finish("inspect")

    # -------------------------------------------------------------- Part 2.3
    def cap(self):
        """Part 2.3 -- residues only: build ACE and NME caps."""
        self._start(
            "cap", "Part 2.3", "Capping the residue",
            "A residue inside a protein is bonded to its neighbours through "
            "peptide bonds. A free amino acid instead has an -NH2 and a "
            "-COOH, which put the wrong charges on the backbone atoms. So we "
            "compute the charges on ACE-X-NME: an acetyl cap on the N side "
            "and an N-methylamide cap on the C side, which mimic the "
            "neighbouring residues.\n\n"
            "solerita removes the atoms the CCD marks as leaving and builds "
            "both caps from standard peptide geometry, in an extended "
            "backbone. It also rebuilds the backbone H and O in the plane of "
            "the new peptide bonds.")
        st = self.structure
        if looks_precapped(st):
            print("  The structure already has caps; renaming them only.")
            rename_existing_caps(st, P.CAP_NAMES)
            scheme = P.CAP_NAMES
        else:
            st, rep = cap_residue(st, P.CAP_NAMES, head=P.HEAD_ATOM,
                                  tail=P.TAIL_ATOM)
            scheme = rep.scheme
            print(f"  removed  {' '.join(rep.deleted) or '-'}")
            print(f"  added    {' '.join(rep.added)}")
            if rep.rebuilt:
                print(f"  rebuilt  {' '.join(rep.rebuilt)}  (moved into the "
                      f"plane of the new peptide bond)")
            print(f"  backbone phi {rep.phi:+.0f}, psi {rep.psi:+.0f} "
                  f"(extended)\n")
            self._check("cap:geometry", "did the cap atoms land where "
                        "intended?", rep.worst_internal < 0.05,
                        f"largest deviation {rep.worst_internal:.2g} from "
                        f"standard peptide geometry")
            self._check("cap:clashes", "do the caps bump into the side "
                        "chain?", not rep.clashes,
                        "no close contacts" if not rep.clashes
                        else "close contacts: " + ", ".join(rep.clashes),
                        guidance="start from a different conformation, or "
                                 "cap the structure yourself", warn_only=True)
            for _, old, new in rep.substitutions:
                print(f"  [ note ] cap atom {old} renamed {new}: the residue "
                      f"already has an atom called {old}")
        self.omit = sorted(set(scheme["ACE"].values())
                           | set(scheme["NME"].values()))
        self._stop_if_failed("cap")
        out = f"{self.resname}.pdb"
        st.write_pdb(self.path(out), resname=self.resname)
        self.structure = st
        self._note_command(
            "Part 2.3: caps were added by solerita (Python), not by a shell "
            "command. The capped structure is " + out, f"# {out}")
        print(f"\n  Wrote {out}: ACE-{self.resname}-NME, {len(st.atoms)} atoms."
              f"\n  Open it in VMD or PyMOL and look at it before going on.")
        self._finish("cap")

    # ---------------------------------------------------------------- Part 3
    def charges(self):
        """Part 3 -- AM1-BCC charges and GAFF2 atom types (antechamber)."""
        self._start(
            "charges", "Part 3", "Charges and atom types",
            "antechamber does two things at once. It assigns every atom a "
            "GAFF2 atom type (a label that says which bonded parameters "
            "apply to it), and it computes AM1-BCC partial charges. For the "
            "charges it runs sqm, a semi-empirical quantum calculation, which "
            "first relaxes the geometry in vacuum. The -nc flag tells it the "
            "net charge; if that is wrong, every charge is wrong.\n\n"
            "This is the slow step: from a few seconds to a few minutes.")
        rn = self.resname
        out = f"{rn}.mol2"
        self._run(["antechamber", "-i", f"{rn}.pdb", "-fi", "pdb",
                   "-o", out, "-fo", "mol2", "-at", P.ATOM_TYPES,
                   "-nc", str(self.charge), "-rn", rn,
                   "-c", P.CHARGE_METHOD, "-s", "2"],
                  "Part 3: AM1-BCC charges and GAFF2 atom types",
                  "antechamber_charges.log")
        print()
        qm = None
        if os.path.exists(self.path("sqm.in")):
            m = re.search(r"qmcharge\s*=\s*(-?\d+)",
                          open(self.path("sqm.in")).read())
            qm = int(m.group(1)) if m else None
        self._check("charges:net", "was the calculation run at the charge "
                    "you asked for?", qm == self.charge,
                    f"you asked for {self.charge:+d}; sqm used "
                    f"{'unknown' if qm is None else f'{qm:+d}'}",
                    guidance="re-run with the correct charge")
        text = (open(self.path("sqm.out"), errors="replace").read()
                if os.path.exists(self.path("sqm.out")) else "")
        self._check("charges:sqm", "did the quantum calculation converge?",
                    "Calculation Completed" in text,
                    "sqm.out says 'Calculation Completed'"
                    if "Calculation Completed" in text
                    else "sqm.out does not report completion",
                    guidance="usually a bad starting structure (atoms too "
                             "close) or a wrong net charge")
        atoms, _ = F.parse_mol2(self.path(out))
        total = sum(a.charge for a in atoms)
        self._check("charges:sum", "do the charges add up to the net charge?",
                    abs(total - self.charge) <= max(P.CHARGE_TOLERANCE,
                                                    2e-3 * len(atoms)),
                    f"{len(atoms)} atoms, sum {total:+.6f}, "
                    f"net charge {self.charge:+d}")
        self._stop_if_failed("charges")
        types = sorted({a.atype for a in atoms})
        print(f"\n  Wrote {out}. GAFF2 types used: {' '.join(types)}")
        self._finish("charges")

    # ---------------------------------------------------------------- Part 4
    def backbone(self):
        """Part 4 -- residues only: ff14SB atom types on the backbone."""
        self._start(
            "backbone", "Part 4", "Backbone atom types",
            "GAFF2 typed every atom, including the backbone. But this "
            "residue has to bond to normal ff14SB residues, and the "
            "parameters for that junction only exist for ff14SB types. So "
            "the six backbone atoms get their ff14SB types back: "
            + " ".join(f"{k}->{v}" for k, v in P.BACKBONE_TYPES.items())
            + ". Only the labels change; the charges stay as they are.\n\n"
            "First the file is converted to .ac format, because that is "
            "what prepgen reads in the next step. The retyping is done on "
            "the .ac file: converting after retyping would undo it.")
        rn = self.resname
        ac = f"{rn}.ac"
        self._run(["antechamber", "-i", f"{rn}.mol2", "-fi", "mol2",
                   "-o", ac, "-fo", "ac", "-at", P.ATOM_TYPES, "-pf", "y"],
                  "Part 4: convert to .ac, the format prepgen reads",
                  "antechamber_ac.log")
        before = f"{rn}_gaff2_only.ac"
        shutil.copy(self.path(ac), self.path(before))
        atoms, _, _ = F.parse_ac(self.path(ac))
        names = {a.name for a in atoms}
        missing = [n for n in P.BACKBONE_TYPES if n not in names]
        if missing:
            self._check("backbone:names", "are the backbone atoms named "
                        "N, H, CA, HA, C, O?", False,
                        f"not found: {' '.join(missing)}",
                        guidance="rename them in your input structure")
            self._stop_if_failed("backbone")
        changed, _ = F.rewrite_ac_types(self.path(before), self.path(ac),
                                        dict(P.BACKBONE_TYPES))
        print()
        for name, old, new in changed:
            print(f"    {name:<3} {old:>3} -> {new}")
        after, _, _ = F.parse_ac(self.path(ac))
        got = {a.name: a.atype for a in after}
        self._check("backbone:types", "does the backbone now carry ff14SB "
                    "types?", all(got.get(n) == t
                                  for n, t in P.BACKBONE_TYPES.items()),
                    " ".join(f"{n}:{got.get(n)}" for n in P.BACKBONE_TYPES))
        q0 = sum(a.charge for a in atoms)
        q1 = sum(a.charge for a in after)
        self._check("backbone:charges", "were the charges left alone?",
                    abs(q0 - q1) < 1e-9, f"sum {q1:+.6f} before and after")
        self._note_command(
            f"Part 4: edit {ac} so the backbone types read "
            + " ".join(P.BACKBONE_TYPES.values())
            + f" (solerita did this; the unedited file is {before})",
            f"# {ac}")
        self._stop_if_failed("backbone")
        self._finish("backbone")

    # ------------------------------------------------------------- Part 5/6
    def library(self):
        """Part 5 (ligand .lib) or Part 6.3 (residue .prepi via prepgen)."""
        if self.kind == "ligand":
            return self._ligand_library()
        return self._residue_library()

    def _ligand_library(self):
        self._start(
            "library", "Part 5", "The ligand library (.lib)",
            "LEaP, the program that builds AMBER topologies, needs the "
            "molecule as a 'unit': atoms, types, charges and bonds. We load "
            "the mol2 and the frcmod into LEaP, let it check that every "
            "parameter is available, and save the unit as a .lib file. We "
            "also save a topology of the ligand alone to test in Part 8.")
        rn = self.resname
        script = (f"source {P.GAFF_LEAPRC}\n\n"
                  f"{rn} = loadmol2 {rn}.mol2\n"
                  f"loadamberparams {rn}.frcmod\n\n"
                  f"charge {rn}\ncheck {rn}\n\n"
                  f"saveoff {rn} {rn}.lib\n"
                  f"saveAmberParmNetCDF {rn} {rn}.prmtop {rn}.ncrst\n"
                  f"savepdb {rn} {rn}_check.pdb\nquit\n")
        self._write("lib.leap", script, "Part 5: the LEaP input")
        log = self._run(["tleap", "-f", "lib.leap"],
                        "Part 5: build and save the library", "tleap_lib.log")
        print()
        self._scan_leap(log, "library", rn)
        self._check("library:written", "was the .lib file written?",
                    os.path.exists(self.path(f"{rn}.lib")), f"{rn}.lib")
        print("\n  [ note ] LEaP says the residue lacks 'connect0/connect1'. "
              "That is expected for a ligand: it is not part of a chain.")
        self._stop_if_failed("library")
        self._finish("library")

    def _residue_library(self):
        self._start(
            "library", "Part 6.3", "The residue library (.prepi, prepgen)",
            "A residue needs more than a ligand: LEaP must know which atom "
            "bonds to the previous residue (the HEAD, N) and which to the "
            "next (the TAIL, C). prepgen writes a .prepi file with that "
            "information. We tell it, in a small 'mainchain' file, which "
            "atoms are head and tail, and which atoms to throw away: the "
            "caps, which were only there to get the charges right.\n\n"
            "prepgen then spreads the small leftover charge of the removed "
            "caps over the residue, so it sums exactly to the net charge.")
        rn = self.resname
        atoms, _, _ = F.parse_ac(self.path(f"{rn}.ac"))
        names = [a.name for a in atoms]
        omit = [n for n in self.omit if n in names]
        mc = [f"HEAD_NAME {P.HEAD_ATOM}", f"TAIL_NAME {P.TAIL_ATOM}"]
        mc += [f"MAIN_CHAIN {a}" for a in P.MAIN_CHAIN]
        mc += [f"OMIT_NAME {n}" for n in omit]
        mc += [f"PRE_HEAD_TYPE {P.PRE_HEAD_TYPE}",
               f"POST_TAIL_TYPE {P.POST_TAIL_TYPE}",
               f"CHARGE {float(self.charge):.1f}"]
        self._write(f"{rn}.mc", "\n".join(mc) + "\n",
                    "Part 6.3: the mainchain file")
        print("  " + f"{rn}.mc:  " + " | ".join(mc[:2]) + f" | omit "
              f"{len(omit)} cap atoms | CHARGE {float(self.charge):.1f}")
        out = f"{rn}.prepi"
        self._run(["prepgen", "-i", f"{rn}.ac", "-o", out, "-m", f"{rn}.mc",
                   "-rn", rn], "Part 6.3: build the prepi", "prepgen.log")
        print()
        p = F.parse_prepi(self.path(out))
        total = p.charge_sum()
        expected = len(names) - len(omit)
        self._check("library:atoms", "did prepgen keep exactly the residue "
                    "atoms?", len(p.real) == expected,
                    f"{len(p.real)} atoms kept "
                    f"({len(names)} minus {len(omit)} cap atoms)")
        self._check("library:charge", "does the residue add up to its net "
                    "charge?", abs(total - self.charge) <= P.PREPI_CHARGE_TOLERANCE,
                    f"sum of the charge column {total:+.6f}, net charge "
                    f"{self.charge:+d}",
                    guidance="the mainchain file is omitting the wrong atoms")
        issues, checked = F.check_prepi_tree(p)
        self._check("library:tree", "is the internal build tree consistent?",
                    not issues, f"{checked} atoms checked" + (
                        "" if not issues else "; " + "; ".join(
                            f"{n} ({t}): {m}" for n, t, m in issues)))
        lines = open(self.path(out)).read().splitlines(keepends=True)
        if len(lines) > 2 and lines[2].lower().startswith("this is a remark"):
            lines[2] = (f"{rn} net charge {total:+.6f} ({self.charge:+d}) | "
                        f"AM1-BCC/GAFF2 + ff14SB backbone | solerita\n")
            open(self.path(out), "w").writelines(lines)
        print(f"\n  [ note ] In {out}, the net charge is the sum of the last "
              f"column. The 0.0000 near the top is a bond cutoff, not the "
              f"charge.")
        self._stop_if_failed("library")
        self._finish("library")

    # ---------------------------------------------------------------- Part 7
    def parameters(self):
        """Part 7 -- the parameters the force field does not have (parmchk2)."""
        if self.kind == "ligand":
            why = ("GAFF2 covers most organic chemistry, but not every "
                   "combination of atom types. parmchk2 looks for every "
                   "bond, angle and torsion in the molecule, finds the ones "
                   "GAFF2 does not have, and proposes the closest analogue, "
                   "with a penalty score saying how far it had to reach. The "
                   "result is the .frcmod file.")
        else:
            why = ("This residue mixes two type families: ff14SB on the "
                   "backbone, GAFF2 on the side chain. So parmchk2 runs "
                   "twice. First against parm10.dat, the ff14SB parameter "
                   "file, keeping only what ff14SB really provides: lines "
                   "marked ATTN are placeholders, and so are the "
                   "zero-periodicity torsion lines next to them. Then "
                   "against GAFF2 for everything else.\n\n"
                   "In LEaP the GAFF2 file is loaded first and the parm10 "
                   "file second, so where both have a value the ff14SB one "
                   "wins.")
        self._start("parameters", "Part 7", "Missing parameters", why)
        rn = self.resname
        if self.kind == "ligand":
            files = [f"{rn}.frcmod"]
            self._run(["parmchk2", "-i", f"{rn}.mol2", "-f", "mol2",
                       "-o", files[0], "-s", str(P.PARMCHK2_GAFF)],
                      "Part 7.1: missing parameters, GAFF2", "parmchk2.log")
        else:
            raw = f"frcmod_parm10.{rn}"
            parm10 = os.path.join(amberhome(), "dat", "leap", "parm",
                                  "parm10.dat")
            self._run(["parmchk2", "-i", f"{rn}.prepi", "-f", "prepi",
                       "-o", raw, "-a", "Y", "-p", parm10],
                      "Part 7.2, pass 1: everything parm10.dat can supply",
                      "parmchk2_parm10.log")
            n_attn, n_null = F.strip_attn(self.path(raw),
                                          self.path(f"frcmod1.{rn}"))
            self._note_command(
                "Part 7.2: keep only what parm10 really provides -- drop the "
                "ATTN placeholders and the zero-periodicity torsion lines "
                "parmchk2 writes next to them",
                f"grep -v ATTN {raw} | awk '/^DIHE/{{d=1}} /^IMPROPER/{{d=0}} "
                f"d && length($0)>11 {{split(substr($0,12),f,\" \"); "
                f"if (f[4]!=\"\" && f[4]+0==0) next}} {{print}}' "
                f"> frcmod1.{rn}")
            print(f"  removed {n_attn} ATTN lines and {n_null} "
                  f"zero-periodicity torsion lines -> frcmod1.{rn}")
            self._run(["parmchk2", "-i", f"{rn}.prepi", "-f", "prepi",
                       "-o", f"frcmod2.{rn}", "-s", str(P.PARMCHK2_GAFF)],
                      "Part 7.2, pass 2: GAFF2 fills in the rest",
                      "parmchk2_gaff2.log")
            files = [f"frcmod2.{rn}", f"frcmod1.{rn}"]     # LEaP load order
        print()
        attn, bands, poor = [], {}, []
        for f in files:
            terms, b = F.parse_frcmod(self.path(f), P.PENALTY_OK,
                                      P.PENALTY_REACH)
            attn += [t for t in terms if t.attn]
            for k, v in b.items():
                bands[k] = bands.get(k, 0) + len(v)
            poor += b.get("poor", [])
        self._check("parameters:attn", "is any parameter still a "
                    "placeholder (ATTN)?", not attn,
                    "none" if not attn else
                    f"{len(attn)}: " + ", ".join(f"{t.section} {t.key}"
                                                 for t in attn[:5]),
                    guidance="a placeholder has no physical value; find the "
                             "parameter in the literature or compute it")
        summary = ", ".join(f"{v} {k}" for k, v in sorted(bands.items()))
        self._check("parameters:penalty", "how close are the analogies "
                    "parmchk2 used?", not poor,
                    summary or "no new parameters were needed",
                    guidance="a penalty above 50 is a distant analogy: look "
                             "at the term and decide whether it matters",
                    warn_only=True,
                    evidence="\n".join(f"{t.section} {t.key} penalty "
                                       f"{t.penalty:g}" for t in poor[:4]))
        self.frcmods = files
        self._stop_if_failed("parameters")
        print(f"\n  Files to load in LEaP, in this order: {' then '.join(files)}")
        self._finish("parameters")

    # ---------------------------------------------------------------- Part 8
    def test(self):
        """Part 8 -- does it build, and does a short minimisation behave?"""
        why = ("Getting files is not the same as getting correct files. Two "
               "tests catch most mistakes. LEaP must build the molecule "
               "without complaining and agree on its total charge. And a "
               "short energy minimisation must give a finite energy that "
               "goes down: a wrong or missing parameter usually shows up "
               "here as a huge energy or NaN.")
        if self.kind == "residue":
            why += ("\n\nFor a residue the real test is whether it joins a "
                    "chain, so we build the tripeptide ALA-X-ALA.")
        self._start("test", "Part 8", "Validation", why)
        rn = self.resname
        if self.kind == "residue":
            loads = "\n".join(f"loadamberparams {f}" for f in self.frcmods)
            script = (f"source {P.PROTEIN_FF}\nsource {P.GAFF_LEAPRC}\n"
                      f"loadamberprep {rn}.prepi\n{loads}\n\n"
                      f"pep = sequence {{ ALA {rn} ALA }}\n"
                      f"charge pep\ncheck pep\n"
                      f"saveAmberParmNetCDF pep {rn}peptide.prmtop "
                      f"{rn}peptide.ncrst\n"
                      f"savepdb pep {rn}peptide.pdb\nquit\n")
            self._write("peptide.leap", script, "Part 8: the test peptide")
            log = self._run(["tleap", "-f", "peptide.leap"],
                            "Part 8: build ALA-X-ALA", "tleap_peptide.log")
            print()
            self._scan_leap(log, "test", f"ALA-{rn}-ALA")
            top, crd = f"{rn}peptide.prmtop", f"{rn}peptide.ncrst"
            self._check("test:built", "did the tripeptide build?",
                        os.path.exists(self.path(top)),
                        f"ALA-{rn}-ALA -> {top}",
                        guidance="LEaP could not join the residue to its "
                                 "neighbours: check HEAD/TAIL in the .mc file")
            self._stop_if_failed("test")
        else:
            top, crd = f"{rn}.prmtop", f"{rn}.ncrst"
        self._write("min.in", P.MIN_INPUT, "Part 8: minimisation input")
        self._run(["sander", "-O", "-i", "min.in", "-p", top, "-c", crd,
                   "-o", "min.out", "-r", "min.ncrst"],
                  "Part 8: short minimisation", "sander.log")
        text = open(self.path("min.out"), errors="replace").read()
        e = [float(m.group(2)) for m in
             re.finditer(r"^\s*(\d+)\s+(-?[\d.]+E[+-]\d+)\s", text, re.M)]
        finite = bool(e) and all(abs(x) < 1e10 for x in e) \
            and "NaN" not in text
        self._check("test:energy", "is the energy finite, and does it go "
                    "down?", finite and e[-1] <= e[0],
                    (f"{e[0]:.2f} -> {e[-1]:.2f} kcal/mol over "
                     f"{len(e)} printed steps") if e
                    else "no energies found in min.out",
                    guidance="NaN or a huge value usually means a placeholder "
                             "parameter or overlapping atoms")
        self._stop_if_failed("test")
        self._finish("test")

    def _scan_leap(self, log, step, label):
        bad = [ln.strip() for ln in log.splitlines()
               if "FATAL" in ln or "Could not find" in ln
               or "does not have a type" in ln]
        self._check(f"{step}:leap", "did LEaP find every parameter?",
                    not bad, "no missing parameters" if not bad
                    else f"{len(bad)} problem(s)", evidence="\n".join(bad),
                    guidance="a parameter is missing from the frcmod files")
        m = re.search(r"Total unperturbed charge:\s*(-?[\d.]+)", log)
        q = float(m.group(1)) if m else None
        self._check(f"{step}:charge", "does LEaP agree on the total charge?",
                    q is not None and abs(q - self.charge) < 1e-3,
                    f"LEaP: {q if q is not None else 'not reported'}, "
                    f"expected {self.charge:+d}")
        self._check(f"{step}:check", f"does LEaP's own check of {label} "
                    f"pass?", "Unit is OK" in log,
                    "'Unit is OK'" if "Unit is OK" in log
                    else "LEaP did not say 'Unit is OK'")

    # ---------------------------------------------------------------- finish
    def run_all(self):
        """Every step in order, then finish(). Used by the command line."""
        self.inspect()
        for step in self.steps[1:]:
            getattr(self, step)()
        return self.finish()

    def finish(self):
        """Write commands.sh and print a summary of what you now have."""
        rn = self.resname
        L = ["#!/usr/bin/env bash",
             f"# Every command solerita ran to parametrize {rn} "
             f"({self.kind}, net charge {self.charge:+d}).",
             "# Read it next to the tutorial to see what each step does.",
             "# Two steps are not shell commands and were done by solerita",
             "# in Python: building the caps (Part 2.3) and retyping the",
             "# backbone in the .ac file (Part 4). To run this script yourself",
             "# you must do those two by hand, as the tutorial describes,",
             "# at the points marked below.",
             "set -e", ""]
        order = {s: i for i, s in enumerate(["intake"] + self.steps[1:])}
        for _, comment, text in sorted(self.commands,
                                       key=lambda c: order.get(c[0], 99)):
            L += [f"# {comment}", text, ""]
        with open(self.path("commands.sh"), "w") as fh:
            fh.write("\n".join(L))
        os.chmod(self.path("commands.sh"), 0o755)
        n = {r: sum(c.result == r for c in self.checks.checks)
             for r in ("PASS", "WARN", "FAIL")}
        keep = ([f"{rn}.lib", f"{rn}.frcmod"] if self.kind == "ligand"
                else [f"{rn}.prepi", f"frcmod2.{rn}", f"frcmod1.{rn}"])
        print("\n" + "=" * WIDTH)
        print(f"DONE   {rn}   {n['PASS']} checks passed, {n['WARN']} to look "
              f"at, {n['FAIL']} failed")
        print("=" * WIDTH)
        print(f"  Folder: {self.workdir}")
        print(f"  The files you need to use {rn} in a simulation:")
        for f in keep:
            print(f"    {f}")
        print("  Every command, in order: commands.sh")
        if n["WARN"]:
            print("\n  Look at these before using the parameters:")
            for c in self.checks.checks:
                if c.result == "WARN":
                    print(f"    - {c.question}  {c.message}")
        return self
