"""Say out loud what was read, and flag what a human must decide.

Nothing here changes the molecule. It produces the evidence the interview is
conducted against, and the findings that become checkpoints. A finding has a
severity:

  error  the workflow does not apply, or the input cannot be used as given
  warn   proceed, but a human has to have seen this
  note   recorded for the report, no action
"""
from __future__ import annotations

from dataclasses import dataclass, field

# GAFF2 covers organic main-group chemistry. Outside this set the protocol
# does not apply and no amount of parmchk2 will rescue it.
GAFF2_ELEMENTS = {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"}

MAX_DEGREE = {"H": 1, "C": 4, "N": 4, "O": 2, "S": 6, "P": 5,
              "F": 1, "Cl": 1, "Br": 1, "I": 1}

CCD_RESIDUE_TYPES = ("PEPTIDE LINKING", "PEPTIDE-LIKE", "DNA LINKING",
                     "RNA LINKING", "SACCHARIDE")


@dataclass
class Finding:
    severity: str
    topic: str
    message: str
    detail: str = ""

    def __str__(self):
        s = f"[{self.severity:<5}] {self.topic}: {self.message}"
        return s + (f"\n          {self.detail}" if self.detail else "")


@dataclass
class Inventory:
    structure: object = None
    findings: list = field(default_factory=list)
    facts: dict = field(default_factory=dict)

    def add(self, severity, topic, message, detail=""):
        self.findings.append(Finding(severity, topic, message, detail))

    @property
    def errors(self):
        return [f for f in self.findings if f.severity == "error"]

    @property
    def warnings(self):
        return [f for f in self.findings if f.severity == "warn"]

    def describe(self):
        st, F = self.structure, self.facts
        L = [f"input: {st.source_path or '(in memory)'}",
             f"  name / title    {st.name}"
             + (f"  -- {st.title}" if st.title else ""),
             f"  atoms           {len(st.atoms)}  "
             f"({F['n_heavy']} heavy, {F['n_h']} H)",
             f"  formula         {st.formula()}   "
             f"mass {st.mass():.2f}",
             f"  bonds           {len(st.bonds)}"
             + ("  (PERCEIVED from distances, not declared)"
                if st.bonds_perceived else "  (declared by the file)"),
             f"  formal charge   {st.formal_charge():+d}  (sum of per-atom "
             f"formal charges in the file)"]
        if F.get("ccd_type"):
            L.append(f"  CCD type        {F['ccd_type']}"
                     f"  -> looks like a {F['inferred_kind']}")
        if F.get("parent_comp_id"):
            L.append(f"  parent residue  {F['parent_comp_id']}")
        if F.get("leaving"):
            L.append(f"  leaving atoms   {' '.join(F['leaving'])}"
                     f"   (CCD flags these as removed on polymerisation)")
        if F.get("stereocentres"):
            L.append("  stereocentres   "
                     + ", ".join(f"{n}:{c}" for n, c in F["stereocentres"]))
        if F.get("fragments", 1) > 1:
            L.append(f"  fragments       {F['fragments']}  (the file contains "
                     f"more than one disconnected molecule)")
        L.append(f"  backbone found  "
                 + ("yes: " + " ".join(F["backbone_atoms"])
                    if F.get("backbone_atoms") else "no alpha-amino-acid "
                    "backbone (N, CA, C, O) detected"))
        return "\n".join(L)


def _fragments(st):
    n = len(st.atoms)
    seen, count = set(), 0
    adj = {i: [] for i in range(n)}
    for a, b, _ in st.bonds:
        adj[a].append(b)
        adj[b].append(a)
    for i in range(n):
        if i in seen:
            continue
        count += 1
        stack = [i]
        while stack:
            k = stack.pop()
            if k in seen:
                continue
            seen.add(k)
            stack.extend(adj[k])
    return count


def take_inventory(st, expect_charge=None, expect_kind=None, reference=None):
    inv = Inventory(structure=st)
    F = inv.facts
    F["n_h"] = st.n_hydrogens()
    F["n_heavy"] = len(st.atoms) - F["n_h"]
    F["elements"] = sorted({a.element for a in st.atoms})
    F["leaving"] = [a.name for a in st.atoms if a.leaving]
    F["stereocentres"] = [(a.name, a.stereo) for a in st.atoms
                          if a.stereo not in ("", "N")]
    F["ccd_type"] = st.meta.get("ccd_type", "")
    F["parent_comp_id"] = st.meta.get("parent_comp_id", "")
    F["fragments"] = _fragments(st) if st.bonds else len(st.atoms)

    ccd_type = F["ccd_type"].upper()
    if ccd_type:
        F["inferred_kind"] = ("residue"
                              if any(t in ccd_type for t in CCD_RESIDUE_TYPES)
                              else "ligand")
    else:
        F["inferred_kind"] = ""

    bb = [n for n in ("N", "CA", "C", "O") if st.index(n) >= 0]
    F["backbone_atoms"] = bb if len(bb) == 4 else []

    # --- elements outside the force field --------------------------------
    outside = [e for e in F["elements"] if e not in GAFF2_ELEMENTS]
    if outside:
        inv.add("error", "elements",
                f"{', '.join(outside)} not covered by GAFF2",
                "This workflow assigns GAFF2 types and AM1-BCC charges. "
                "Neither applies to metals or open-shell centres; a bonded "
                "or non-bonded metal model is a different procedure.")

    # --- hydrogens ---------------------------------------------------------
    if F["n_h"] == 0:
        inv.add("error", "hydrogens",
                "no hydrogens in the input",
                "AM1-BCC charges are computed for the structure you hand in. "
                "Add all hydrogens explicitly, in the protonation state you "
                "intend to simulate, before continuing.")
    elif F["n_h"] < 0.5 * max(1, sum(1 for a in st.atoms
                                     if a.element in ("C", "N", "O"))):
        inv.add("warn", "hydrogens",
                f"only {F['n_h']} hydrogens for {F['n_heavy']} heavy atoms",
                "This is low for an organic molecule. Check that none are "
                "missing, especially on heteroatoms.")

    # --- names ------------------------------------------------------------
    dup = st.duplicate_names()
    if dup:
        inv.add("error", "atom names",
                f"duplicate atom names: {', '.join(dup)}",
                "prepgen and LEaP match atoms by name; duplicates are "
                "ambiguous and will be rejected or silently mismatched.")
    long_names = [a.name for a in st.atoms if len(a.name) > 4]
    if long_names:
        inv.add("error", "atom names",
                f"atom names longer than 4 characters: "
                f"{', '.join(long_names[:5])}",
                "The PDB format allows four columns for an atom name.")

    # --- geometry ----------------------------------------------------------
    clash = st.close_contacts(0.7)
    if clash:
        inv.add("error", "geometry",
                f"{len(clash)} atom pairs closer than 0.7 A",
                "e.g. " + ", ".join(f"{a}-{b}" for a, b in clash[:5])
                + ". sqm will not converge on this.")
    if F["fragments"] > 1:
        inv.add("warn", "connectivity",
                f"the file contains {F['fragments']} disconnected fragments",
                "A parametrization input should contain exactly one "
                "molecule. Counter-ions, waters or a second copy must go.")
    if st.bonds_perceived:
        inv.add("note", "connectivity",
                "connectivity was inferred from interatomic distances",
                "The input format carried no bond block. A distorted "
                "geometry can therefore produce wrong bonds; a format that "
                "declares them (mol2, sdf, CCD cif) is safer.")

    # --- valence ------------------------------------------------------------
    deg = {i: 0 for i in range(len(st.atoms))}
    for a, b, _ in st.bonds:
        deg[a] += 1
        deg[b] += 1
    bad = []
    for i, a in enumerate(st.atoms):
        cap = MAX_DEGREE.get(a.element)
        if cap is not None and deg[i] > cap:
            bad.append(f"{a.name}({a.element}) has {deg[i]} bonds")
        if deg[i] == 0 and len(st.atoms) > 1:
            bad.append(f"{a.name} is bonded to nothing")
    if bad:
        inv.add("error", "valence", "impossible valences",
                "; ".join(bad[:6])
                + ". antechamber's acdoctor will refuse this input.")

    # --- charge --------------------------------------------------------------
    ccd_q = st.meta.get("ccd_formal_charge", "")
    if ccd_q not in ("", None):
        try:
            F["ccd_formal_charge"] = int(float(ccd_q))
        except ValueError:
            pass

    # --- what this file was edited from, if anything ---------------------
    implied, implied_basis = None, ""
    if reference is not None:
        ref_el = {a.name: a.element for a in reference.atoms}
        cur_el = {a.name: a.element for a in st.atoms}
        removed = [n for n in ref_el if n not in cur_el]
        added = [n for n in cur_el if n not in ref_el]
        rem_h = [n for n in removed if ref_el[n] == "H"]
        rem_hv = [n for n in removed if ref_el[n] != "H"]
        add_h = [n for n in added if cur_el[n] == "H"]
        add_hv = [n for n in added if cur_el[n] != "H"]
        try:
            ref_q = int(float(reference.meta.get("ccd_formal_charge") or 0))
        except (TypeError, ValueError):
            ref_q = 0
        F["reference_code"] = reference.name
        F["reference_charge"] = ref_q
        # An edited PDB has lost the component metadata. The entry it came
        # from still has it, so the parent residue and the component type
        # need not be asked again just because the file was edited.
        for key in ("parent_comp_id", "ccd_type", "one_letter"):
            if not F.get(key) and reference.meta.get(key):
                F[key] = reference.meta[key]
        ccd_type = F["ccd_type"].upper()
        if not F.get("inferred_kind") and ccd_type:
            F["inferred_kind"] = ("residue"
                                  if any(t in ccd_type
                                         for t in CCD_RESIDUE_TYPES)
                                  else "ligand")
        F["removed_vs_reference"] = sorted(removed)
        F["added_vs_reference"] = sorted(added)
        if removed or added:
            inv.add("note", "edited from CCD",
                    f"compared with CCD entry {reference.name}: "
                    + (f"removed {' '.join(sorted(removed))}" if removed
                       else "nothing removed")
                    + (f"; added {' '.join(sorted(added))}" if added else ""))
        if (rem_h or add_h) and not rem_hv and not add_hv:
            # Only protons changed, so the charge follows arithmetically.
            implied = ref_q - len(rem_h) + len(add_h)
            implied_basis = (f"CCD entry {reference.name} is {ref_q:+d}; "
                             + (f"{len(rem_h)} H removed "
                                f"({' '.join(sorted(rem_h))})"
                                if rem_h else "")
                             + ("; " if rem_h and add_h else "")
                             + (f"{len(add_h)} H added "
                                f"({' '.join(sorted(add_h))})"
                                if add_h else ""))
            F["implied_charge"] = implied
            F["implied_basis"] = implied_basis
        elif rem_hv or add_hv:
            inv.add("note", "net charge",
                    "heavy atoms differ from the CCD entry, so the formal "
                    "charge cannot be worked out by proton arithmetic; the "
                    "charge you declare is taken on trust")

    # --- declared charge versus whatever the input can actually attest ---
    if expect_charge is not None:
        if implied is not None:
            if implied != expect_charge:
                inv.add("error", "net charge",
                        f"you declared {expect_charge:+d}, but removing and "
                        f"adding those protons gives {implied:+d}",
                        f"{implied_basis}. Either declare {implied:+d}, or "
                        f"go back and edit the structure so that it really "
                        f"is the species you mean -- declaring a charge "
                        f"changes nothing in the file.")
            else:
                inv.add("note", "net charge",
                        f"declared {expect_charge:+d}, consistent with the "
                        f"protons present",
                        implied_basis)
        elif "ccd_formal_charge" in F:
            if F["ccd_formal_charge"] != expect_charge:
                inv.add("warn", "net charge",
                        f"you declared {expect_charge:+d} but the CCD entry "
                        f"is {F['ccd_formal_charge']:+d}",
                        "The CCD entry is very often the neutral species "
                        "(for a phosphorylated residue it is the fully "
                        "protonated form). If you mean the ionised state you "
                        "must remove the protons yourself; declaring the "
                        "charge alone does not change the molecule.")
        elif any(a.formal_charge for a in st.atoms):
            if st.formal_charge() != expect_charge:
                inv.add("warn", "net charge",
                        f"you declared {expect_charge:+d} but the per-atom "
                        f"formal charges in the file sum to "
                        f"{st.formal_charge():+d}")
        else:
            inv.add("note", "net charge",
                    f"declared {expect_charge:+d}; nothing in this file can "
                    f"confirm it",
                    "A PDB carries no formal charges, and no CCD entry was "
                    "found to diff against, so the declared charge is taken "
                    "on trust. It is still checked downstream: antechamber "
                    "is given this charge, and the final parameters are "
                    "verified to sum to it.")

    # --- kind ----------------------------------------------------------------
    if expect_kind and F["inferred_kind"] and expect_kind != F["inferred_kind"]:
        inv.add("warn", "kind",
                f"you said {expect_kind}, the CCD type "
                f"{F['ccd_type']!r} suggests {F['inferred_kind']}")

    # --- stereochemistry -------------------------------------------------------
    if F["stereocentres"]:
        inv.add("note", "stereochemistry",
                f"{len(F['stereocentres'])} stereocentre(s): "
                + ", ".join(f"{n} ({c})" for n, c in F["stereocentres"]),
                "Each is a different molecule with different parameters. "
                "Confirm the structure is the isomer you mean and record it.")
    return inv
