"""Reading a molecule in, and saying out loud what was read.

Part 2 of the protocol is the step people rush. Nothing downstream can repair
a mistake made here, and most such mistakes are visible in ten seconds -- to
someone who looks. This module makes the tool look, and print what it saw, so
that the human judgement that follows is made against evidence rather than
against an assumption.

Supported inputs: PDB, mol2, SDF/MOL, and PDB chemical-component CIF, either
from a local file or fetched by three-character code. The CCD is the richest
source by a wide margin: besides coordinates it carries the authoritative
atom names, the formal charge, the parent residue of a modified amino acid,
and a per-atom flag marking exactly which atoms leave when the component is
polymerised. We use all four.
"""
from __future__ import annotations

import hashlib
import math
import os
import re
import urllib.request
from dataclasses import dataclass, field

CCD_URL = "https://files.rcsb.org/ligands/download/{code}.cif"

# Covalent radii (Angstrom), Cordero et al. 2008. Used only when the input
# format carries no connectivity (i.e. PDB).
COVALENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "I": 1.39, "Mo": 1.54, "Ru": 1.46, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
}

MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
          "Na": 22.990, "Mg": 24.305, "P": 30.974, "S": 32.06, "Cl": 35.45,
          "K": 39.098, "Ca": 40.078, "Fe": 55.845, "Zn": 65.38,
          "Se": 78.971, "Br": 79.904, "I": 126.904}

# Atoms that make an alpha-amino-acid fragment a free amino acid rather than
# a chain unit. Used only as a fallback: the CCD leaving-atom flag is better.
FREE_TERMINUS_ATOMS = ("H2", "HN2", "OXT", "HXT", "HO", "H3", "HN3")


@dataclass
class Atom:
    name: str
    element: str
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    formal_charge: int = 0
    resname: str = ""
    resseq: int = 1
    aromatic: bool = False
    leaving: bool = False        # CCD pdbx_leaving_atom_flag == Y
    stereo: str = ""             # CCD pdbx_stereo_config: R / S / N

    def copy(self):
        return Atom(**self.__dict__)


@dataclass
class Structure:
    atoms: list = field(default_factory=list)
    bonds: list = field(default_factory=list)   # (i, j, order) 0-based
    name: str = ""
    title: str = ""
    meta: dict = field(default_factory=dict)
    source_path: str = ""
    source_sha256: str = ""
    bonds_perceived: bool = False               # True == inferred by distance

    # ------------------------------------------------------------------
    def index(self, name):
        for i, a in enumerate(self.atoms):
            if a.name == name:
                return i
        return -1

    def get(self, name):
        i = self.index(name)
        return self.atoms[i] if i >= 0 else None

    def names(self):
        return [a.name for a in self.atoms]

    def neighbours(self, i):
        out = []
        for a, b, _ in self.bonds:
            if a == i:
                out.append(b)
            elif b == i:
                out.append(a)
        return sorted(out)

    def formula(self):
        counts = {}
        for a in self.atoms:
            counts[a.element] = counts.get(a.element, 0) + 1
        order = ["C", "H"] + sorted(e for e in counts if e not in ("C", "H"))
        return " ".join(f"{e}{counts[e]}" for e in order if e in counts)

    def mass(self):
        return sum(MASSES.get(a.element, 0.0) for a in self.atoms)

    def n_hydrogens(self):
        return sum(1 for a in self.atoms if a.element == "H")

    def formal_charge(self):
        return sum(a.formal_charge for a in self.atoms)

    def duplicate_names(self):
        seen, dup = set(), []
        for a in self.atoms:
            if a.name in seen:
                dup.append(a.name)
            seen.add(a.name)
        return sorted(set(dup))

    def close_contacts(self, cutoff=0.7):
        out = []
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                if _dist(self.atoms[i], self.atoms[j]) < cutoff:
                    out.append((self.atoms[i].name, self.atoms[j].name))
        return out

    def write_pdb(self, path, resname=None, resseq=1):
        resname = resname or self.name or "UNL"
        with open(path, "w") as fh:
            for n, a in enumerate(self.atoms, 1):
                nm = a.name
                # PDB atom-name column convention: element right-justified in
                # columns 13-14 unless the name is 4 characters.
                pname = (f"{nm:<4s}" if len(nm) >= 4
                         else f" {nm:<3s}")
                # cols 13-16 name, 17 altLoc, 18-20 resName, 22 chainID,
                # 23-26 resSeq, 31-38 x
                fh.write(
                    f"ATOM  {n:5d} {pname} {resname:>3s} A{resseq:4d}    "
                    f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00"
                    f"          {a.element:>2s}\n")
            # CONECT records preserve perceived or declared connectivity so
            # antechamber does not have to re-guess it from geometry.
            adj = {}
            for i, j, _ in self.bonds:
                adj.setdefault(i, []).append(j)
                adj.setdefault(j, []).append(i)
            for i in sorted(adj):
                partners = sorted(adj[i])
                for k in range(0, len(partners), 4):
                    chunk = partners[k:k + 4]
                    fh.write("CONECT" + f"{i + 1:5d}"
                             + "".join(f"{p + 1:5d}" for p in chunk) + "\n")
            fh.write("END\n")
        return path


def _dist(a, b):
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# minimal CIF reader (enough for a chemical-component file)
# ---------------------------------------------------------------------------
def _cif_tokens(line):
    out, i, n = [], 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c in "'\"":
            j = i + 1
            while j < n and line[j] != c:
                j += 1
            out.append(line[i + 1:j])
            i = j + 1
        else:
            j = i
            while j < n and line[j] not in " \t":
                j += 1
            out.append(line[i:j])
            i = j
    return out


def parse_cif(text):
    """Return (singles, loops) from a small CIF block."""
    singles, loops = {}, []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        raw = lines[i]
        s = raw.strip()
        if not s or s.startswith("#") or s.startswith("data_"):
            i += 1
            continue
        if s == "loop_":
            i += 1
            tags = []
            while i < len(lines) and lines[i].strip().startswith("_"):
                tags.append(lines[i].strip())
                i += 1
            rows = []
            while i < len(lines):
                ln = lines[i]
                st = ln.strip()
                if not st or st.startswith("#"):
                    i += 1
                    if st.startswith("#"):
                        break
                    continue
                if st.startswith("_") or st == "loop_" or st.startswith("data_"):
                    break
                if st.startswith(";"):
                    buf = [st[1:]]
                    i += 1
                    while i < len(lines) and not lines[i].startswith(";"):
                        buf.append(lines[i])
                        i += 1
                    i += 1
                    rows.append(["\n".join(buf)])
                    continue
                rows.append(_cif_tokens(ln))
                i += 1
            fixed = []
            for r in rows:
                if len(r) == len(tags):
                    fixed.append(dict(zip(tags, r)))
            loops.append((tags, fixed))
            continue
        if s.startswith("_"):
            toks = _cif_tokens(s)
            if len(toks) >= 2:
                singles[toks[0]] = " ".join(toks[1:])
            else:
                # value on the following line(s)
                i += 1
                if i < len(lines) and lines[i].strip().startswith(";"):
                    buf = [lines[i].strip()[1:]]
                    i += 1
                    while i < len(lines) and not lines[i].startswith(";"):
                        buf.append(lines[i])
                        i += 1
                    singles[toks[0]] = "\n".join(buf).strip()
                elif i < len(lines):
                    v = _cif_tokens(lines[i].strip())
                    singles[toks[0]] = " ".join(v)
            i += 1
            continue
        i += 1
    return singles, loops


def _clean(v):
    return "" if v in (None, "?", ".", "") else v


def read_ccd_cif(path_or_text, prefer="ideal"):
    """Parse a PDB chemical-component CIF."""
    if os.path.exists(str(path_or_text)):
        text = open(path_or_text, errors="replace").read()
        src = str(path_or_text)
    else:
        text, src = str(path_or_text), ""
    singles, loops = parse_cif(text)
    st = Structure(source_path=src)
    st.name = _clean(singles.get("_chem_comp.id", ""))
    st.title = _clean(singles.get("_chem_comp.name", ""))
    st.meta["ccd_type"] = _clean(singles.get("_chem_comp.type", ""))
    st.meta["ccd_formula"] = _clean(singles.get("_chem_comp.formula", ""))
    st.meta["ccd_formal_charge"] = _clean(
        singles.get("_chem_comp.pdbx_formal_charge",
                    singles.get("_chem_comp.formal_charge", "")))
    st.meta["parent_comp_id"] = _clean(
        singles.get("_chem_comp.mon_nstd_parent_comp_id", ""))
    st.meta["one_letter"] = _clean(
        singles.get("_chem_comp.one_letter_code", ""))

    atom_rows, bond_rows = [], []
    for tags, rows in loops:
        if any(t.startswith("_chem_comp_atom.") for t in tags):
            atom_rows = rows
        elif any(t.startswith("_chem_comp_bond.") for t in tags):
            bond_rows = rows

    xk = ("_chem_comp_atom.pdbx_model_Cartn_x_ideal",
          "_chem_comp_atom.model_Cartn_x")
    if prefer != "ideal":
        xk = xk[::-1]
    for r in atom_rows:
        key = next((k for k in xk if _clean(r.get(k))), None)
        pref = key.rsplit("_", 1)[0] if key else None
        if pref and key.endswith("_ideal"):
            kx, ky, kz = (key, key.replace("_x_", "_y_"),
                          key.replace("_x_", "_z_"))
        elif key:
            kx, ky, kz = key, key.replace("_x", "_y"), key.replace("_x", "_z")
        else:
            kx = ky = kz = None
        try:
            xyz = ((float(r[kx]), float(r[ky]), float(r[kz]))
                   if kx else (0.0, 0.0, 0.0))
        except (KeyError, ValueError):
            xyz = (0.0, 0.0, 0.0)
        chg = _clean(r.get("_chem_comp_atom.charge", "0")) or "0"
        st.atoms.append(Atom(
            name=_clean(r.get("_chem_comp_atom.atom_id", "")),
            element=_clean(r.get("_chem_comp_atom.type_symbol", "")).capitalize(),
            x=xyz[0], y=xyz[1], z=xyz[2],
            formal_charge=int(float(chg)),
            resname=st.name,
            aromatic=_clean(
                r.get("_chem_comp_atom.pdbx_aromatic_flag", "N")) == "Y",
            leaving=_clean(
                r.get("_chem_comp_atom.pdbx_leaving_atom_flag", "N")) == "Y",
            stereo=_clean(r.get("_chem_comp_atom.pdbx_stereo_config", "")),
        ))
    order_map = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 4, "QUAD": 4}
    for r in bond_rows:
        i = st.index(_clean(r.get("_chem_comp_bond.atom_id_1", "")))
        j = st.index(_clean(r.get("_chem_comp_bond.atom_id_2", "")))
        if i >= 0 and j >= 0:
            o = order_map.get(
                _clean(r.get("_chem_comp_bond.value_order", "SING")), 1)
            if _clean(r.get("_chem_comp_bond.pdbx_aromatic_flag", "N")) == "Y":
                o = 4
            st.bonds.append((i, j, o))
    return st


def fetch_ccd(code, cache_dir=".solerita-cache", prefer="ideal", timeout=60):
    """Download a chemical-component file, caching it on disk."""
    code = code.strip().upper()
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{code}.cif")
    if not os.path.exists(path):
        url = CCD_URL.format(code=code)
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = resp.read()
        if b"_chem_comp." not in data:
            raise ValueError(f"{url} did not return a chemical component; "
                             f"is {code!r} a valid CCD code?")
        with open(path, "wb") as fh:
            fh.write(data)
    st = read_ccd_cif(path, prefer=prefer)
    st.source_path = path
    st.source_sha256 = sha256_file(path)
    return st


# ---------------------------------------------------------------------------
# PDB / mol2 / SDF
# ---------------------------------------------------------------------------
def read_pdb(path, perceive=True):
    st = Structure(source_path=str(path), source_sha256=sha256_file(path))
    serial_to_idx = {}
    conect = []
    with open(path, errors="replace") as fh:
        for line in fh:
            rec = line[:6]
            if rec in ("ATOM  ", "HETATM"):
                try:
                    serial = int(line[6:11])
                except ValueError:
                    serial = len(st.atoms) + 1
                name = line[12:16].strip()
                altloc = line[16:17].strip()
                if altloc not in ("", "A", "1"):
                    continue
                resname = line[17:20].strip()
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    resseq = 1
                x, y, z = (float(line[30:38]), float(line[38:46]),
                           float(line[46:54]))
                el = line[76:78].strip().capitalize()
                if not el:
                    el = _guess_element(name)
                fc = 0
                raw = line[78:80].strip()
                if raw:
                    m = re.match(r"^(\d)([+-])$|^([+-])(\d)$", raw)
                    if m:
                        g = m.groups()
                        mag = int(g[0] or g[3])
                        sign = g[1] or g[2]
                        fc = mag if sign == "+" else -mag
                serial_to_idx[serial] = len(st.atoms)
                st.atoms.append(Atom(name=name, element=el, x=x, y=y, z=z,
                                     formal_charge=fc, resname=resname,
                                     resseq=resseq))
            elif rec == "CONECT":
                nums = [line[i:i + 5].strip() for i in range(6, len(line.rstrip()), 5)]
                nums = [int(n) for n in nums if n.isdigit()]
                if len(nums) >= 2:
                    conect.append(nums)
            elif rec.startswith("TITLE"):
                st.title = (st.title + " " + line[10:].strip()).strip()
    if st.atoms:
        st.name = st.atoms[0].resname
    seen = set()
    for nums in conect:
        a = serial_to_idx.get(nums[0])
        for s in nums[1:]:
            b = serial_to_idx.get(s)
            if a is not None and b is not None and a != b:
                key = (min(a, b), max(a, b))
                if key not in seen:
                    seen.add(key)
                    st.bonds.append((key[0], key[1], 1))
    if not st.bonds and perceive:
        perceive_bonds(st)
    return st


def _guess_element(name):
    s = name.strip()
    s = re.sub(r"^[0-9]+", "", s)
    if not s:
        return "X"
    two = s[:2].capitalize()
    if two in COVALENT_RADII and two not in ("Na", "Ca") or two in (
            "Cl", "Br", "Zn", "Fe", "Mg", "Mn", "Se", "Si"):
        return two
    return s[0].upper()


def perceive_bonds(st, tolerance=0.45):
    """Distance-based connectivity, for inputs that carry none."""
    st.bonds = []
    n = len(st.atoms)
    for i in range(n):
        ri = COVALENT_RADII.get(st.atoms[i].element, 0.77)
        for j in range(i + 1, n):
            rj = COVALENT_RADII.get(st.atoms[j].element, 0.77)
            d = _dist(st.atoms[i], st.atoms[j])
            if 0.4 < d < ri + rj + tolerance:
                if st.atoms[i].element == "H" and st.atoms[j].element == "H":
                    continue
                st.bonds.append((i, j, 1))
    st.bonds_perceived = True
    return st


def read_mol2(path):
    st = Structure(source_path=str(path), source_sha256=sha256_file(path))
    section = None
    with open(path, errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                section = s.split(">", 1)[1].strip().upper()
                continue
            if not s or s.startswith("#"):
                continue
            if section == "MOLECULE" and not st.name:
                st.name = s.split()[0]
            elif section == "ATOM":
                f = s.split()
                if len(f) < 6:
                    continue
                atype = f[5]
                el = atype.split(".")[0].capitalize()
                st.atoms.append(Atom(
                    name=f[1], element=el,
                    x=float(f[2]), y=float(f[3]), z=float(f[4]),
                    resname=(f[7] if len(f) > 7 else st.name),
                    resseq=int(f[6]) if len(f) > 6 and f[6].isdigit() else 1))
            elif section == "BOND":
                f = s.split()
                if len(f) >= 4:
                    o = {"1": 1, "2": 2, "3": 3, "ar": 4, "am": 1}.get(f[3], 1)
                    st.bonds.append((int(f[1]) - 1, int(f[2]) - 1, o))
    return st


def read_sdf(path):
    st = Structure(source_path=str(path), source_sha256=sha256_file(path))
    lines = open(path, errors="replace").read().splitlines()
    if len(lines) < 4:
        raise ValueError(f"{path}: not an SDF/MOL file")
    st.title = lines[0].strip()
    counts = lines[3]
    na, nb = int(counts[0:3]), int(counts[3:6])
    elcount = {}
    for k in range(4, 4 + na):
        ln = lines[k]
        x, y, z = float(ln[0:10]), float(ln[10:20]), float(ln[20:30])
        el = ln[31:34].strip().capitalize()
        elcount[el] = elcount.get(el, 0) + 1
        st.atoms.append(Atom(name=f"{el.upper()}{elcount[el]}", element=el,
                             x=x, y=y, z=z))
    for k in range(4 + na, 4 + na + nb):
        ln = lines[k]
        i, j, o = int(ln[0:3]) - 1, int(ln[3:6]) - 1, int(ln[6:9])
        st.bonds.append((i, j, 4 if o == 4 else o))
    for ln in lines[4 + na + nb:]:
        if ln.startswith("M  CHG"):
            f = ln.split()
            for p in range(4, len(f), 2):
                st.atoms[int(f[p - 1]) - 1].formal_charge = int(f[p])
    return st


def ccd_reference(st, cache_dir=None):
    """Best-effort: the CCD entry this structure appears to be derived from.

    The common workflow is to fetch a component, edit it (delete a proton,
    change a protonation state) and feed the edited file back. A PDB carries
    no formal charges, so on its own the edited file cannot say what species
    it is -- but diffed against the entry it came from, it says exactly:
    these atoms were removed, those were added. Returns None if the residue
    name is not a plausible component code or the fetch fails; this is an
    enrichment, never a requirement.
    """
    names = [a.resname for a in st.atoms if a.resname]
    if not names:
        return None
    code = max(set(names), key=names.count).strip().upper()
    if not code or not code.isalnum() or len(code) > 5:
        return None
    try:
        return fetch_ccd(code, cache_dir=cache_dir)
    except Exception:
        return None


def adopt_reference_flags(st, reference):
    """Copy per-atom annotations a hand-edited file cannot carry.

    A PDB has no leaving-atom flags and no stereo configuration. The CCD
    entry the file was edited from has both, and the atoms are matched by
    name, so they can be carried across for every atom that survived the
    edit. Returns (n_leaving, n_stereo) actually copied.
    """
    if reference is None:
        return (0, 0)
    ref = {a.name: a for a in reference.atoms}
    nl = ns = 0
    for a in st.atoms:
        r = ref.get(a.name)
        if r is None:
            continue
        if r.leaving and not a.leaving:
            a.leaving = True
            nl += 1
        if r.stereo and r.stereo != "N" and not a.stereo:
            a.stereo = r.stereo
            ns += 1
    return (nl, ns)


def n_rotatable(st):
    """Single bonds between heavy atoms that can actually turn.

    Both ends must carry another heavy atom (so terminal groups and lone
    hydroxyl hydrogens do not count), and the bond must not be in a ring.
    This is a flexibility indicator for the report, not a torsion list: it
    says how much the starting conformation is likely to matter.
    """
    heavy = {i for i, a in enumerate(st.atoms) if a.element != "H"}
    adj = {}
    for i, j, _o in st.bonds:
        if i in heavy and j in heavy:
            adj.setdefault(i, set()).add(j)
            adj.setdefault(j, set()).add(i)
    n = 0
    for i, j, order in st.bonds:
        if order != 1 or i not in heavy or j not in heavy:
            continue
        if len(adj.get(i, ())) < 2 or len(adj.get(j, ())) < 2:
            continue
        # In a ring iff j is still reachable from i without this bond.
        seen, stack = {i}, [i]
        reach = False
        while stack and not reach:
            cur = stack.pop()
            for nb in adj.get(cur, ()):
                if (cur, nb) in ((i, j), (j, i)):
                    continue
                if nb == j:
                    reach = True
                    break
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        if not reach:
            n += 1
    return n


READERS = {".pdb": read_pdb, ".ent": read_pdb, ".mol2": read_mol2,
           ".sdf": read_sdf, ".mol": read_sdf, ".cif": read_ccd_cif}


def read_structure(path):
    ext = os.path.splitext(str(path))[1].lower()
    if ext not in READERS:
        raise ValueError(f"{path}: unsupported format {ext!r}; "
                         f"use one of {sorted(READERS)}")
    st = READERS[ext](path)
    if not st.source_sha256:
        st.source_sha256 = sha256_file(path)
    st.source_path = str(path)
    return st
