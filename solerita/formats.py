"""Reading and editing the AmberTools file formats, without disturbing them.

Two jobs. Reading: pull out the few numbers each checkpoint needs -- charges,
atom types, topology letters, penalty scores -- so a check can be an assertion
instead of a human squinting at a log. Editing: change exactly one field of
one line and leave every other byte alone, because these formats are
column-sensitive and a reformatted file is a new source of error.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

_TRAILING = re.compile(r"^(.*?)(\s+)(\S+)\s*$")


def _replace_last_token(line, new):
    """Swap the final whitespace-delimited token, preserving line width."""
    m = _TRAILING.match(line.rstrip("\n"))
    if not m:
        return line
    prefix, sep, old = m.groups()
    width = len(sep) + len(old)
    return prefix + new.rjust(width) + "\n"


def _replace_token(line, index, new):
    """Swap the index-th whitespace-delimited token, preserving its width."""
    out, pos, n = [], 0, 0
    text = line.rstrip("\n")
    for m in re.finditer(r"\S+", text):
        if n == index:
            width = m.end() - m.start()
            return (text[:m.start()] + new.rjust(width)
                    + text[m.end():] + "\n")
        n += 1
    return line


# ---------------------------------------------------------------------------
# mol2
# ---------------------------------------------------------------------------
@dataclass
class Mol2Atom:
    serial: int
    name: str
    x: float
    y: float
    z: float
    atype: str
    resseq: int
    resname: str
    charge: float
    lineno: int = -1


def parse_mol2(path):
    atoms, section = [], None
    lines = open(path, errors="replace").read().splitlines(keepends=True)
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("@<TRIPOS>"):
            section = s.split(">", 1)[1].strip().upper()
            continue
        if section == "ATOM" and s:
            f = s.split()
            if len(f) >= 9:
                atoms.append(Mol2Atom(int(f[0]), f[1], float(f[2]),
                                      float(f[3]), float(f[4]), f[5],
                                      int(f[6]), f[7], float(f[8]), i))
    return atoms, lines


def rewrite_mol2_types(path, out_path, type_by_name):
    atoms, lines = parse_mol2(path)
    changed = []
    for a in atoms:
        new = type_by_name.get(a.name)
        if new and new != a.atype:
            lines[a.lineno] = _replace_token(lines[a.lineno], 5, new)
            changed.append((a.name, a.atype, new))
    with open(out_path, "w") as fh:
        fh.writelines(lines)
    return changed


# ---------------------------------------------------------------------------
# ac
# ---------------------------------------------------------------------------
@dataclass
class AcAtom:
    serial: int
    name: str
    resname: str
    resseq: int
    x: float
    y: float
    z: float
    charge: float
    atype: str
    lineno: int


# The .ac ATOM record is column-oriented: the atom name occupies columns
# 14-17 (0-based 13:17) and the residue name 18-20. Splitting on whitespace
# works only while every atom name is three characters or fewer; a four-
# character name (HOP2, HOP3, and every other phosphate or long-branch PTM)
# butts straight against the residue name, the token count drops by one, and
# every field after it silently shifts. Parse by column.
_AC_NAME = slice(13, 17)
_AC_RESNAME = slice(17, 20)
_AC_REST = 20


def parse_ac(path):
    atoms, bonds = [], []
    lines = open(path, errors="replace").read().splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith("ATOM"):
            name = line[_AC_NAME].strip()
            resname = line[_AC_RESNAME].strip()
            rest = line[_AC_REST:].split()
            if not name or len(rest) < 6:
                f = line.split()           # tolerate a non-standard writer
                name, resname, rest = f[2], f[3], f[4:]
            atoms.append(AcAtom(int(line[4:11]), name, resname,
                                int(rest[0]), float(rest[1]),
                                float(rest[2]), float(rest[3]),
                                float(rest[4]), rest[5], i))
        elif line.startswith("BOND"):
            f = line.split()
            bonds.append((int(f[2]), int(f[3]), int(f[4])))
    return atoms, bonds, lines


def rewrite_ac_types(path, out_path, type_by_name):
    atoms, _, lines = parse_ac(path)
    changed, missing = [], set(type_by_name)
    for a in atoms:
        new = type_by_name.get(a.name)
        if new is None:
            continue
        missing.discard(a.name)
        if new != a.atype:
            lines[a.lineno] = _replace_last_token(lines[a.lineno], new)
            changed.append((a.name, a.atype, new))
    with open(out_path, "w") as fh:
        fh.writelines(lines)
    return changed, sorted(missing)


def rename_ac_atoms(path, out_path, rename):
    atoms, _, lines = parse_ac(path)
    done = []
    for a in atoms:
        new = rename.get(a.name)
        if new and new != a.name:
            if len(new) > 4:
                raise ValueError(f"atom name {new!r} does not fit the "
                                 f"four-column name field")
            ln = lines[a.lineno]
            lines[a.lineno] = ln[:13] + new.ljust(4) + ln[17:]
            done.append((a.name, new))
    with open(out_path, "w") as fh:
        fh.writelines(lines)
    return done


# ---------------------------------------------------------------------------
# prepi
# ---------------------------------------------------------------------------
@dataclass
class PrepAtom:
    index: int
    name: str
    atype: str
    topo: str
    na: int
    nb: int
    nc: int
    r: float
    theta: float
    phi: float
    charge: float
    lineno: int


@dataclass
class Prepi:
    atoms: list = field(default_factory=list)
    impropers: list = field(default_factory=list)
    loops: list = field(default_factory=list)
    resname: str = ""
    lines: list = field(default_factory=list)

    @property
    def real(self):
        return [a for a in self.atoms if a.atype != "DU"]

    def charge_sum(self):
        return sum(a.charge for a in self.real)


def parse_prepi(path):
    p = Prepi()
    p.lines = open(path, errors="replace").read().splitlines(keepends=True)
    section = "head"
    for i, line in enumerate(p.lines):
        s = line.strip()
        if s in ("LOOP",):
            section = "loop"
            continue
        if s in ("IMPROPER",):
            section = "improper"
            continue
        if s in ("DONE", "STOP"):
            section = "head"
            continue
        f = s.split()
        if section == "head" and len(f) >= 11 and f[0].isdigit():
            try:
                p.atoms.append(PrepAtom(
                    int(f[0]), f[1], f[2], f[3], int(f[4]), int(f[5]),
                    int(f[6]), float(f[7]), float(f[8]), float(f[9]),
                    float(f[10]), i))
            except ValueError:
                pass
        elif section == "improper" and len(f) == 4:
            p.impropers.append(tuple(f))
        elif section == "loop" and len(f) == 2:
            p.loops.append(tuple(f))
        elif (section == "head" and len(f) >= 3
              and f[1].upper() in ("INT", "XYZ") and not p.resname):
            p.resname = f[0]
    return p


def set_prepi_charges(p, charge_by_name, out_path):
    """Write the prepi back with new charges, touching nothing else."""
    lines = list(p.lines)
    changed = []
    for a in p.atoms:
        q = charge_by_name.get(a.name)
        if q is None or abs(q - a.charge) < 1e-9:
            continue
        lines[a.lineno] = _replace_last_token(lines[a.lineno], f"{q:.6f}")
        changed.append((a.name, a.charge, q))
    with open(out_path, "w") as fh:
        fh.writelines(lines)
    return changed


def check_prepi_tree(p):
    """Verify the topology letter of every atom against the build tree.

    The letter counts how many children an atom has in the tree that the
    internal coordinates define -- hydrogens included -- and adding an atom
    changes its parent's letter. This is the error the protocol warns about
    at the sulfoxide sulfur (S gains OE and must go from S to B), and it is
    exactly the kind of thing a human misses and nothing else checks.

    Main-chain atoms are marked M regardless of their child count, so their
    letter is not derivable and is reported as unchecked.
    """
    children = {a.index: 0 for a in p.atoms}
    for a in p.atoms:
        if a.na in children and a.index > 3:
            children[a.na] += 1
    expect = {0: "E", 1: "S", 2: "B", 3: "3", 4: "4"}
    issues, checked = [], 0
    for a in p.atoms:
        if a.atype == "DU":
            continue
        n = children.get(a.index, 0)
        if a.topo == "M":
            continue
        checked += 1
        want = expect.get(n)
        if want is None:
            issues.append((a.name, a.topo, f"{n} children (unexpected)"))
        elif a.topo != want:
            issues.append((a.name, a.topo,
                           f"has {n} children, so the letter should be "
                           f"{want!r}, not {a.topo!r}"))
    return issues, checked


# ---------------------------------------------------------------------------
# frcmod
# ---------------------------------------------------------------------------
_PENALTY = re.compile(r"penalty score\s*=\s*([-\d.]+)", re.I)
SECTIONS = ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON", "LJEDIT",
            "CMAP")


@dataclass
class FrcTerm:
    section: str
    key: str
    rest: str
    comment: str
    penalty: float = 0.0
    attn: bool = False
    line: str = ""

    @property
    def band(self):
        if self.attn:
            return "ATTN"
        if self.penalty == 0:
            return "exact"
        if self.penalty <= 10:
            return "acceptable"
        if self.penalty <= 50:
            return "reach"
        return "poor"


def parse_frcmod(path, acceptable=10.0, reach=50.0):
    terms, section = [], None
    for line in open(path, errors="replace"):
        s = line.rstrip("\n")
        head = s.strip().split()[0].upper() if s.strip() else ""
        if head in SECTIONS and len(s.strip().split()) <= 2:
            section = head
            continue
        if not s.strip() or section is None:
            continue
        m = _PENALTY.search(s)
        pen = float(m.group(1)) if m else 0.0
        attn = "ATTN" in s.upper()
        # The key is the leading type field: for BOND/ANGLE/DIHE it is the
        # dash-joined atom types, which may contain spaces around the dashes.
        key = s[:11].strip() if section != "MASS" else s.split()[0]
        comment = ""
        for marker in ("same as", "Using", "ATTN"):
            idx = s.find(marker)
            if idx >= 0:
                comment = s[idx:].strip()
                break
        terms.append(FrcTerm(section, key, s[11:].strip(), comment, pen,
                             attn, s))
    bands = {}
    for t in terms:
        if t.penalty == 0 and not t.attn:
            b = "exact"
        elif t.attn:
            b = "ATTN"
        elif t.penalty <= acceptable:
            b = "acceptable"
        elif t.penalty <= reach:
            b = "reach"
        else:
            b = "poor"
        t.penalty_band = b
        bands.setdefault(b, []).append(t)
    return terms, bands


def strip_attn(in_path, out_path):
    """The protocol's `grep -v ATTN`, plus the orphans grep leaves behind.

    `grep -v ATTN` is not sufficient, and the residue it leaves is not
    harmless. For a torsion it cannot parametrize, parmchk2 writes TWO
    lines: the flagged placeholder, and an unflagged companion with idivf 0
    and periodicity 0.000 --

        h1-c3-ss-sh   1    0.000   0.000   2.000   ATTN, need revision
        h1-c3-ss-sh   0    0.000   0.000   0.000

    Dropping only the first leaves the second, which LEaP reads as a real
    and *specific* torsion type. A specific type beats the generic
    X -b-c-X  that GAFF2 would otherwise have supplied, so the barrier is
    silently replaced by zero rather than falling back to the force field.
    Measured on S-mercaptocysteine, that zeroed 8 of 97 torsions in the
    built topology -- every rotation barrier involving a CB hydrogen.

    A torsion of periodicity zero has no meaning in AMBER's functional
    form, so removing these is safe regardless of how they arose, and lets
    the generic parameter apply as intended.

    Returns (n_attn, n_null).
    """
    kept, n_attn, n_null = [], 0, 0
    section = None
    for line in open(in_path, errors="replace"):
        head = line.strip().split()[0].upper() if line.strip() else ""
        if head in SECTIONS and len(line.strip().split()) <= 2:
            section = head
            kept.append(line)
            continue
        if "ATTN" in line.upper():
            n_attn += 1
            continue
        if section == "DIHE" and len(line) > 11:
            f = line[11:].split()
            try:
                if float(f[3]) == 0.0:
                    n_null += 1
                    continue
            except (IndexError, ValueError):
                pass
        kept.append(line)
    with open(out_path, "w") as fh:
        fh.writelines(kept)
    return n_attn, n_null
