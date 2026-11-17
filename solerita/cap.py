"""Capping a residue fragment as ACE-X-NME, without a GUI.

Part 2.3 of the protocol asks you to open PyMOL, delete three atoms and build
two caps by hand. That is the one manual step in the residue branch that is
pure geometry: the atoms to delete are declared by the chemical component
dictionary, and the cap atoms are placed from standard peptide internal
coordinates. Doing it in code makes it reproducible and removes the commonest
source of silent error -- not the geometry but the atom naming that prepgen
later has to disambiguate.

One subtlety that makes this more than bolt-on-twelve-atoms. A free amino
acid has a pyramidal -NH2 and a -COOH; a chain unit has two planar amides.
So the backbone H on N and the carbonyl O on C are *not* in the right place
once the caps are attached, and leaving them where the input had them puts
the cap carbonyl on top of the amine hydrogen. Both are therefore rebuilt in
the amide plane together with the caps.

Caps are built directly with the names the profile requires (CY/OY/CAY/HY1-3
and NZ/HZ/CAT/HT1-3), so no renaming step follows. A structure the user
capped themselves is brought to the same naming by `rename_existing_caps`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .geom import place, distance, angle as g_angle, dihedral as g_dihedral
from .structure import Atom, Structure, FREE_TERMINUS_ATOMS

# Standard ff14SB peptide geometry (the numbers that appear in amino12.in).
R_N_C, A_C_N_CA = 1.335, 121.9        # N-CY bond, CY-N-CA angle
R_N_H, A_H_N_CA = 1.010, 119.8        # the amide H we rebuild on N
R_C_N, A_N_C_CA = 1.335, 116.6        # C-NZ bond, NZ-C-CA angle
R_C_O, A_O_C_CA = 1.229, 120.5        # the carbonyl O we rebuild on C
R_CO, A_O_C_N = 1.229, 122.9          # cap carbonyl
R_CC, A_CH3_C_N = 1.522, 116.6        # ACE methyl
R_NC, A_CH3_N_C = 1.449, 121.9        # NME methyl
R_NH, A_H_N_C = 1.010, 119.8          # NME amide H
R_CH, A_H_C_X = 1.090, 109.5          # methyl hydrogens

# An extended backbone: AM1-BCC charges are computed at the AM1-optimised
# geometry reached from the input, and a curled-up start stays curled up.
DEFAULT_PHI = -135.0
DEFAULT_PSI = 135.0
# Fallbacks tried, in order, if the default produces a steric clash with the
# sidechain conformation the input happens to carry.
PHI_PSI_FALLBACKS = [(-135.0, 135.0), (-150.0, 150.0), (-120.0, 120.0),
                     (-60.0, 135.0), (-135.0, 60.0), (-60.0, -45.0),
                     (180.0, 180.0), (-90.0, 90.0)]


class CappingError(ValueError):
    pass


def resolve_cap_scheme(existing_names, scheme):
    """A cap naming that collides with nothing in the residue.

    The protocol's default names are CY/OY/CAY/HY1-3 and NZ/HZ/CAT/HT1-3.
    Those are fine for most residues and wrong for a whole family of common
    ones: every lysine-derived modification (MLY, MLZ, M3L, ALY, ...) already
    owns NZ, and several own HZ. prepgen matches atoms by name and cannot
    disambiguate, so the cap must move, not the residue -- the residue's
    names have to match the PDB files you will later load.

    Returns (scheme, substitutions) where substitutions is a list of
    (role, default, used) for the report.
    """
    taken = {n.upper() for n in existing_names}
    out, subs = {}, []
    for resname in ("ACE", "NME"):
        out[resname] = {}
        for role, default in scheme[resname].items():
            name = default
            if name.upper() in taken:
                stem = "".join(ch for ch in default if not ch.isdigit())
                for suffix in list("123456789") + ["A", "B", "C"]:
                    cand = f"{stem[:3]}{suffix}"
                    if cand.upper() not in taken and len(cand) <= 4:
                        name = cand
                        break
                else:
                    raise CappingError(
                        f"cannot find a free cap name for {resname} {role!r}; "
                        f"the residue already uses {sorted(taken)}")
                subs.append((f"{resname}.{role}", default, name))
            taken.add(name.upper())
            out[resname][role] = name
    return out, subs


@dataclass
class CapReport:
    deleted: list = field(default_factory=list)
    added: list = field(default_factory=list)
    rebuilt: list = field(default_factory=list)
    scheme: dict = field(default_factory=dict)
    substitutions: list = field(default_factory=list)
    phi: float = DEFAULT_PHI
    psi: float = DEFAULT_PSI
    tried: int = 1
    measured: list = field(default_factory=list)
    clashes: list = field(default_factory=list)
    messages: list = field(default_factory=list)

    @staticmethod
    def _dev(label, want, got):
        """Deviation, comparing torsions circularly so -180 == +180."""
        if label.endswith(" dih"):
            return abs(((want - got + 180.0) % 360.0) - 180.0)
        return abs(want - got)

    @property
    def deviations(self):
        return [(l, w, g, self._dev(l, w, g)) for l, w, g in self.measured]

    @property
    def worst_internal(self):
        d = self.deviations
        return max((x[3] for x in d), default=0.0)

    def describe(self):
        L = [f"capping: ACE-X-NME, phi {self.phi:+.1f} psi {self.psi:+.1f}"
             + (f"  (default clashed; {self.tried} settings tried)"
                if self.tried > 1 else ""),
             f"  deleted   {' '.join(self.deleted) or '(nothing)'}",
             f"  added     {' '.join(self.added)}",
             f"  rebuilt   {' '.join(self.rebuilt) or '(nothing)'}"
             + "   (free amino acid -> amide geometry)",
             f"  geometry  {len(self.measured)} internals rebuilt, worst "
             f"deviation from the requested value {self.worst_internal:.3g}"]
        if self.substitutions:
            L.append("  RENAMED   default cap names collide with the "
                     "residue's own atoms: "
                     + ", ".join(f"{d}->{u}" for _, d, u
                                 in self.substitutions))
        if self.clashes:
            L.append(f"  CLASHES   {', '.join(self.clashes)}")
        L += [f"  {m}" for m in self.messages]
        return "\n".join(L)


def find_leaving_atoms(st, fallback=FREE_TERMINUS_ATOMS):
    flagged = [a.name for a in st.atoms if a.leaving]
    if flagged:
        return flagged, "CCD pdbx_leaving_atom_flag"
    present = [n for n in fallback if st.index(n) >= 0]
    return present, "name-based fallback (no CCD leaving flags available)"


def looks_precapped(st):
    resnames = {a.resname.upper() for a in st.atoms}
    if {"ACE", "NME"} & resnames:
        return True
    return st.index("CY") >= 0 and st.index("NZ") >= 0


def rename_existing_caps(st, scheme, report=None):
    """Give hand-built ACE/NME caps the unambiguous names prepgen needs."""
    renamed = []
    for resname in ("ACE", "NME"):
        mapping = scheme[resname]
        members = [a for a in st.atoms if a.resname.upper() == resname]
        h_seen = 0
        for a in members:
            old, key = a.name, a.name
            up = old.upper()
            if up.startswith("HH3") or (up.startswith("H") and len(up) > 1
                                        and resname == "ACE"):
                h_seen += 1
                key = f"H{h_seen}"
            elif up.startswith("HH3") or (up.startswith("H") and len(up) > 1):
                h_seen += 1
                key = f"H{h_seen}"
            new = mapping.get(key, mapping.get(old))
            if new and new != old:
                a.name = new
                renamed.append(f"{resname}:{old}->{new}")
    if report is not None:
        report.messages.append(
            f"renamed {len(renamed)} cap atoms: {', '.join(renamed)}"
            if renamed else "cap atoms already carry the protocol names")
    return renamed


def _strip(st, delete):
    keep_names = [a.name for a in st.atoms if a.name not in delete]
    new = Structure(name=st.name, title=st.title, meta=dict(st.meta),
                    source_path=st.source_path,
                    source_sha256=st.source_sha256,
                    bonds_perceived=st.bonds_perceived)
    new.atoms = [a.copy() for a in st.atoms if a.name in set(keep_names)]
    kept = set(keep_names)
    for i, j, o in st.bonds:
        ni, nj = st.atoms[i].name, st.atoms[j].name
        if ni in kept and nj in kept:
            new.bonds.append((new.index(ni), new.index(nj), o))
    return new


def _build_once(base, scheme, phi, psi, head, alpha, tail, carbonyl_o,
                amide_h):
    """Place caps and rebuild the junction H/O for one (phi, psi)."""
    P = {a.name: (a.x, a.y, a.z) for a in base.atoms}
    N, CA, C = P[head], P[alpha], P[tail]
    s = scheme
    new_pos, rebuilt = {}, {}

    # --- N-terminal side: ACE, and the amide H that replaces the amine ---
    CY = place(N, CA, C, R_N_C, A_C_N_CA, phi)
    new_pos[s["ACE"]["C"]] = CY
    if amide_h:
        # Third substituent on a now-planar N: 180 degrees from CY about CA-N.
        rebuilt[amide_h] = place(N, CA, C, R_N_H, A_H_N_CA, phi + 180.0)
    new_pos[s["ACE"]["O"]] = place(CY, N, CA, R_CO, A_O_C_N, 0.0)
    CAY = place(CY, N, CA, R_CC, A_CH3_C_N, 180.0)
    new_pos[s["ACE"]["CH3"]] = CAY
    for k, d in (("H1", 60.0), ("H2", 180.0), ("H3", 300.0)):
        new_pos[s["ACE"][k]] = place(CAY, CY, N, R_CH, A_H_C_X, d)

    # --- C-terminal side: NME, and the carbonyl O that replaces -COOH ----
    NZ = place(C, CA, N, R_C_N, A_N_C_CA, psi)
    new_pos[s["NME"]["N"]] = NZ
    if carbonyl_o:
        rebuilt[carbonyl_o] = place(C, CA, N, R_C_O, A_O_C_CA, psi + 180.0)
    new_pos[s["NME"]["H"]] = place(NZ, C, CA, R_NH, A_H_N_C, 0.0)
    CAT = place(NZ, C, CA, R_NC, A_CH3_N_C, 180.0)
    new_pos[s["NME"]["CH3"]] = CAT
    for k, d in (("H1", 60.0), ("H2", 180.0), ("H3", 300.0)):
        new_pos[s["NME"][k]] = place(CAT, NZ, C, R_CH, A_H_C_X, d)
    return new_pos, rebuilt


def _assemble(base, scheme, new_pos, rebuilt, head, tail):
    st = Structure(name=base.name, title=base.title, meta=dict(base.meta),
                   source_path=base.source_path,
                   source_sha256=base.source_sha256,
                   bonds_perceived=base.bonds_perceived)
    s = scheme
    ace_order = [s["ACE"][k] for k in ("C", "O", "CH3", "H1", "H2", "H3")]
    nme_order = [s["NME"][k] for k in ("N", "H", "CH3", "H1", "H2", "H3")]

    def mk(name, resname, resseq):
        x, y, z = new_pos[name]
        el = "H" if name.startswith("H") else name[0]
        return Atom(name=name, element=el, x=x, y=y, z=z,
                    resname=resname, resseq=resseq)

    mid = []
    for a in base.atoms:
        b = a.copy()
        b.resseq = 1
        b.leaving = False
        if b.name in rebuilt:
            b.x, b.y, b.z = rebuilt[b.name]
        mid.append(b)
    st.atoms = ([mk(n, "ACE", 0) for n in ace_order] + mid
                + [mk(n, "NME", 2) for n in nme_order])
    for i, j, o in base.bonds:
        st.bonds.append((st.index(base.atoms[i].name),
                         st.index(base.atoms[j].name), o))
    for a, b, o in (
            (head, s["ACE"]["C"], 1),
            (s["ACE"]["C"], s["ACE"]["O"], 2),
            (s["ACE"]["C"], s["ACE"]["CH3"], 1),
            (s["ACE"]["CH3"], s["ACE"]["H1"], 1),
            (s["ACE"]["CH3"], s["ACE"]["H2"], 1),
            (s["ACE"]["CH3"], s["ACE"]["H3"], 1),
            (tail, s["NME"]["N"], 1),
            (s["NME"]["N"], s["NME"]["H"], 1),
            (s["NME"]["N"], s["NME"]["CH3"], 1),
            (s["NME"]["CH3"], s["NME"]["H1"], 1),
            (s["NME"]["CH3"], s["NME"]["H2"], 1),
            (s["NME"]["CH3"], s["NME"]["H3"], 1)):
        st.bonds.append((st.index(a), st.index(b), o))
    return st


def _clashes(st, moved_names, cutoff=1.35):
    bonded = {(min(i, j), max(i, j)) for i, j, _ in st.bonds}
    one_three = set()
    for i in range(len(st.atoms)):
        nb = st.neighbours(i)
        for a in nb:
            for b in nb:
                if a < b:
                    one_three.add((a, b))
    out = []
    for i, a in enumerate(st.atoms):
        if a.name not in moved_names:
            continue
        for j, b in enumerate(st.atoms):
            if j == i:
                continue
            key = (min(i, j), max(i, j))
            if key in bonded or key in one_three:
                continue
            if distance((a.x, a.y, a.z), (b.x, b.y, b.z)) < cutoff:
                out.append(f"{a.name}-{b.name}")
    return sorted(set(out))


def cap_residue(st, scheme, phi=None, psi=None, head="N", tail="C",
                alpha="CA", delete=None, auto_relieve=True):
    """Return (capped Structure, CapReport). The input is not modified."""
    rep = CapReport()
    for nm in (head, alpha, tail):
        if st.index(nm) < 0:
            raise CappingError(
                f"atom {nm!r} not found, so this is not an alpha-amino-acid "
                f"fragment the tool knows how to cap. Cap it yourself and set "
                f"capping='precapped'. Atoms present: "
                f"{' '.join(st.names()[:20])}")

    if delete is None:
        delete, how = find_leaving_atoms(st)
        rep.messages.append(f"leaving atoms identified by {how}")
    delete = [d for d in delete if st.index(d) >= 0]
    if set(delete) & {head, tail, alpha}:
        raise CappingError(f"refusing to delete a backbone atom: {delete}")
    rep.deleted = list(delete)
    base = _strip(st, set(delete))

    scheme, subs = resolve_cap_scheme(base.names(), scheme)
    rep.scheme, rep.substitutions = scheme, subs

    hi = base.index(head)
    h_on_n = [base.atoms[k].name for k in base.neighbours(hi)
              if base.atoms[k].element == "H"]
    if len(h_on_n) > 1:
        raise CappingError(
            f"{head} still carries {len(h_on_n)} hydrogens ({h_on_n}) after "
            f"removing {delete}. A chain nitrogen carries one. If this is a "
            f"protonated N-terminus, neutralise it first.")
    amide_h = h_on_n[0] if h_on_n else None
    if amide_h is None:
        rep.messages.append(f"{head} carries no hydrogen (proline-like); "
                            f"none rebuilt")

    ti = base.index(tail)
    o_on_c = [base.atoms[k].name for k in base.neighbours(ti)
              if base.atoms[k].element == "O"]
    if len(o_on_c) != 1:
        raise CappingError(
            f"{tail} carries {len(o_on_c)} oxygens ({o_on_c}) after removing "
            f"{delete}; expected exactly the carbonyl O.")
    carbonyl_o = o_on_c[0]
    rep.rebuilt = [x for x in (amide_h, carbonyl_o) if x]

    trials = ([(phi if phi is not None else DEFAULT_PHI,
                psi if psi is not None else DEFAULT_PSI)]
              + (PHI_PSI_FALLBACKS if auto_relieve and phi is None
                 and psi is None else []))
    s = scheme
    cap_names = set(s["ACE"].values()) | set(s["NME"].values())
    moved = cap_names | set(rep.rebuilt)
    best = None
    for n, (ph, ps) in enumerate(trials, 1):
        new_pos, rebuilt = _build_once(base, scheme, ph, ps, head, alpha,
                                       tail, carbonyl_o, amide_h)
        cand = _assemble(base, scheme, new_pos, rebuilt, head, tail)
        cl = _clashes(cand, moved)
        if best is None or len(cl) < len(best[2]):
            best = (cand, (ph, ps), cl, n)
        if not cl:
            break
    cand, (ph, ps), cl, n = best
    rep.phi, rep.psi, rep.clashes, rep.tried = ph, ps, cl, n
    dup = cand.duplicate_names()
    if dup:
        raise CappingError(f"capping produced duplicate atom names {dup}; "
                           f"this is a bug in the cap naming, not in your "
                           f"input")
    rep.added = [s["ACE"][k] for k in ("C", "O", "CH3", "H1", "H2", "H3")] + \
                [s["NME"][k] for k in ("N", "H", "CH3", "H1", "H2", "H3")]

    P = {a.name: (a.x, a.y, a.z) for a in cand.atoms}
    N, CA, C = P[head], P[alpha], P[tail]
    checks = [(s["ACE"]["C"], (head, alpha, tail), R_N_C, A_C_N_CA, ph),
              (s["NME"]["N"], (tail, alpha, head), R_C_N, A_N_C_CA, ps),
              (s["ACE"]["CH3"], (s["ACE"]["C"], head, alpha), R_CC,
               A_CH3_C_N, 180.0),
              (s["ACE"]["O"], (s["ACE"]["C"], head, alpha), R_CO,
               A_O_C_N, 0.0),
              (s["NME"]["CH3"], (s["NME"]["N"], tail, alpha), R_NC,
               A_CH3_N_C, 180.0),
              (s["NME"]["H"], (s["NME"]["N"], tail, alpha), R_NH,
               A_H_N_C, 0.0)]
    if amide_h:
        checks.append((amide_h, (head, alpha, tail), R_N_H, A_H_N_CA,
                       ph + 180.0))
    checks.append((carbonyl_o, (tail, alpha, head), R_C_O, A_O_C_CA,
                   ps + 180.0))
    for nm, (rb, rc, rd), r, th, dh in checks:
        want = ((dh + 180.0) % 360.0) - 180.0
        rep.measured += [
            (f"{nm} r", r, distance(P[nm], P[rb])),
            (f"{nm} angle", th, g_angle(P[nm], P[rb], P[rc])),
            (f"{nm} dih", want, g_dihedral(P[nm], P[rb], P[rc], P[rd]))]
    return cand, rep
