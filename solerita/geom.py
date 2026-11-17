"""Placing an atom from internal coordinates.

Used by the capping step: given three placed atoms and a (bond, angle,
dihedral) triple, return the Cartesian position of the fourth. This is the
same construction the prepi format itself uses, which is convenient -- the
numbers we put in are directly comparable with the numbers prepgen writes out.
"""
from __future__ import annotations

import math


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a):
    n = math.sqrt(_dot(a, a))
    if n < 1e-12:
        raise ValueError("cannot normalise a zero-length vector; three "
                         "reference atoms are collinear or coincident")
    return _scale(a, 1.0 / n)


def place(c, b, a, bond, angle_deg, dihedral_deg):
    """Position of atom D given C, B, A and the internals D-C, D-C-B, D-C-B-A.

    c, b, a are (x, y, z) tuples; the new atom is bonded to c.
    """
    # place() already builds in the IUPAC sense; dihedral() below carries the
    # compensating sign so that measuring a built atom returns what was asked.
    th, ph = math.radians(angle_deg), math.radians(dihedral_deg)
    bc = _norm(_sub(c, b))
    ab = _sub(b, a)
    n = _cross(ab, bc)
    if math.sqrt(_dot(n, n)) < 1e-8:
        # A, B, C collinear: pick any perpendicular so the result is defined
        tmp = (1.0, 0.0, 0.0) if abs(bc[0]) < 0.9 else (0.0, 1.0, 0.0)
        n = _cross(tmp, bc)
    n = _norm(n)
    m = _cross(n, bc)
    d2 = (-bond * math.cos(th),
          bond * math.cos(ph) * math.sin(th),
          bond * math.sin(ph) * math.sin(th))
    return _add(c, (bc[0] * d2[0] + m[0] * d2[1] + n[0] * d2[2],
                    bc[1] * d2[0] + m[1] * d2[1] + n[1] * d2[2],
                    bc[2] * d2[0] + m[2] * d2[1] + n[2] * d2[2]))


def distance(p, q):
    return math.sqrt(sum((p[i] - q[i]) ** 2 for i in range(3)))


def angle(p, q, r):
    v1, v2 = _sub(p, q), _sub(r, q)
    c = _dot(v1, v2) / (math.sqrt(_dot(v1, v1)) * math.sqrt(_dot(v2, v2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def dihedral(p, q, r, s):
    b1, b2, b3 = _sub(q, p), _sub(r, q), _sub(s, r)
    n1, n2 = _cross(b1, b2), _cross(b2, b3)
    m = _cross(n1, _norm(b2))
    x, y = _dot(n1, n2), _dot(m, n2)
    return -math.degrees(math.atan2(y, x))
