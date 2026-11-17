"""The protocol, as plain settings.

Every choice the tutorial makes is written down here, once, with the reason
for it. Nothing else in solerita decides anything scientific: the other files
only carry these choices out and check the results.

If you want to change the protocol, this is the only file to edit. If you do,
say so when you report results -- parameters made with different settings are
not directly comparable.
"""

# --- Part 3: charges and atom types ------------------------------------------

# AM1-BCC: a fast semi-empirical calculation (AM1) followed by bond-charge
# corrections fitted to reproduce higher-level (HF/6-31G*) electrostatics.
# The usual alternative, RESP, needs a full quantum-chemistry program; AM1-BCC
# comes with AmberTools. It is the standard choice for GAFF/GAFF2 ligands.
CHARGE_METHOD = "bcc"

# GAFF2: the General AMBER Force Field, version 2. It has atom types for
# almost any organic molecule (lower-case names: c3, ss, oh, ...), which is why
# we use it for the parts the protein force field does not cover.
ATOM_TYPES = "gaff2"

# --- Part 4: modified residues keep protein types on the backbone -------------

# GAFF2 types the whole molecule, backbone included. But a residue must join a
# protein built with ff14SB, so its backbone has to look exactly like every
# other residue's backbone. These are the ff14SB types for an amino-acid
# backbone. The alpha hydrogen is H1 (an aliphatic hydrogen next to one
# electronegative atom), as in every standard residue in amino12.lib.
BACKBONE_TYPES = {"N": "N", "H": "H", "CA": "CX", "HA": "H1",
                  "C": "C", "O": "O"}

# --- Part 2.3 and Part 6: capping and joining a chain -------------------------

# A residue is parametrized as ACE-X-NME: capped on both sides so its ends
# look like they are inside a protein. The cap atoms are given these names so
# that prepgen can be told exactly which atoms to throw away afterwards.
CAP_NAMES = {
    "ACE": {"C": "CY", "O": "OY", "CH3": "CAY",
            "H1": "HY1", "H2": "HY2", "H3": "HY3"},
    "NME": {"N": "NZ", "H": "HZ", "CH3": "CAT",
            "H1": "HT1", "H2": "HT2", "H3": "HT3"},
}
HEAD_ATOM = "N"            # bonds to the previous residue
TAIL_ATOM = "C"            # bonds to the next residue
MAIN_CHAIN = ["CA"]        # atoms on the path from HEAD to TAIL
PRE_HEAD_TYPE = "C"        # atom type of the previous residue's C
POST_TAIL_TYPE = "N"       # atom type of the next residue's N

# --- Part 7: missing parameters ------------------------------------------------

PARMCHK2_GAFF = 2          # parmchk2 -s 2 means "use GAFF2", not GAFF

# parmchk2 marks every analogy with a penalty score. Up to 10 is a close
# analogue; up to 50 is a reach worth looking at; above that, check it.
PENALTY_OK = 10.0
PENALTY_REACH = 50.0

# --- Part 8: validation -------------------------------------------------------

PROTEIN_FF = "leaprc.protein.ff14SB"
GAFF_LEAPRC = "leaprc.gaff2"
CHARGE_TOLERANCE = 1.0e-4      # |sum of charges - net charge| must be below
# The .prepi stores charges to six decimals and prepgen spreads the removed
# caps' charge over the residue, so its sum carries a rounding residual of a
# few 1e-4 e. The residue library check therefore allows 1e-3.
PREPI_CHARGE_TOLERANCE = 1.0e-3
MIN_INPUT = ("Minimisation\n &cntrl\n"
             "  imin=1, maxcyc=500, ncyc=200, ntb=0, igb=0, cut=999.0,\n /\n")
