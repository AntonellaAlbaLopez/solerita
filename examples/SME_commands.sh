#!/usr/bin/env bash
# Every command solerita ran to parametrize SME (residue, net charge +0).
# Read it next to the tutorial to see what each step does.
# Two steps are not shell commands and were done by solerita
# in Python: building the caps (Part 2.3) and retyping the
# backbone in the .ac file (Part 4). To run this script yourself
# you must do those two by hand, as the tutorial describes
# at the points marked below.
set -e

# Part 2: the CCD entry this was built from
wget https://files.rcsb.org/ligands/download/SME.cif

# Part 2.3: caps were added by solerita (Python), not by a shell command. The capped structure is SME.pdb
# SME.pdb

# Part 3: AM1-BCC charges and GAFF2 atom types
antechamber -i SME.pdb -fi pdb -o SME.mol2 -fo mol2 -at gaff2 -nc 0 -rn SME -c bcc -s 2

# Part 4: convert to .ac, the format prepgen reads
antechamber -i SME.mol2 -fi mol2 -o SME.ac -fo ac -at gaff2 -pf y

# Part 4: edit SME.ac so the backbone types read N H CX H1 C O (solerita did this; the unedited file is SME_gaff2_only.ac)
# SME.ac

# Part 6.3: the mainchain file
cat > SME.mc <<'EOF'
HEAD_NAME N
TAIL_NAME C
MAIN_CHAIN CA
OMIT_NAME CAT
OMIT_NAME CAY
OMIT_NAME CY
OMIT_NAME HT1
OMIT_NAME HT2
OMIT_NAME HT3
OMIT_NAME HY1
OMIT_NAME HY2
OMIT_NAME HY3
OMIT_NAME HZ
OMIT_NAME NZ
OMIT_NAME OY
PRE_HEAD_TYPE C
POST_TAIL_TYPE N
CHARGE 0.0
EOF

# Part 6.3: build the prepi
prepgen -i SME.ac -o SME.prepi -m SME.mc -rn SME

# Part 7.2, pass 1: everything parm10.dat can supply
parmchk2 -i SME.prepi -f prepi -o frcmod_parm10.SME -a Y -p /home/aalba/.claude-science/orgs/f90077d4-a9df-43a5-838b-4ea76472f1d5/workspaces/010d434a-d677-48c9-bd23-f9ab0db16732/nbtest/content/amber/dat/leap/parm/parm10.dat

# Part 7.2: keep only what parm10 really provides -- drop the ATTN placeholders and the zero-periodicity torsion lines parmchk2 writes next to them
grep -v ATTN frcmod_parm10.SME | awk '/^DIHE/{d=1} /^IMPROPER/{d=0} d && length($0)>11 {split(substr($0,12),f," "); if (f[4]!="" && f[4]+0==0) next} {print}' > frcmod1.SME

# Part 7.2, pass 2: GAFF2 fills in the rest
parmchk2 -i SME.prepi -f prepi -o frcmod2.SME -s 2

# Part 8: the test peptide
cat > peptide.leap <<'EOF'
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
EOF

# Part 8: build ALA-X-ALA
tleap -f peptide.leap

# Part 8: minimisation input
cat > min.in <<'EOF'
Minimisation
 &cntrl
  imin=1, maxcyc=500, ncyc=200, ntb=0, igb=0, cut=999.0,
 /
EOF

# Part 8: short minimisation
sander -O -i min.in -p SMEpeptide.prmtop -c SMEpeptide.ncrst -o min.out -r min.ncrst
