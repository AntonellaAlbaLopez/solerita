#!/usr/bin/env bash
# Every command solerita ran to parametrize CER (ligand, net charge +0).
# Read it next to the tutorial to see what each step does.
# Two steps are not shell commands and were done by solerita
# in Python: building the caps (Part 2.3) and retyping the
# backbone in the .ac file (Part 4). To run this script yourself
# you must do those two by hand, as the tutorial describes
# at the points marked below.
set -e

# Part 2: the CCD entry this was built from
wget https://files.rcsb.org/ligands/download/CER.cif

# Part 3: AM1-BCC charges and GAFF2 atom types
antechamber -i CER.pdb -fi pdb -o CER.mol2 -fo mol2 -at gaff2 -nc 0 -rn CER -c bcc -s 2

# Part 7.1: missing parameters, GAFF2
parmchk2 -i CER.mol2 -f mol2 -o CER.frcmod -s 2

# Part 5: the LEaP input
cat > lib.leap <<'EOF'
source leaprc.gaff2

CER = loadmol2 CER.mol2
loadamberparams CER.frcmod

charge CER
check CER

saveoff CER CER.lib
saveAmberParmNetCDF CER CER.prmtop CER.ncrst
savepdb CER CER_check.pdb
quit
EOF

# Part 5: build and save the library
tleap -f lib.leap

# Part 8: minimisation input
cat > min.in <<'EOF'
Minimisation
 &cntrl
  imin=1, maxcyc=500, ncyc=200, ntb=0, igb=0, cut=999.0,
 /
EOF

# Part 8: short minimisation
sander -O -i min.in -p CER.prmtop -c CER.ncrst -o min.out -r min.ncrst
