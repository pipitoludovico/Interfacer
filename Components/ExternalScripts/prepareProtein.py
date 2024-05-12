#!/usr/bin/env python
#Giuseppe Deganutti @ Coventry University 2020
from htmd.ui import *
import sys

memb = sys.argv[2]

if memb == 'membrane':
	tick = 32.0

if memb == 'cytosol':
	tick = None

os.system("sed -i 's/HETATM/ATOM  /g' %s" %sys.argv[1])
prot = Molecule(sys.argv[1])

prot_H = systemPrepare(prot,  hydrophobic_thickness=tick, ignore_ns_errors=True, hold_nonpeptidic_bonds=True)
prot_H.write('protein_H.pdb')
os.system("sed -i 's/CYX/CYS/g' protein_H.pdb")
