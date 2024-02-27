import sys, os
from pmx import *

from pmx.geometry import Rotation, Rotation2

import time

R = Rotation([0,0,0],[1,0,0])

m = Model("./protLig_benchmark/cdk2/protein_amber/protein.pdb")

t1 = time.time()

for i, r in enumerate(m.residues[1:297]):
    r.set_phi(-139, True)
    r.set_psi(135, True)
## for i in range(100):
##     for atom in m.atoms:
##         atom.x = R.apply( atom.x, 60*pi/180.)


t2 = time.time()
print(t2-t1)
m.write("out2.pdb")
