import sys, os
from pmx import *

c = Chain().create('AAAAAAAAAAAAAAAAAAAAA')
for r in c.residues:
    r.set_phi( -57, True )
    r.set_psi( -47, True )
c.write('hel.pdb')
