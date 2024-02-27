import sys, os
from pmx import *
from pmx.ndx import *

m = Model( "./protLig_benchmark/cdk2/ligands_gaff2/lig_1h1q/mol_gmx.pdb" )

import pytest
pytest.skip(allow_module_level=True)
# skip because don't have an example .ndx file

ndx = IndexFile("index.ndx")

print(ndx)
print(ndx['Backbone'])

atoms = ndx['Backbone'].select_atoms( m )
del ndx['Backbone']
grp = IndexGroup( "Backbone", atoms = atoms )

ndx.add_group( grp )



print(ndx)
#for atom in atoms:
#    print atom
