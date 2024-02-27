from pmx import *
from pmx.rotamer import *

c = Chain().create("ALARYTK")
print('done')
c.add_nterm_cap()
c.add_cterm_cap()

c.rename_atoms()
print(c.sequence())
c.write('x.pdb')
## bbdep = load_bbdep()

## res = c.residues[3]

## phi = res.get_phi(True)
## psi = res.get_psi(True)

## rot = get_rotamers( bbdep, "ARG", phi, psi)

## for i, r in enumerate(rot):
##     res.set_conformation( r )
##     c.write("rot%d.pdb" % i )
