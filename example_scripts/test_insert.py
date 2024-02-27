import sys, os
from pmx import Model, Chain

def test_insert():
    c = Chain().create("ALKIRTS")
    m = Model("./protLig_benchmark/cdk2/protein_amber/protein.pdb")

    chB = m.chdic["B"]
    chB.insert_chain(2, c )
    m.write("ins.pdb")

    c = Chain().create("ALT")
    print(c.atoms[0].name)
    del c.atoms[0]                    # delete list item, first atom is gone
    print(c.atoms[0].name)                  # first atom is a different one now
    atom = c.residues[0].atoms[0]     # select first atom from the first residue
    print(atom.name)                        # should be the same as above
