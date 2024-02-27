import sys, os
from pmx import Model, Chain


def test_edit():
    seq1 = 'klrtsfcvnme'*2
    seq2 = 'irtiervcywq'*2


    c1 = Chain().create( seq1.upper() )
    c2 = Chain().create( seq2.upper() )
    c1.set_chain_id('A')
    c2.set_chain_id('B')
    c2.translate( [20,0,0] )
    m = Model(chains = [c1,c2] )
    m.write('protein.pdb')

    m = Model('protein.pdb')

    for atom in m.atoms:
        print(atom.id, atom.name, atom.resname)

    for resi in m.residues:
        if resi.resname in ['ALA','SER']:     # select some residues
            print(resi)
            for atom in resi.atoms:
                print(atom.bfac)                     # print some properties


    for c in m.chains:
        print(c.id, len(c.residues), len(c.atoms)) # print chain id, number of
        # residues and number of atoms

    chainB = m.chdic['B']            # select chain
    for atom in chainB.atoms:
        print(atom.occ)


    resl = m.fetch_residues(["LEU","PHE"])
    for r in resl:
        print('leu_or_phe:', r)

    resl = m.fetch_residues(["LEU","PHE"], inv = True)
    for r in resl:
        print('not leu_or_phe:', r)
