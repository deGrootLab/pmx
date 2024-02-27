from pmx import Molecule, Chain
from pmx.builder import build_dna_strand


def test_build():
    # start a peptide chain from an amino acid
    mol = Molecule().new_aa( "ALA" )      # build alanine
    c = Chain( residues = [ mol ] )       # initialize chain
    c.nbuild("T")                         # extend at n-terminus
    c.cbuild("R")                         # extend at c-terminus
    c.write('TAR.pdb')                    # write pdb file

    # build a peptide chain from sequence

    sequence = "FRTLKNCWQ"
    c = Chain().create( sequence )
    c.write("extended.pdb")

    for resi in c.residues:
        resi.set_phi( -57, True )
        resi.set_psi( -47, True )
    c.write("helix_pep.pdb")

    # fusing two chains
    seq = "LMNFRTS"
    c2 = Chain().create(seq)
    c.fuse( c2 )
    c.write("fused_pep.pdb")


    # create repeating seuqences



    rep_seq = "HEATLLFT"
    c = Chain().create( rep_seq )  # start with new chain
    new_chain = c.copy()

    for i in range(5):
        c.fuse( new_chain.copy() )
    c.write("rep_pep.pdb")
    print(c.sequence())
    helical_residues = c.fetch_residues(["HIS","GLU","ALA","THR","LEU"])
    for resi in helical_residues:
        resi.set_phi( -57, True )
        resi.set_psi( -47, True )
    c.write("repeat_helix.pdb")


    # building a DNA strand
    model = build_dna_strand("ACGTGTCA")
    model.write("dna.pdb")
