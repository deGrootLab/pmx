#!/usr/bin/env python
from pmx.model import Model
from pmx.forcefield import Topology
from pmx.alchemy import mutate, gen_hybrid_top
from pmx.gmx import set_gmxlib
from pmx.ligand_alchemy import _merge_FF_files
from filecmp import cmp

# set GMXLIB variable
set_gmxlib()

def test_mutate(gf, tmpdir):
    # load, mutate, and save the PDB file
    m = Model(gf('alchemy/protein.pdb'), rename_atoms=True)
    m2 = mutate(m=m, mut_resid=6, mut_resname='F', ff='amber99sb-star-ildn-mut')
    # reference output file
    ref_output = gf('alchemy/mutant.pdb')
    # save mutant and compare
    outfile = str(tmpdir.mkdir("test_alchemy").join("mutant.pdb"))
    m2.write(outfile)
    cmp(ref_output, outfile)


def test_gen_hybrid_top(gf, tmpdir):
    # load topology, fill B states, then write a new topology file
    topol = Topology(gf('alchemy/topol.top'))
    pmxtop, _ = gen_hybrid_top(topol)
    # reference output file
    ref_output = gf('alchemy/pmxtop.top')
    # save top and compare
    outfile = str(tmpdir.mkdir("test_alchemy").join("pmxtop.top"))
    pmxtop.write(outfile)
    cmp(ref_output, outfile)

    # test case where one topology includes itps containing hybrid residues
    # ---------------------------------------------------------------------
    topol = Topology(gf('alchemy/inc_topol.top'))
    # fill B states for hybrid residues present in topol2.itp
    pmxtop, pmxitps = gen_hybrid_top(topol=topol, recursive=True)
    # reference output file
    ref_top = gf('alchemy/inc_pmxtop.top')
    ref_itp = gf('alchemy/inc_pmxitp_0.itp')
    # save top and compare
    out_top = str(tmpdir.join("test_alchemy/pmx.top"))
    out_itp = str(tmpdir.join("test_alchemy/pmx.itp"))
    pmxtop.write(out_top)
    pmxitps[0].write(out_itp)
    # compare
    cmp(ref_top, out_top)
    cmp(ref_itp, out_itp)

def test_merging_itps(gf, tmpdir):
    #input files
    itp1=gf('alchemy/atomtypes/ffmol1.itp')
    itp2=gf('alchemy/atomtypes/ffmol2.itp')
    itp3=gf('alchemy/atomtypes/ffmol3.itp')
    #output file
    out_itp = str(tmpdir.join("ffmol_out.itp"))
    # reference output file
    ref_itp = gf('alchemy/atomtypes/ffmol_ref.itp')
    #run merge
    _merge_FF_files(out_itp, [itp1, itp2, itp3])
    #raise(Exception(out_itp))
    # compare
    assert cmp(ref_itp, out_itp), "does not match extected output."
