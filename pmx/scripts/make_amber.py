#!/usr/bin/env python
# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2011 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

# pmx script to make pdb files ready for use with amber
# dseelig@gwdg.de

import sys,os
from pmx import *


dna_res = ['DA','DG','DT','DC']
rna_res = ['RA','RU','RC','RG']

def chain_type(ch):
    if ch.residues[0].resname in list(library._one_letter.keys()):
        return 'pep'
    elif ch.residues[0].resname in dna_res:
        return 'dna'
    elif ch.residues[0].resname in rna_res:
        return 'rna'
    else: return 'unk'


# describe what the script does

help_text =("This script reads a structure file",
      "and outputs a new structure file",
      "with some renamed residues.",
      "The N-terminal and C-terminal residues",
      "are renamed such that we can use the",
      "AMBER force field. We also check for",
      "cysteine residues. So far the script works for"
      "proteins only. \n"
      "Take an all-atom protonated input pdb file (pdb2gmx -ff oplsaa)."
      "Renaming of residues (e.g. HIS) will be done based on"
      "the protons in the structure file.")

# define input/output files

files= [
    FileOption("-f","r",["pdb"],"protein.pdb","Input pdb file"),
    FileOption("-o","w",["pdb","gro"],"amber.pdb","Input pdb/gro file"),
]


# define options

options=[]

# pass options, files and the command line to pmx
cmdl = Commandline( sys.argv, options = options,
                    fileoptions = files,
                    program_desc = help_text,
                    check_for_existing_files = False )


# read structure file
model = Model(cmdl['-f'])


# now check cysteines
print('\nChecking cys....')
cysl = model.fetch_residues('CYS')     # get a list with all cysteines

# we do a simple check. If a HG is there it's CYS, else it's CYS2

for res in cysl:
    hg = res.fetch_atoms('HG')        # select HG atom from residue
    sg1 = res.fetch_atoms('SG')[0]
    if not hg:    # no hydrogen
        ss_bond = False
        for r in cysl:
            if r!=res:
                sg2 = r.fetch_atoms('SG')[0]
                d = sg1 - sg2
                if d < .25:
                    ss_bond = True
                    break
        if ss_bond:
            # terminal cys2 is ccyx
            rr = 'CYS2'
            print('Residue %d-%s (chain %s) will become %s' % (res.id, res.resname,res.chain_id,rr))
            res.set_resname(rr)
        else:
            print('Residue %d-%s (chain %s) will become %s' % (res.id, res.resname,res.chain_id,'CYM'))
            res.set_resname('CYM')

    else:
        res.set_resname('CYN')
        print('Residue %d-%s (chain %s) will become %s' % (res.id, res.resname, res.chain_id, res.resname))


# lysine
print('Checking lys....')
lysl = model.fetch_residues('LYS')
for res in lysl:
    at = res.fetch('HZ3')
    at2 = res.fetch('HZ2')
    if at or not at2:
        res.set_resname('LYP')
    print('Residue %d-%s (chain %s) will become %s' % (res.id, 'LYS', res.chain_id, res.resname))


# histidine
print('Checking his......')
hisl = model.fetch_residues('HIS')
for res in hisl:
    bHE2 = False
    bHD1 = False
    he2 = res.fetch('HE2')
    if he2: bHE2 = True
    hd1 = res.fetch('HD1')
    if hd1: bHD1 = True
    if hd1 and he2:
        res.set_resname('HIP')
    elif hd1 and not he2:
        res.set_resname('HID')
    elif he2 and not hd1:
        res.set_resname('HIE')
    else:
        res.set_resname('HID')
    print('Residue %d-%s (chain %s) will become %s' % (res.id, 'HIS', res.chain_id, res.resname))

print('Checking asp......')
aspl = model.fetch_residues('ASP')
for res in aspl:
    bHD2 = False
    hd2 = res.fetch('HD2')
    if hd2:
        res.set_resname('ASH')
    print('Residue %d-%s (chain %s) will become %s' % (res.id, 'ASP', res.chain_id, res.resname))

print('Checking glu......')
glul = model.fetch_residues('GLU')
for res in glul:
    bHD2 = False
    hd2 = res.fetch('HE2')
    if hd2:
        res.set_resname('GLH')
    print('Residue %d-%s (chain %s) will become %s' % (res.id, 'GLU', res.chain_id, res.resname))



print('Checking termini.....')
for chain in model.chains:
    print('Processing chain %s' % chain.id)
    ct = chain_type(chain)
    print('Chain type of chain %s: %s' % (chain.id, ct.upper()))
    first = chain.residues[0]      # first residue
    last = chain.residues[-1]      # last residue
    if ct == 'pep':
        if first.resname in list(library._one_letter.keys()):
            first.set_resname('N'+first.resname) # rename e.g. ALA to NALA
        if last.resname in list(library._one_letter.keys()):
            if last.resname == 'CYS2':
                last.set_resname('CCYX')   # rename e.g. ARG to CARG
            else:
                last.set_resname('C'+last.resname)   # rename e.g. ARG to CARG
        elif last.resname in library._ions: # find last protein residue
            print('Searching last peptide residue in chain %s' % chain.id)
            found = False
            idx = chain.residues.index(last)-1
            while not found:
                r = chain.residues[idx]
                if r.resname in list(library._one_letter.keys()):
                    if r.resname == 'CYS2':
                        r.set_resname('CCYX')
                    else:
                        r.set_resname('C'+r.resname)
                    last = r
                    found = True
                else:
                    idx-=1
        try:
            o1,o2 = last.fetchm(['O1','O2'])
            o1.name = 'OC1'
            o2.name = 'OC2'
        except:
            try:
                o1,o2 = last.fetchm(['O','OXT'])
                o1.name = 'OC1'
                o2.name = 'OC2'
            except:
                print('Error: No terminal O1, O2 atoms found in chain %s' % chain.id)
                print('       In pdb2gmx generated structures these should be there.')
                print('       Exiting')
                sys.exit(1)
    elif ct in ['dna','rna']:
        first.set_resname(first.resname+'5')
        last.set_resname(last.resname+'3')
        for r in chain.residues:
            try:
                o1p = r.fetch_atoms('OP1')[0]
                o1p.name = 'O1P'
                o2p = r.fetch_atoms('OP2')[0]
                o2p.name = 'O2P'
                o3p = r.fetch_atoms('OP3')[0]
                o3p.name = 'O3P'
            except:
                pass





# hack to get pdb file with 4 character residue names
#for atom in model.atoms:
#    atom.chain_id = ' '
#    if len(atom.resname)==4:
#        atom.chain_id = atom.resname[-1]+atom.chain_id
model.write(cmdl['-o'])
## fp = open(cmdl['-o'],'w')
## print 'Writing ouput to: %s' % cmdl['-o']
## for atom in model.atoms:
##     print >>fp, atom
