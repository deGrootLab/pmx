# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
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
__doc__="""
Some functions to deal with tCNC data

"""

import sys, os
from .parser import *
from numpy import *
from . import library

def read_atom_types(f):
    if hasattr(f,"read"):
        fp = f
    else:
        fp = open(f,'r')
    l = fp.readlines()
    l = kickOutComments(l,'#')
    keys = []
    # we search for [ XXX ]
    for line in l:
        if line.startswith('['):
            entr = line.split()[1]
            keys.append(entr)
    # we turn each entry into a dic
    dic = {}
    for key in keys:
        dic[key] = {}
        sec = readSection(l,'[ '+key+' ]',"[")
        for line in sec:
            entr = line.split()
            if len(entr)==3:
                hyb = entr[2]
            else:
                hyb = ''
            pdic = {entr[0]:{'type':entr[1],
                             'hyb':hyb}}
            dic[key].update(pdic)
    return keys, dic

def make_lib_dic(f):
    keys, dic = read_atom_types(f)
    print('atom_types = {')
    for key in keys:
        val = dic[key]
        print("\t'%s': {" % key)
        for name, entr in list(val.items()):
            print("\t\t\"%s\" : {" % name)
            print("\t\t\t'type':'%s'," % entr['type'])
            print("\t\t\t'hyb':'%s'" % entr['hyb'])
            print("\t\t},")
        print("\t},")
    print("}")


def assign_types(model,verbose=False):
    atom_types = library._atom_types
    # default first
    dic = atom_types['DEFAULT']
    for atom in model.atoms:
        name = atom.long_name.strip()
        if name in dic:
            atom.atype = dic[name]['type']
            atom.hyb = dic[name]['hyb']
    #if isinstance(model,pmx.Model) or isinstance(model,pmx.Chain):
    if hasattr( model, "residues" ):
        for r in model.residues:
            key = r.resname
            if key in atom_types:
                dic = atom_types[key]
                for atom in r.atoms:
                    name = atom.long_name.strip()
                    if name in dic:
                        atom.atype = dic[name]['type']
                        atom.hyb = dic[name]['hyb']
    else:
    #elif isinstance(model,pmx.Molecule):
        key = model.resname
        if key in atom_types:
            dic = atom_types[key]
            for atom in model.atoms:
                name = atom.long_name.strip()
                if name in dic:
                    atom.atype = dic[name]['type']
                    atom.hyb = dic[name]['hyb']

    # check if we got all
    # and do generic
    dic = atom_types['GENERIC']
    for atom in model.atoms:
        if atom.atype == '':
            if atom.symbol in dic:
                atom.atype = dic[atom.symbol]['type']
                if verbose:
                    print('Using generic atom type for atom %d-%s/%d-%s (%s)' %\
                          (atom.id, atom.name, atom.resnr, atom.resname, atom.long_name))
            else:
                print('Could not assign atom type to atom %d-%s/%d-%s' %\
                          (atom.id, atom.name, atom.resnr, atom.resname))


def assign_radii(model):
    try:
        lst = open('Atomradii.dat').readlines()
    except:
        p = os.environ.get('CNCLIB')
        lst = open(os.path.join(p,'Atomradii.dat')).readlines()
    lst = kickOutComments(lst,';')
    tps = readSection(lst,'[ TYPES ]','[')
    tps = parseList('sff',tps)
    types = [a[0] for a in tps]
    dic = {}
    pdic = {}
    for i, line in enumerate(tps):
        dic[line[0]] = [line[1],line[2]]
        pdic[line[0]] = i
    for atom in model.atoms:
        atom.vdw = dic[atom.atype][0]
        atom.vdw14 = dic[atom.atype][1]
        atom.ptype = pdic[atom.atype]

    if hasattr(model,"chains"):
        comb = readSection(lst,'[ COMBINATIONS ]','[')
        comb = parseList('ssf',comb)
        comb14 = readSection(lst,'[ 14_COMBINATIONS ]','[')
        comb14 = parseList('ssf',comb14)

        size = len(tps)
        table = zeros((size,size))
        table14 = zeros((size,size))
        for i in range(size):
            for k in range(size):
                table[i][k] = tps[i][1]+tps[k][1]
                table14[i][k] = tps[i][2]+tps[k][2]
        for c in comb:
            idx1 = types.index(c[0])
            idx2 = types.index(c[1])
            table[idx1][idx2] = c[2]
            table[idx2][idx1] = c[2]

        for c in comb14:
            idx1 = types.index(c[0])
            idx2 = types.index(c[1])
            table14[idx1][idx2] = c[2]
            table14[idx2][idx1] = c[2]

        model.vdwtab = table
        model.vdw14tab = table14
