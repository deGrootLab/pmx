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
from .parser import *
import re
import sys
#---------------------------------------------------
class IndexGroup:
    """ CLASS TO DEAL WITH GROMACS INDEX FILES"""
    def __init__(self, name = '', ids = [], atoms = []):
        self.ids = []
        if atoms:
            self.ids = [a.id for a in atoms]
        else:
            self.ids = ids
        self.name = name
        self.ids.sort()

    def __str__(self):
        s = ''
        s+= '[ %s ]\n' % self.name
        count = 0
        for i in range(len(self.ids)):
            if len(str(self.ids[i])) > 6:
                fs = "%d "
            else:
                fs = "%6d "
            s+=fs % self.ids[i]
            count+=1
            if count % 15 == 0: s+='\n'
        return s

    def read_index_group(self,name,lines):
        if '[' in name:
            name = name[1:-1].strip()
        self.name = name
        key = '[ '+name+' ]'
        ll = readSection(lines,key,'[')
        self.ids = []
        for line in ll:
            self.ids.extend([int(x) for x in line.split()])
        return self

    def select_atoms( self, atom_sel):
        atoms = []
        for idx in self.ids:
            for atom in atom_sel.atoms:
                if atom.id == idx:
                    atoms.append( atom )
        return atoms

#---------------------------------------------------

class IndexFile:

    def __init__(self, fname = None, names = [], groups = []):
        self.groups = []
        self.names = []
        self.dic = {}
        if fname is not None:
            self.parse(fname)
        if groups:
            for g in groups:
                self.add_group( g )

    def parse(self,fp):
        if hasattr(fp,"read"):
            f = fp.read()
        else:
            f = open(fp).read()
        names = self.__get_index_groups(f)
        lines = f.split('\n')
        for name in names:
            idx = IndexGroup().read_index_group(name,lines)
            self.add_group( idx )
        return self

    def __get_index_groups(self,f):
        x = re.compile('\[ .* \]\n')
        names = x.findall(f)
        r = []
        for name in names:
            r.append(name.strip())
        return r

    def write(self,fn=None,fp=None):
        if not fn and not fp:
            fp = sys.stdout
        if fn:
            fp = open(fn,'w')
        for gr in self.groups:
            print(str(gr)+'\n', file=fp)


    def __str__(self):
        s = 'Gromacs index file (%d index groups):\n' % len(self.groups)
        for i, group in enumerate( self.groups ):
            s+=' %d  %-20s   :    %d atoms\n' % (i, group.name, len(group.ids) )
        return s

    def __getitem__(self, item):
        return self.dic[item]

    def add_group( self, group ):
        if group.name in self.names:
            print("IndexFile has group %s !! " % group.name, file=sys.stderr)
            print("Group %s will be replaced !!" % group.name, file=sys.stderr)
            self.delete_group( group.name )
        self.names.append( group.name )
        self.groups.append( group )
        self.dic[group.name] = group

    def delete_group( self, name ):
        idx = -1
        for i, group in enumerate(self.groups):
            if group.name == name:
                idx = i
                break
        if idx == -1: return
        del self.groups[idx]
        del self.dic[name]
        self.names.remove( name )

    def __delitem__(self, item ):
        self.delete_group( item )



#---------------------------------------------------


def get_index(atom_list = None, residue_list = None, chain_list = None):
    """ return atom indices from a list of atoms/residues/chains"""
    if not atom_list and not residue_list and not chain_list:
        print('Error: Need list~')
        sys.exit(1)
    if atom_list:
        lst = [a.id for a in atom_list]
        return lst
    if residue_list:
        al = []
        list(map(lambda r: al.extend(r.atoms), residue_list))
        return get_index(atom_list = al)
    if chain_list:
        al = []
        list(map(lambda c: al.extend(c.atoms), chain_list))
        return get_index(atom_list = al)


def make_index_group(atomlist, name):
    lst = get_index(atomlist)
    g = IndexGroup(ids = lst, name = name)
    return g
