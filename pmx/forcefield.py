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
Functions to read gromacs forcefield files
"""
import sys,os,re, copy
from .parser import *
from . import cpp
from .atom import Atom
from .odict import *
from .library import _aliases
from .ffparser import *
from  . import _pmx as _p

def get_bond_param(type1,type2,bond_lib):
    for entr in bond_lib:
        if (type1==entr[0] and type2==entr[1]) or \
           (type2==entr[0] and type1==entr[1]):
            return entr[2:]
    return None


def get_angle_param(type1,type2,type3,ang_lib):

    for entr in ang_lib:
        if (type1 == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2]) or \
           (type1 == entr[2] and \
            type2 == entr[1] and \
            type3 == entr[0]):
            return entr[3:]
    return None

def get_dihedral_param(type1,type2,type3,type4,dih_lib, func):
    for entr in dih_lib:
        if (type1 == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2] and \
            type4 == entr[3] and func==entr[4]) or \
           (type1 == entr[3] and \
            type2 == entr[2] and \
            type3 == entr[1] and \
            type4 == entr[0] and func==entr[4]):
            return entr[4:]
    for entr in dih_lib:
        if ('X' == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2] and \
            type4 == entr[3] and func==entr[4]) or \
           (type1 == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2] and \
            'X' == entr[3] and func==entr[4]):
            return entr[4:]
        if ('X' == entr[3] and \
            type2 == entr[2] and \
            type3 == entr[1] and \
            type4 == entr[0] and func==entr[4]) or \
           (type1 == entr[3] and \
            type2 == entr[2] and \
            type3 == entr[1] and \
            'X' == entr[0] and func==entr[4]):
            return entr[4:]
    for entr in dih_lib:
        if ('X' == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2] and \
            'X' == entr[3] and func==entr[4]) or \
           ('X' == entr[3] and \
            type2 == entr[2] and \
            type3 == entr[1] and \
            'X' == entr[0] and func==entr[4]):
            return entr[4:]
    for entr in dih_lib:
        if ('X' == entr[0] and \
            'X' == entr[1] and \
            type3 == entr[2] and \
            type4 == entr[3] and func==entr[4]) or \
           (type1 == entr[3] and \
            type2 == entr[2] and \
            'X' == entr[1] and \
            'X' == entr[0] and func==entr[4]):
            return entr[4:]
    for entr in dih_lib:
        if (typ1 == entr[0] and \
            'X' == entr[1] and \
            'X' == entr[2] and \
            type4 == entr[3] and func==entr[4]) or \
           (type1 == entr[3] and \
            'X' == entr[2] and \
            'X' == entr[1] and \
            type4 == entr[0] and func==entr[4]):
            return entr[4:]
    return None



class ITPFile:

    def __init__(self, fname = None, ff = 'amber99sb'):
        self.name = 'MOL'
        self.nrexcl = 3
        self.atoms = []
        self.pairs = []
        self.angles = []
        self.dihedrals = []
        self.atomtypes = None
        self.virtual_sites1 = []
        self.virtual_sites2 = []
        self.virtual_sites3 = []
        self.virtual_sites4 = []
        self.has_vsites1 = False
        self.has_vsites2 = False
        self.has_vsites3 = False
        self.has_vsites4 = False
        if fname:
            self.read(fname, ff = ff)

    def read(self,fname, ff = None):
        if not hasattr(fname,"readlines"):
            lines = open(fname).readlines()
        else:
            lines = fname.readlines()
        lines = kickOutComments(lines,';')
        lines = kickOutComments(lines,'#')
        self.name, self.nrexcl = read_moleculetype(lines)
        self.atoms = read_itp_atoms(lines)
        self.bonds = read_itp_bonds(lines)
        self.pairs = read_itp_pairs(lines)
        self.angles = read_itp_angles(lines)
        self.dihedrals = read_itp_dihedrals(lines)
        self.atomtypes = read_atomtypes(lines,ff)
        self.read_vsites2(lines)


    def write(self,fname):
        if not hasattr(fname,"write"):
            fp = open(fname,"w")
        else:
            fp = fname
        write_itp_moleculetype(fp,self.name,self.nrexcl)
        write_itp_atoms(fp, self.atoms)
        if self.bonds:
            write_itp_bonds(fp, self.bonds)
        if self.pairs:
            write_itp_pairs(fp, self.pairs)
        if self.angles:
            write_itp_angles(fp, self.angles)
        if self.dihedrals:
            write_itp_dihedrals(fp, self.dihedrals)
        if self.has_vsites2:
            self.write_itp_vsites2(fp)

    def write_itp_vsites2(self, fp ):
        print('[ virtual_sites2 ]', file=fp)
        for v in self.virtual_sites2:
            print("%8d %8d %8d %s %s" % (v[0].id, v[1].id, v[2].id, v[3], v[4]), file=fp)

    def set_name(self, name):
        self.name = name
        for atom in self.atoms:
            atom.resname = name

    def as_rtp(self):
        for i, bond in enumerate(self.bonds):
            id1 = bond[0]
            id2 = bond[1]
            self.bonds[i][0] = self.atoms[id1-1]
            self.bonds[i][1] = self.atoms[id2-1]
        for i, angle in enumerate(self.angles):
            id1 = angle[0]
            id2 = angle[1]
            id3 = angle[2]
            self.angles[i][0] = self.atoms[id1-1]
            self.angles[i][1] = self.atoms[id2-1]
            self.angles[i][2] = self.atoms[id3-1]

        for i, dih in enumerate(self.dihedrals):
            id1 = dih[0]
            id2 = dih[1]
            id3 = dih[2]
            id4 = dih[3]
            self.dihedrals[i][0] = self.atoms[id1-1]
            self.dihedrals[i][1] = self.atoms[id2-1]
            self.dihedrals[i][2] = self.atoms[id3-1]
            self.dihedrals[i][3] = self.atoms[id4-1]

        for i, vs in enumerate(self.virtual_sites2):
            id1 = dih[0]
            id2 = dih[1]
            id3 = dih[2]
            self.virtual_sites2[i][0] = self.atoms[id1-1]
            self.virtual_sites2[i][1] = self.atoms[id2-1]
            self.virtual_sites2[i][2] = self.atoms[id3-1]



    def id2atoms(self):

        for i, bond in enumerate(self.bonds):
            id1 = bond[0]
            id2 = bond[1]
            self.bonds[i][0] = self.atoms[id1-1]
            self.bonds[i][1] = self.atoms[id2-1]

        for i, pairs in enumerate(self.pairs):
            id1 = pairs[0]
            id2 = pairs[1]
            self.pairs[i][0] = self.atoms[id1-1]
            self.pairs[i][1] = self.atoms[id2-1]

        for i, angle in enumerate(self.angles):
            id1 = angle[0]
            id2 = angle[1]
            id3 = angle[2]
            self.angles[i][0] = self.atoms[id1-1]
            self.angles[i][1] = self.atoms[id2-1]
            self.angles[i][2] = self.atoms[id3-1]

        for i, dih in enumerate(self.dihedrals):
            id1 = dih[0]
            id2 = dih[1]
            id3 = dih[2]
            id4 = dih[3]
            self.dihedrals[i][0] = self.atoms[id1-1]
            self.dihedrals[i][1] = self.atoms[id2-1]
            self.dihedrals[i][2] = self.atoms[id3-1]
            self.dihedrals[i][3] = self.atoms[id4-1]



    def write_rtp(self, filename ='mol.rtp'):
        fp = open(filename,'w')
        print('[ %s ]' % self.name, file=fp)
        print(' [ atoms ]', file=fp)
        for atom in self.atoms:
            print("%8s %-12s %8.6f %5d" % \
                  (atom.name,atom.atomtype,atom.q,atom.cgnr), file=fp)

        print('\n [ bonds ]', file=fp)
        for bond in self.bonds:
            if len(bond)<=3:
                print("%8s %8s "% \
                      (bond[0].name, bond[1].name), file=fp)
            else:
                print("%8s %8s %8.4f %8.4f "% \
                      (bond[0].name, bond[1].name, bond[3], bond[4]), file=fp)


        print('\n [ angles ]', file=fp)
        for angle in self.angles:
            if len(angle)<=4:
                print("%8s %8s %8s "% \
                      (angle[0].name, angle[1].name,angle[2].name), file=fp)
            elif angle[3]==5: # U-B
                print("%8s %8s %8s %8.4f %8.4f %8.4f %8.4f "% \
                      (angle[0].name, angle[1].name,angle[2].name,
                        angle[4],angle[5],angle[6],angle[7]), file=fp)
            else:
                print("%8s %8s %8s %8.4f %8.4f "% \
                      (angle[0].name, angle[1].name,angle[2].name,angle[4],angle[5]), file=fp)


        print('\n [ dihedrals ]', file=fp)
        for dih in self.dihedrals:
            if len(dih)<=5: # no parameters
                print("%8s %8s %8s %s "% \
                      (dih[0].name, dih[1].name,dih[2].name, dih[3].name), file=fp)
            elif dih[4]==3:
                print("%8s %8s %8s %s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f "% \
                      (dih[0].name, dih[1].name,dih[2].name, dih[3].name,
                       dih[5], dih[6], dih[7], dih[8], dih[9], dih[10]), file=fp)
            elif (dih[4]==1) or (dih[4]==4) or (dih[4]==9):
                print("%8s %8s %8s %s %8.4f %8.4f %8.4f "% \
                      (dih[0].name, dih[1].name,dih[2].name, dih[3].name,
                       dih[5], dih[6], dih[7]), file=fp)
            elif (dih[4]==2) or (dih[4]==11):
                print("%8s %8s %8s %s %8.4f %8.4f "% \
                      (dih[0].name, dih[1].name,dih[2].name, dih[3].name,
                       dih[5], dih[6]), file=fp)

#        print >>fp, '\n [ impropers ]'
#        for dih in self.impropers:
#            if len(dih)<=5: # no parameters
#                print >>fp, "%8s %8s %8s %s "% \
#                      (dih[0].name, dih[1].name,dih[2].name, dih[3].name)
#            elif dih[4]==2:
#                print >>fp, "%8s %8s %8s %s %8.4f %8.4f "% \
#                      (dih[0].name, dih[1].name,dih[2].name, dih[3].name, dih[5], dih[6])
#            elif dih[4]==4:
#                print >>fp, "%8s %8s %8s %s %8.4f %8.4f %8.4f "% \
#                      (dih[0].name, dih[1].name,dih[2].name, dih[3].name, dih[5], dih[6], dih[7])


    def read_vsites2(self, lines):
        starts = []
        dih = []
        for i, line in enumerate(lines):
            if line.strip().startswith('[ virtual_sites2 ]'):
                starts.append(i)
        if starts:
            self.has_vsites2 = True
        for s in starts:
            lst = readSection(lines[s:],'[ virtual_sites2 ]','[')
            for line in lst:
                entr = line.split()
                idx = [int(x) for x in entr[:3]]

                func = int(entr[3])
                try:
                    rest = ' '.join(entr[4:])
                except:
                    rest = ''
                self.virtual_sites2.append([self.atoms[idx[0]-1],\
                                            self.atoms[idx[1]-1],\
                                            self.atoms[idx[2]-1],\
                                            func,rest])


    def read_vsites3(self, lines):
        starts = []
        dih = []
        for i, line in enumerate(lines):
            if line.strip().startswith('[ virtual_sites3 ]'):
                starts.append(i)
        if starts:
            self.has_vsites3 = True
        for s in starts:
            lst = readSection(lines[s:],'[ virtual_sites3 ]','[')
            for line in lst:
                entr = line.split()
                idx = [int(x) for x in entr[:4]]

                func = int(entr[4])
                try:
                    rest = ' '.join(entr[5:])
                except:
                    rest = ''
                self.virtual_sites3.append([self.atoms[idx[0]-1],\
                                            self.atoms[idx[1]-1],\
                                            self.atoms[idx[2]-1],\
                                            self.atoms[idx[3]-1],\
                                            func,rest])

    def read_vsites4(self, lines):
        starts = []
        dih = []
        for i, line in enumerate(lines):
            if line.strip().startswith('[ virtual_sites4 ]'):
                starts.append(i)
        if starts:
            self.has_vsites4 = True
        for s in starts:
            lst = readSection(lines[s:],'[ virtual_sites4 ]','[')
            for line in lst:
                entr = line.split()
                idx = [int(x) for x in entr[:5]]

                func = int(entr[5])
                try:
                    rest = ' '.join(entr[6:])
                except:
                    rest = ''
                self.virtual_sites4.append([self.atoms[idx[0]-1],\
                                            self.atoms[idx[1]-1],\
                                            self.atoms[idx[2]-1],\
                                            self.atoms[idx[3]-1],\
                                            self.atoms[idx[4]-1],\
                                            func,rest])


class Topology:

    def __init__(self, filename = None, ff = 'amber99sb', itp=False, top = None):
        self.filename = filename
        self.is_itp = itp
        if self.is_itp and top == None:
            print("Error:You have to provide the .top file if you read a .itp")
            sys.exit(1)
        if top is not None:
            self.topfile = top
        else:
            self.topfile = self.filename
        self.atoms = []
        self.bonds = []
        self.constraints = []
        self.have_constraints = False
        self.pairs = []
        self.angles = []
        self.dihedrals = []
        self.virtual_sites3 = []
        self.virtual_sites4 = []
        self.has_vsites3 = False
        self.has_vsites4 = False
        self.molecules = []
        self.system = ''
        self.qA = 0.
        self.qB = 0.
        if filename is not None:
            self.read_top(filename, ff = ff)
            l = cpp_parse_file(self.topfile)
            l = kickOutComments(l,'#')
            l = kickOutComments(l,';')
            self.BondedParams = BondedParser( l )
            self.NBParams = NBParser( l )
#            self.types, self.bond_lib, self.ang_lib, self.dih_lib = \
#                        read_ff(self.topfile,ff=ff)
            self.assign_fftypes()

    def read_top(self, fname, ff = 'amber99sb'):
        if not hasattr(fname,"readlines"):
            lines = open(fname).readlines()
        else:
            lines = fname.readlines()
        lines = kickOutComments(lines,';')
        if not self.is_itp:
            self.read_header(lines)
            self.read_footer(lines)
        lines = kickOutComments(lines,'#')
        self.read_moleculetype(lines)
        self.read_atoms(lines)
        if not self.atoms:
            self.no_itp = False
        else:
            self.read_bonds(lines)
            self.read_constraints(lines)
            self.read_pairs(lines)
            self.read_angles(lines)
            self.read_dihedrals(lines)
            self.read_vsites3(lines)
            self.read_vsites4(lines)
            self.no_itp = True
        if not self.is_itp:
            self.read_system(lines)
            self.read_molecules(lines)

    def assign_forcefield_parameters(self, eval_cpp = True):
        if eval_cpp:
            proc = cpp.PreProcessor()
            proc(self.filename)
            self.cpp_dic = proc.cpp_namespace
            for d in self.dihedrals:
                if len(d) == 6:
                    if not hasattr(d[5],"append"):
                        if d[5] in self.cpp_dic:
                            d[5] = [float(x) for x in self.cpp_dic[d[5]].split()]
                elif len(d) == 7:
                    if not hasattr(d[5],"append"):
                        if d[5] in self.cpp_dic:
                            d[5] = [float(x) for x in self.cpp_dic[d[5]].split()]
                    if not hasattr(d[6],"append"):
                        if d[6] in self.cpp_dic:
                            d[6] = [float(x) for x in self.cpp_dic[d[6]].split()]
        self.make_bond_params()
        self.make_angle_params()
        self.make_dihedral_params()


    def get_qA(self):
        qA = 0
        for atom in self.atoms:
            qA+=atom.q
        return round(qA,3)

    def get_qB(self):
        qB = 0
        for atom in self.atoms:
            if atom.atomtypeB is not None:
                qB+=atom.qB
            else:
                qB+=atom.q
        return round(qB,3)

    def make_Bstates(self, subset = None):
        if subset:
            atomlist = subset
        else:
            atomlist = self.atoms
        for atom in atomlist:
            atom.atomtypeB = atom.atomtype
            atom.mB = atom.m
            atom.qB = atom.q

        for i, b in enumerate(self.bonds):
            if len(b) > 3:
                if self.bonds[i][0] in atomlist or \
                       self.bonds[i][1] in atomlist:
                    param = copy.deepcopy(b[-1])
                    self.bonds[i].append(param)
        for i, ang in enumerate(self.angles):
            if len(ang) > 4:
                if self.angles[i][0] in atomlist or \
                       self.angles[i][1] in atomlist or \
                       self.angles[i][2] in atomlist:
                    param = copy.deepcopy(ang[-1])
                    self.angles[i].append(param)
        for i, dih in enumerate(self.dihedrals):
            if len(dih) > 5:
                if self.dihedrals[i][0] in atomlist or \
                   self.dihedrals[i][1] in atomlist or \
                   self.dihedrals[i][2] in atomlist or \
                   self.dihedrals[i][3] in atomlist:
                    param = copy.deepcopy(dih[-1])
                    self.dihedrals[i].append(param)

    def set_Astate_zero(self):
        for i, b in enumerate(self.bonds):
            for k in range(len(self.bonds[i][-2])):
                self.bonds[i][-2][k] = 0
        for i, ang in enumerate(self.angles):
            for k in range(len(self.angles[i][-2])):
                self.angles[i][-2][k] = 0
        for i, dih in enumerate(self.dihedrals):
            if dih[4] == 3:
                for k in range(len(self.dihedrals[i][-2])):
                    self.dihedrals[i][-2][k] = 0
            elif dih[4] == 1:
                self.dihedrals[i][-2][1] = 0

    def read_molecules(self,lines):
        lst = readSection(lines,'[ molecules ]','[')
        self.molecules = []
        for line in lst:
            entr = line.split()
            self.molecules.append([entr[0],int(entr[1])])

    def read_system(self,lines):
        lst = readSection(lines,'[ system ]','[')
        self.system = lst[0].strip()



    def read_atoms(self,lines):
        lst = readSection(lines,'[ atoms ]','[')
        self.atoms = []
        for line in lst:
            a = topline2atom(line)
            self.atoms.append(a)

    def read_bonds(self,lines):
        lst = readSection(lines,'[ bonds ]','[')
        self.bonds = []
        for line in lst:
            entries = line.split()
            if len(entries) == 3:
                idx = [int(x) for x in line.split()]
                self.bonds.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], idx[2]])
            elif len(entries) == 5:
                idx = [int(x) for x in entries[:3]]
                l = float(entries[3])
                k = float(entries[4])
                self.bonds.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], idx[2], [l,k]])

            elif len(entries) == 7:
                idx = [int(x) for x in entries[:3]]
                lA = float(entries[3])
                kA = float(entries[4])
                lB = float(entries[5])
                kB = float(entries[6])
                self.bonds.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], idx[2], [lA,kA],[lB,kB]])

    def read_pairs(self,lines):
        lst = readSection(lines,'[ pairs ]','[')
        self.pairs = []
        for line in lst:
            idx = [int(x) for x in line.split()]
            self.pairs.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], idx[2]])

    def read_constraints(self,lines):
        lst = readSection(lines,'[ constraints ]','[')
        self.constraints = []
        for line in lst:
            idx = [int(x) for x in line.split()]
            self.constraints.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], idx[2]])
        if self.constraints:
            self.have_constraints = True

    def read_angles(self, lines):
        lst = readSection(lines,'[ angles ]','[')
        angles = []
        for line in lst:
            entries = line.split()
            if len(entries) == 4:
                idx = [int(x) for x in line.split()]
                self.angles.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], \
                                    self.atoms[idx[2]-1], idx[3]])
            elif len(entries) == 6:
                idx = [int(x) for x in entries[:4]]
                l = float(entries[4])
                k = float(entries[5])
                self.angles.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], \
                                    self.atoms[idx[2]-1], idx[3], [l,k]])
            elif len(entries) == 8:
                idx = [int(x) for x in entries[:4]]
                lA = float(entries[4])
                kA = float(entries[5])
                lB = float(entries[6])
                kB = float(entries[7])
                self.angles.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], \
                                    self.atoms[idx[2]-1], idx[3], [lA,kA],[lB,kB]])

    def read_dihedrals(self, lines):
        starts = []
        dih = []
        for i, line in enumerate(lines):
            if line.strip().startswith('[ dihedrals ]'):
                starts.append(i)
        for s in starts:
            lst = readSection(lines[s:],'[ dihedrals ]','[')
            for line in lst:
                entr = line.split()
                idx = [int(x) for x in entr[:4]]

                func = int(entr[4])
                try:
                    rest = ' '.join(entr[5:])
                except:
                    rest = ''
                self.dihedrals.append([self.atoms[idx[0]-1],\
                                       self.atoms[idx[1]-1],\
                                       self.atoms[idx[2]-1],\
                                       self.atoms[idx[3]-1],\
                                       func,rest])
    def read_vsites3(self, lines):
        starts = []
        dih = []
        for i, line in enumerate(lines):
            if line.strip().startswith('[ virtual_sites3 ]'):
                starts.append(i)
        if starts:
            self.has_vsites3 = True
        for s in starts:
            lst = readSection(lines[s:],'[ virtual_sites3 ]','[')
            for line in lst:
                entr = line.split()
                idx = [int(x) for x in entr[:4]]

                func = int(entr[4])
                try:
                    rest = ' '.join(entr[5:])
                except:
                    rest = ''
                self.virtual_sites3.append([self.atoms[idx[0]-1],\
                                            self.atoms[idx[1]-1],\
                                            self.atoms[idx[2]-1],\
                                            self.atoms[idx[3]-1],\
                                            func,rest])

    def read_vsites4(self, lines):
        starts = []
        dih = []
        for i, line in enumerate(lines):
            if line.strip().startswith('[ virtual_sites4 ]'):
                starts.append(i)
        if starts:
            self.has_vsites4 = True
        for s in starts:
            lst = readSection(lines[s:],'[ virtual_sites4 ]','[')
            for line in lst:
                entr = line.split()
                idx = [int(x) for x in entr[:5]]

                func = int(entr[5])
                try:
                    rest = ' '.join(entr[6:])
                except:
                    rest = ''
                self.virtual_sites4.append([self.atoms[idx[0]-1],\
                                            self.atoms[idx[1]-1],\
                                            self.atoms[idx[2]-1],\
                                            self.atoms[idx[3]-1],\
                                            self.atoms[idx[4]-1],\
                                            func,rest])

    def read_header(self, lines):
        ret = []
        for line in lines:
            if not line.strip().startswith('[') and \
                   not line.strip().startswith('#ifdef POSRES'):
                ret.append(line.rstrip())
            else:
                break
        self.header = ret

    def read_footer(self, lines):
        for line in lines:
            if line.strip().startswith('#ifdef POSRES'):
                idx = lines.index(line)
                self.footer =  [l.rstrip() for l in lines[idx:]]
                break
        idx = self.footer.index('[ system ]')
        self.footer = self.footer[:idx]

    def read_moleculetype(self, lines):
        l = readSection(lines,'[ moleculetype ]','[')
        if l:
            self.name, self.nexcl =  l[0].split()[0], int(l[0].split()[1])

    def set_molecule(self, molname, n):
        mol_exists = False
        for i, mol in enumerate(self.molecules):
            if mol[0] == molname:
                self.molecules[i][1] = n
                mol_exists = True
        if not mol_exists:
            self.molecules.append([molname,n])

    def del_molecule(self, molname):
        if not hasattr(molname,"append"):
            molname = [molname]
        new = []
        for m in self.molecules:
            if m[0] not in molname:
                new.append(m)
        self.molecules = new

    def write_top(self, outfile, stateBonded = 'AB', stateTypes = 'AB', stateQ = 'AB',
                  scale_mass = False, dummy_qA = 'on', dummy_qB = 'on', target_qB = [],
                  full_morphe = True):

        fp = open(outfile,'w')
        if not self.is_itp:
            self.write_header(fp)
        if self.no_itp:
            self.write_moleculetype(fp)
            self.write_atoms(fp, charges = stateQ, atomtypes = stateTypes, dummy_qA = dummy_qA,
                             dummy_qB = dummy_qB, scale_mass = scale_mass,
                             target_qB = target_qB, full_morphe = full_morphe)
            self.write_bonds(fp, state = stateBonded)
            if self.have_constraints:
                self.write_constraints(fp)
            self.write_pairs(fp)
            self.write_angles(fp, state = stateBonded)
            self.write_dihedrals(fp, state = stateBonded)
            if self.has_vsites3:
                self.write_vsites3(fp)
            if self.has_vsites4:
                self.write_vsites4(fp)
        if not self.is_itp:
            self.write_footer(fp)
            self.write_system(fp)
            self.write_molecules(fp)
        fp.close()

    def write_header(self,fp):
        for line in self.header:
            print(line, file=fp)

    def write_footer(self,fp):
        for line in self.footer:
            print(line, file=fp)

        def write_moleculetype(self, fp):
            print('[ moleculetype ]', file=fp)
            print('; Name        nrexcl', file=fp)
            print('%s  %d' % (self.name,self.nrexcl), file=fp)

    def __atoms_morphe( self, atoms ):
        for atom in atoms:
            if atom.atomtypeB is not None and (atom.q!=atom.qB or atom.m != atom.mB): return True
        return False

    def __atomtypes_morphe(self, atoms):
        for atom in atoms:
            if atom.atomtypeB is not None and atom.atomtype != atom.atomtypeB: return True
        return False

    def __is_perturbed_residue( self, residue ):
        if self.__atoms_morphe(residue.atoms): return True
        return False

    def __last_perturbed_atom(self, r):

        max_order = 0
        last_atom = None
        for atom in r.atoms:
            if self.__atoms_morphe([atom]) and atom.name not in ['N','CA','C','O','H']:
                if not atom.atomtype.startswith('DUM') and not atom.atomtypeB.startswith('DUM'):
                    last_atom = atom
        if last_atom == None:
            print('Error: Could not find a perturbed atom to put rest charges on !', file=sys.stderr)
            sys.exit(1)
        return last_atom

    def check_special_dihedrals( self ):
        for d in self.dihedrals:
            A, B = self.check_case( d[:4] )
            if ('D' not in A and 'D' in B) or ('D' in A and 'D' not in B):
                print(d[0].name, d[1].name, d[2].name, d[3].name, d[4:], A, B)

    def write_atoms(self, fp, charges = 'AB', atomtypes = 'AB', dummy_qA = 'on',\
                dummy_qB = 'on',  scale_mass=True, target_qB = [], full_morphe = True):

        self.qA = 0
        self.qB = 0
        for r in self.residues:
            if self.__is_perturbed_residue(r):
                target_chargeB = target_qB.pop(0)
                print('Making target charge %g for residue %s' % (round(target_chargeB,5), r.resname))
                for atom in r.atoms:
                    if self.__atoms_morphe([atom]):
                        if charges == 'AB':      # we move the charges from state A to state B
                            atom.qqA = atom.q
                            atom.qqB = atom.qB
                            if not full_morphe and (atom.q*atom.qB < 0 or atom.atomtype!=atom.atomtypeB):   # we change a charge from + to - or vice versa
                                atom.qqB = 0
                                atom.to_be_morphed = True
                            else:
                                atom.qqB = atom.qB
                        elif charges == 'AA':        # we keep the charges
                            atom.qqA = atom.q
                            atom.qqB = atom.q

                        elif charges == 'BB':        # take charges of state B
                            if not full_morphe:
                                if hasattr(atom,"contQ"):
                                    atom.qqA = atom.contQ
                                    atom.qqB = atom.qqA
                                if hasattr(atom,"to_be_morphed"): # this a big q morphe. has been set to zero before
                                    if atomtypes == 'BB':
                                        atom.qqA = 0
                                        atom.qqB = atom.qB
                                    elif atomtypes == 'AB':
                                        atom.qqA = 0
                                        atom.qqB = 0
                                elif not hasattr(atom,"contQ") and not hasattr(atom,"to_be_morphed") :
                                    atom.qqA = atom.qB
                                    atom.qqB = atom.qB
                            else:
                                atom.qqA = atom.qB
                                atom.qqB = atom.qB
                        if atom.atomtype.startswith('DUM') or atom.atomtypeB.startswith('DUM'):
                            if dummy_qA == 'off':
                                atom.qqA = 0.
                            if dummy_qB == 'off':
                                atom.qqB = 0.
                    else:
                        atom.qqA = atom.q
                        atom.qqB = atom.q
                qA_tot = sum([a.qqA for a in r.atoms])
                qB_tot = sum([a.qqB for a in r.atoms])
                if qB_tot != target_chargeB:
                    print('State B has total charge of %g' % round(qB_tot,5))
                    print('Applying charge correction to ensure integer charges')
                    latom = self.__last_perturbed_atom(r)
                    print('Selecting atom %d-%s (%s) as perturbed atom with highest order' % (latom.id,latom.name, latom.resname))
                    newqB = latom.qqB-(qB_tot-target_chargeB)
                    print('Changing chargeB of atom %s from %g to %g' % (latom.name, latom.qqB,newqB))
                    latom.qqB = newqB
                    qB_tot = sum([a.qqB for a in r.atoms])
                    print('New total charge of B-state is %g' % round(qB_tot,5))
                else:
                    print('No corrections applied to ensure integer charges')


        print('\n [ atoms ]', file=fp)
        print(';   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB', file=fp)
        al = self.atoms
        for atom in al:
            if self.__atoms_morphe([atom]):

                if atomtypes == 'AB':
                    atA = atom.atomtype
                    atB = atom.atomtypeB
                    mA = atom.m
                    mB = atom.mB
                elif atomtypes == 'AA':
                    atA = atom.atomtype
                    atB = atom.atomtype
                    mA = atom.m
                    mB = atom.m
                elif atomtypes == 'BB':
                    atA = atom.atomtypeB
                    atB = atom.atomtypeB
                    mA = atom.mB
                    mB = atom.mB
                if scale_mass:
                    if atA.startswith('DUM'):
                        mA = 1.
                    if atB.startswith('DUM'):
                        mB = 1.
                if hasattr(atom,"qqB"):
                    qqB = atom.qqB
                    if hasattr(atom,"contQ") and not full_morphe:
                        qqA = atom.contQ
                    else:
                        qqA = atom.qqA
                else:
                    qqA = atom.q
                    qqB = atom.qB
                print('%6d%11s%7d%7s%7s%7d%11.6f%11.4f%11s%11.6f%11.4f' % \
                      (atom.id, atA, atom.resnr, atom.resname, atom.name, \
                       atom.cgnr, qqA, mA, atB, qqB, mB), file=fp)
                self.qA+=qqA
                self.qB+=qqB
            else:
                print('%6d%11s%7d%7s%7s%7d%11.6f%11.4f' % \
                      (atom.id, atom.atomtype, atom.resnr, atom.resname, atom.name, \
                       atom.cgnr, atom.q, atom.m), file=fp)
                self.qA+=atom.q
                self.qB+=atom.q
        # write qB of latom to qA
        if not full_morphe:
            try:
                latom.contQ = latom.qqB
            except:
                pass


##     def write_atoms(self,fp):
##         print >>fp, '[ atoms ]'
##         print >>fp,'; nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB'
##         for atom in self.atoms:
##             if atom.atomtypeB is not None:
##                 print >>fp , '%6d%11s%7d%7s%7s%7d%11.6f%11.4f%11s%11.6f%11.4f' % \
##                       (atom.id, atom.atomtype, atom.resnr, atom.resname, atom.name, \
##                        atom.cgnr, atom.q, atom.m, atom.atomtypeB, atom.qB, atom.mB)
##             else:
##                 print >>fp , '%6d%11s%7d%7s%7s%7d%11.6f%11.4f' % \
##                       (atom.id, atom.atomtype, atom.resnr, atom.resname, atom.name, \
##                        atom.cgnr, atom.q, atom.m)

    def write_bonds(self,fp, state = 'AB'):

        print('\n [ bonds ]', file=fp)
        print(';  ai    aj funct            c0            c1            c2            c3', file=fp)
        for b in self.bonds:
            if len(b) == 3:
                print('%6d %6d %6d' % (b[0].id, b[1].id, b[2]), file=fp)
            elif len(b) == 4:
                s = '   '+'   '.join([str(x) for x in b[3]])
                print('%6d %6d %6d %s' % (b[0].id, b[1].id, b[2], s), file=fp)
            else:
                lA = b[3][1]
                kA = b[3][2]
                lB = b[4][1]
                kB = b[4][2]
                if state == 'AB':
                    print('%6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % \
                          (b[0].id, b[1].id, b[2],lA,kA, lB, kB), file=fp)
                elif state == 'AA':
                    print('%6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % \
                          (b[0].id, b[1].id, b[2],lA, kA, lA, kA), file=fp)
                elif state == 'BB':
                    print('%6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % \
                          (b[0].id, b[1].id, b[2],lB, kB, lB, kB), file=fp)


    def write_pairs(self, fp):
        # CHECK HOW THIS GOES WITH B-STATES
        print('\n [ pairs ]', file=fp)
        print(';  ai    aj funct            c0            c1            c2            c3', file=fp)
        for p in self.pairs:
            print('%6d %6d %6d' % (p[0].id, p[1].id, p[2]), file=fp)

    def write_constraints(self, fp):
        # CHECK HOW THIS GOES WITH B-STATES
        print('\n [ constraints ]', file=fp)
        print(';  ai    aj funct            c0            c1            c2            c3', file=fp)
        for p in self.constraints:
            print('%6d %6d %6d' % (p[0].id, p[1].id, p[2]), file=fp)

    def write_angles(self,fp, state='AB'):
        print('\n [ angles ]', file=fp)
        print(';  ai    aj    ak funct            c0            c1            c2            c3', file=fp)
        for ang in self.angles:
            if len(ang) == 4:
                print('%6d %6d %6d %6d' % (ang[0].id, ang[1].id, ang[2].id,ang[3]), file=fp)
            else:
                if state == 'AB':
                    print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1], \
                           ang[4][2], ang[5][1], ang[5][2], ang[0].name, ang[1].name, ang[2].name), file=fp)
                elif state == 'AA':
                    print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1], \
                           ang[4][2], ang[4][1], ang[4][2], ang[0].name, ang[1].name, ang[2].name), file=fp)
                elif state == 'BB':
                    print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[5][1], \
                           ang[5][2], ang[5][1], ang[5][2], ang[0].name, ang[1].name, ang[2].name), file=fp)


    def write_dihedrals(self, fp, state='AB'):
        print('\n [ dihedrals ]', file=fp)
        print(';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5', file=fp)
        for d in self.dihedrals:
            if len(d) == 5:
                print("%6d %6d %6d %6d %4d" % ( d[0].id, d[1].id, d[2].id, d[3].id, d[4]), file=fp)
            elif len(d) == 6:
                print("%6d %6d %6d %6d %4d %s" % ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], d[5]), file=fp)
            elif len(d) == 7:
                A, B = self.check_case(d[:4])
                ast = d[5]
                bs = d[6]
                if ast == None or bs == None:
                    print(d[0].name, d[1].name, d[2].name, d[3].name, d[0].atomtype, d[1].atomtype, d[2].atomtype, d[3].atomtype, d[0].atomtypeB, d[1].atomtypeB, d[2].atomtypeB, d[3].atomtypeB)
                    print(d[0].type, d[1].type, d[2].type, d[3].type, d[0].typeB, d[1].typeB, d[2].typeB, d[3].typeB)
                if ast == 'NULL':
                    if d[4] == 3: # Ryckaert-Bellemans
                        ast = ' '.join(["%g" % x for x in [0,0,0,0,0,0]])
                    elif d[4] == 1:
                        ast = ' '.join(["%g" % x for x in [0,0,0]])
                elif ast != 'NULL' and hasattr(ast,"append"):
                    ast = ' '.join(["%g" % x for x in d[5][1:]])
                if bs == 'NULL':
                    if d[4] == 3:
                        bs = ' '.join(["%g" % x for x in [0,0,0,0,0,0]])
                    elif d[4] == 1:
                        bs = ' '.join(["%g" % x for x in [0,0,0]])
                elif bs !='NULL' and hasattr(bs,"append"):
                    bs = ' '.join(["%g" % x for x in d[6][1:]])
                if state == 'AB':
                    print("%6d %6d %6d %6d %4d %s %s ; %s %s %s %s %s %s %s %s (%s->%s)" % \
                          ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], ast, bs, d[0].name,d[1].name,d[2].name,d[3].name, \
                            d[0].type,d[1].type,d[2].type,d[3].type,A,B), file=fp)
                elif state == 'AA':
                    print("%6d %6d %6d %6d %4d %s %s ; %s %s %s %s %s %s %s %s (%s->%s)" % \
                          ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], ast, ast, d[0].name,d[1].name,d[2].name,d[3].name, \
                            d[0].type,d[1].type,d[2].type,d[3].type, A,B), file=fp)
                elif state == 'BB':
                    print("%6d %6d %6d %6d %4d %s %s ; %s %s %s %s %s %s %s %s (%s->%s)" % \
                          ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], bs, bs, d[0].name,d[1].name,d[2].name,d[3].name, \
                            d[0].type,d[1].type,d[2].type,d[3].type, A,B), file=fp)


    def check_case(self, atoms):
        A = ''
        B = ''
        for a in atoms:
            if a.atomtype.startswith('DUM'): A += 'D'
            else: A += 'A'
            if a.atomtypeB is not None:
                if a.atomtypeB.startswith('DUM'): B += 'D'
                else: B += 'A'
            else: B += 'A'
        return A, B




    def write_vsites3(self, fp):
        print('\n [ virtual_sites3 ]', file=fp)
        print(';  ai    aj    ak    al funct            c0            c1', file=fp)
        for vs in self.virtual_sites3:
            if len(vs) == 6:
                print("%6d %6d %6d %6d %4d" % ( vs[0].id, vs[1].id, vs[2].id, vs[3].id, vs[4]), file=fp)
            else:
                sys.stderr.write('EEK! Something went wrong while writing virtual_sites3!!!!\n')
                print(vs)
                sys.exit(1)

    def write_vsites4(self, fp):
        print('\n [ virtual_sites4 ]', file=fp)
        print(';  ai    aj    ak    al    am  funct            c0            c1          c2', file=fp)
        for vs in self.virtual_sites4:
            if len(vs) == 7:
                print("%6d %6d %6d %6d %6d %4d" % ( vs[0].id, vs[1].id, vs[2].id, vs[3].id, vs[4].id, vs[5]), file=fp)
            else:
                sys.stderr.write('EEK! Something went wrong while writing virtual_sites4!!!!\n')
                print(vs)
                sys.exit(1)


    def write_system(self,fp):
        print('[ system ]', file=fp)
        print(self.system, file=fp)

    def write_molecules(self,fp):
        print('[ molecules ]', file=fp)
        for mol, num in self.molecules:
            print("%s %d" % (mol,num), file=fp)

    def assign_fftypes(self):
        for atom in self.atoms:
            atom.type = self.NBParams.atomtypes[atom.atomtype]['bond_type']
            if atom.atomtypeB is not None:
                atom.typeB = self.NBParams.atomtypes[atom.atomtypeB]['bond_type']
            else:
                atom.typeB = atom.type

    def make_bond_params(self):
        for i, (at1,at2,func) in enumerate(self.bonds):
            param = get_bond_param(at1.type,at2.type,self.bond_lib)
            if param is None:
                print('Error! No bonded parameters found! (%s-%s)' % \
                      (at1.type, at2.type))
                sys.exit(1)
            self.bonds[i].append(param[1:])

    def make_angle_params(self):
        for i, (at1, at2, at3, func) in enumerate(self.angles):
            param = get_angle_param(at1.type, at2.type, at3.type, self.ang_lib)
            if param is None:
                print('Error! No angle parameters found! (%s-%s-%s)' % \
                      (at1.type, at2.type, at3.type))
                sys.exit(1)
            self.angles[i].append(param[1:])

    def make_dihedral_params(self):
        for i, d in enumerate(self.dihedrals):
            if d[5]!='': # we have a prefefined dihedral
                continue
            else:
                at1, at2, at3, at4, func, dih = d
                param = get_dihedral_param(at1.type, at2.type, \
                                           at3.type, at4.type, \
                                           self.dih_lib, func)
                if param is None:
                    print('Error! No dihedral parameters found! (%s-%s-%s-%s)' % \
                          (at1.type, at2.type, at3.type, at4.type))
                    print(func, dih)
                    sys.exit(1)
                del self.dihedrals[i][-1]
                self.dihedrals[i].append(param[1:])




def cpp_parse_file(fn,cpp_defs=[],cpp_path=[os.environ.get('GMXLIB')] ):

    defs = []
    incs = []
    for d in cpp_defs:
        defs.append('-D%s' % d)
    for i in cpp_path:
        incs.append('-I%s' % i)
    cmd = 'cpp -traditional %s %s %s ' % (' '.join(defs),' '.join(incs),fn)
    return os.popen(cmd,'r').readlines()

def read_atp(fn):
    l = open(fn).readlines()
    l = kickOutComments(l,';')
    l = parseList('sf',l)
    dic = {}
    for type,mass in l:
        dic[type] = mass
    return dic

def read_atomtypes(l, ff= 'amber99sb'):
    lst = readSection(l,'[ atomtypes ]','[')
    if ff == 'oplsaa':
        try:
            lst = parseList('ssiffsff',lst)
        except:
            lst = parseList('sffsff',lst)
    elif ff in ['amber03','amber99sb']:
        lst = parseList('ssffsff',lst)
    elif ('gaff' in ff) or ('cgenff' in ff):
        lst = parseList('sffsff',lst)
    dic = {}
    for line in lst:
        dic[line[0]]=line[1:]
    return dic

def read_bondtypes(l):
    res = []
    starts = []
    for i, line in enumerate(l):
        if line.strip().startswith('[ bondtypes ]'):
            starts.append(i)
    for s in starts:
        lst = readSection(l[s:],'[ bondtypes ]','[')
        lst = parseList('ssiff',lst)
        res.extend(lst)
    return res

def read_angletypes(l):
    res = []
    starts = []
    for i, line in enumerate(l):
        if line.strip().startswith('[ angletypes ]'):
            starts.append(i)
    for s in starts:
        lst = readSection(l[s:],'[ angletypes ]','[')
        lst = parseList('sssiff',lst)
        res.extend(lst)
    return res

def read_dihedraltypes(l):
    res = []
    starts = []
    for i, line in enumerate(l):
        if line.strip().startswith('[ dihedraltypes ]'):
            starts.append(i)
    for s in starts:
        lst = readSection(l[s:],'[ dihedraltypes ]','[')
        try:
            lst = parseList('ssssiffffff',lst)
        except:
            try:
                lst = parseList('ssssiffi',lst)
            except:
                lst = parseList('ssiffi',lst)
        res.extend(lst)
    return res


def read_ff(fn,ff='amber99sb',cpp_defs = [], cpp_path = [os.environ.get('GMXLIB')]):
    l = cpp_parse_file(fn,cpp_defs,cpp_path)
    l = kickOutComments(l,'#')
    l = kickOutComments(l,';')
    atomtypes = read_atomtypes(l,ff=ff)
    bt = read_bondtypes(l)
    at = read_angletypes(l)
    dt = read_dihedraltypes(l)
    return (atomtypes,bt,at,dt)


def __get_rtp_resnames( lines ):
    keys = []
    for line in lines:
        if line.strip().startswith('['):
            if line.strip()[1:-1].strip() not in \
                   ['atoms','bonds','dihedrals','impropers','bondedtypes']:
                keys.append( line.strip()[1:-1].strip() )
    return keys

def __get_rtp_entry( key, lines ):
    r = []
    for line in lines:
        if line.strip()[1:-1].strip() == key:
            idx = lines.index( line )
    for line in lines[idx+1:]:
        if line.strip().startswith('['):
            if line.strip()[1:-1].strip() not in \
               ['atoms','bonds','dihedrals','impropers']:
                break
            else:
                r.append(line)
        else:
            r.append(line)
    return r

def __read_rtp_atoms(resname, lines ):
    atoms = []
    for line in lines:
        entr = line.split()
        if resname in _aliases :
            if entr[0] in _aliases[resname]:
                entr[0] = _aliases[resname][entr[0]]
        entr[2] = float(entr[2])
        entr[3] = int(entr[3])
        atoms.append(entr)
    return atoms

def __read_rtp_bonds( resname, lines ):
    bonds = []
    for line in lines:
        entr = line.split()
        if entr[0] not in ['+','-']:
            if resname in _aliases :
                if entr[0] in _aliases[resname]:
                    entr[0] = _aliases[resname][entr[0]]
        if entr[1] not in ['+','-']:
            if resname in _aliases :
                if entr[1] in _aliases[resname]:
                    entr[1] = _aliases[resname][entr[1]]
        bonds.append(entr)
    return bonds

def __read_rtp_dihedrals( resname, lines ):
    diheds = []
    for line in lines:
        entr = line.split()
        if entr[0] not in ['+','-']:
            if resname in _aliases :
                if entr[0] in _aliases[resname]:
                    entr[0] = _aliases[resname][entr[0]]
        if entr[1] not in ['+','-']:
            if resname in _aliases :
                if entr[1] in _aliases[resname]:
                    entr[1] = _aliases[resname][entr[1]]
        if entr[2] not in ['+','-']:
            if resname in _aliases :
                if entr[2] in _aliases[resname]:
                    entr[2] = _aliases[resname][entr[2]]
        if entr[3] not in ['+','-']:
            if resname in _aliases :
                if entr[3] in _aliases[resname]:
                    entr[3] = _aliases[resname][entr[3]]
        if len(entr) == 5:
            diheds.append(entr)
        else:
            diheds.append(entr+[''])
    return diheds

def __read_rtp_impropers( resname, lines ):
    improps = []
    for line in lines:
        entr = line.split()
        if entr[0] not in ['+','-']:
            if resname in _aliases :
                if entr[0] in _aliases[resname]:
                    entr[0] = _aliases[resname][entr[0]]
        if entr[1] not in ['+','-']:
            if resname in _aliases :
                if entr[1] in _aliases[resname]:
                    entr[1] = _aliases[resname][entr[1]]
        if entr[2] not in ['+','-']:
            if resname in _aliases :
                if entr[2] in _aliases[resname]:
                    entr[2] = _aliases[resname][entr[2]]
        if entr[3] not in ['+','-']:
            if resname in _aliases :
                if entr[3] in _aliases[resname]:
                    entr[3] = _aliases[resname][entr[3]]
        if len(entr) == 5:
            improps.append(entr)
        else:
            improps.append(entr+[''])
    return improps


def read_rtp( filename ):
    rtp_entries = {}
    l = open(filename,'r').readlines()
    lines = kickOutComments(l,';')
    keys = __get_rtp_resnames(lines)
    for key in keys:
        sys.stderr.write('%s -> %4s\r' % (filename, key))
        rtp_lines = __get_rtp_entry( key, lines)
        # read atoms
        al = readSection(rtp_lines,'[ atoms ]','[')
        atoms = __read_rtp_atoms(key, al )
        # read bonds
        bl = readSection(rtp_lines,'[ bonds ]','[')
        bonds = __read_rtp_bonds( key, bl )
        # read dihedrals
        dl = readSection(rtp_lines,'[ dihedrals ]','[')
        diheds = __read_rtp_dihedrals(key, dl)
        # read impropers
        il = readSection(rtp_lines,'[ impropers ]','[')
        improps = __read_rtp_impropers(key, il)
        rtp_entries[key] = {
            'atoms': atoms,
            'bonds': bonds,
            'diheds': diheds,
            'improps': improps
            }
    sys.stderr.write('\ndone...\n' )
    return rtp_entries



def get_rtp_entry(key, filename = 'ffamber99sb.rtp'):
    l = open(filename,'r').readlines()
    lines = kickOutComments(l,';')

    r = []
    idx = 0
    for line in lines:
        if line.strip()[1:-1].strip() == key:
            idx = lines.index(line)
            break
    for line in lines[idx+1:]:
        if line.strip().startswith('['):
            if line.strip()[1:-1].strip() not in \
               ['atoms','bonds','dihedrals','impropers']:
                break
            else:
                r.append(line)
        else:
            r.append(line)
    # read atoms
    atoms = []
    al = readSection(r,'[ atoms ]','[')
    for line in al:
        entr = line.split()
        entr[2] = float(entr[2])
        entr[3] = int(entr[3])
        atoms.append(entr)
    # read bonds
    bonds = []
    bl = readSection(r,'[ bonds ]','[')
    for line in bl:
        entr = line.split()
        bonds.append(entr)
    # read dihedrals
    diheds = []
    dl = readSection(r,'[ dihedrals ]','[')
    for line in dl:
        entr = line.split()
        if len(entr) == 5:
            diheds.append(entr)
        else:
            diheds.append(entr+[''])
    # read impropers
    improps = []
    il = readSection(r,'[ impropers ]','[')
    for line in il:
        entr = line.split()
        if len(entr)==5:
            improps.append(entr)
        else:
            improps.append(entr+[''])
    return atoms, bonds, diheds, improps

def topline2atom(line):
    entr = line.split()
    idx = int(entr[0])
    atomtype = entr[1]
    resnr = int(entr[2])
    resname = entr[3]
    name = entr[4]
    cgnr = int(entr[5])
    q = float(entr[6])
    m = float(entr[7])
    try:
        atomtypeB = entr[8]
        qB = float(entr[9])
        mB = float(entr[10])
    except:
        atomtypeB = None
        qB = None
        mB = None
    a = Atom(id=idx,atomtype=atomtype,\
        resnr = resnr, resname = resname,\
        name = name, cgnr = cgnr, q = q, \
        m = m, atomtypeB = atomtypeB, \
        qB = qB, mB = mB)
    return a

def read_itp_atoms(lines):
    lst = readSection(lines,'[ atoms ]','[')
    al = []
    for line in lst:
        a = topline2atom(line)
        al.append(a)
    return al

def write_itp_atoms(fp,al):
    print('[ atoms ]', file=fp)
    print('; nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB', file=fp)
    for atom in al:
        if atom.atomtypeB is not None:
            print('%6d %11s %7d%7s%7s%7d%11.6f%11.4f %11s%11.6f%11.4f' % \
                 (atom.id, atom.atomtype, atom.resnr, atom.resname, atom.name, \
                  atom.cgnr, atom.q, atom.m, atom.atomtypeB, atom.qB, atom.mB), file=fp)
        else:
            print('%6d %11s %7d%7s%7s%7d%11.6f%11.4f' % \
                  (atom.id, atom.atomtype, atom.resnr, atom.resname, atom.name, \
                   atom.cgnr, atom.q, atom.m), file=fp)

def write_itp_bonds(fp,bonds):
    print('[ bonds ]', file=fp)
    for b in bonds:
        if isinstance(b[0],Atom):
            tmp = [b[0].id, b[1].id]+b[2:]
        else:
            tmp = b
        s = [str(item) for item in tmp]
        out = ''
        for item in s:
            if len(item) > 6:
                form="%%%ds " % len(item)
                out += form % item
            elif len(item) > 0:
                out += "%6s " % item
        print(out, file=fp)

def write_itp_pairs(fp,pairs):
    print('[ pairs ]', file=fp)
    for p in pairs:
        if isinstance(p[0],Atom):
            tmp = [p[0].id, p[1].id]+p[2:]
        else:
            tmp = p
        s = [str(item) for item in tmp]
        out = ''
        for item in s:
            if len(item) > 6:
                form="%%%ds " % len(item)
                out += form % item
            elif len(item) > 0:
                out += "%6s " % item
        print(out, file=fp)

def write_itp_angles(fp,angles):
    print('[ angles ]', file=fp)
    for a in angles:
        if isinstance(a[0],Atom):
            tmp = [a[0].id, a[1].id, a[2].id]+a[3:]
        else:
            tmp = a
        s = [str(item) for item in tmp]
        out = ''
        for item in s:
            if len(item) > 6:
                form="%%%ds " % len(item)
                out += form % item
            elif len(item) > 0:
                out += "%6s " % item
        print(out, file=fp)


def write_itp_dihedrals(fp,dihedrals):
    print('[ dihedrals ]', file=fp)
    for d in dihedrals:
        if isinstance(d[0],Atom):
            tmp = [d[0].id, d[1].id, d[2].id, d[3].id]+d[4:]
        else:
            tmp = d
        s = [str(item) for item in tmp]
        out = ''
        for item in s:
            if len(item) > 6:
                form="%%%ds " % len(item)
                out += form % item
            elif len(item) > 0:
                out += "%6s " % item
        print(out, file=fp)


def write_itp_moleculetype(fp,name,nrexcl):
    print('[ moleculetype ]', file=fp)
    print('; Name        nrexcl', file=fp)
    print('%s  %d' % (name,nrexcl), file=fp)

def read_itp_bonds(lines):
    lst = readSection(lines,'[ bonds ]','[')
    bonds = []
    for line in lst:
        entr = line.split()
        b0 = int(entr[0])
        b1 = int(entr[1])
        bt = int(entr[2])
        if len(entr) > 3:
            params = [float(x) for x in entr[3:]]
        else:
            params = []
        bonds.append([b0,b1,bt]+params)
    return bonds

def read_itp_pairs(lines):
    lst = readSection(lines,'[ pairs ]','[')
    pairs = []
    for line in lst:
        entr = line.split()
        b0 = int(entr[0])
        b1 = int(entr[1])
        bt = int(entr[2])
        pairs.append([b0,b1,bt])
    return pairs

def read_itp_angles(lines):
    lst = readSection(lines,'[ angles ]','[')
    angles = []
    for line in lst:
        entr = line.split()
        b0 = int(entr[0])
        b1 = int(entr[1])
        b2 = int(entr[2])
        bt = int(entr[3])
        if len(entr) > 4:
            params = [float(x) for x in entr[4:]]
        else:
            params = []
        angles.append([b0,b1,b2,bt]+params)
    return angles

def read_itp_dihedrals(lines):
    starts = []
    dihedrals = []
    for i, line in enumerate(lines):
        if line.strip().startswith('[ dihedrals ]'):
            starts.append(i)
    for s in starts:
        lst = readSection(lines[s:],'[ dihedrals ]','[')
        for line in lst:
            entr = line.split()
            b0 = int(entr[0])
            b1 = int(entr[1])
            b2 = int(entr[2])
            b3 = int(entr[3])
            bt = int(entr[4])
            if len(entr) > 5:
                try:
                    params = [float(x) for x in entr[5:]]
                except:
                    params = entr[5:] # if we have defines
            else:
                params = []
            dihedrals.append([b0,b1,b2,b3,bt]+params)
    return dihedrals



def read_moleculetype(lines):
    l = readSection(lines,'[ moleculetype ]','[')
    if l:
        return l[0].split()[0], int(l[0].split()[1])
    else:
        return None, None

def read_gaff_top(fname):
    """ this function reads topology files from gaff """

    lines = open(fname).readlines()
    lines = kickOutComments(lines,';')
    itp = ITPFile(fname, ff = 'amber03')
#    atypes = read_atomtypes(lines, ff = 'amber03')
#    itp.atomtypes = atypes
    return itp

class MDPError(Exception):
    def __init__(self, s):
        self.s = s
    def __str__(self):
        return repr(self.s)


class MDP:

    def __init__(self, fn = None):
        self.parameters = OrderedDict([
        ['include'                  ,''],
        ['define'                   ,''],
        ['integrator'               , 'md'],
        ['tinit'                    , 0],
        ['dt'                       , 0.002],
        ['nsteps'                   , 25000],
        ['simulation_part'          , 1],
        ['init_step'                , 0],
        ['comm-mode'                , 'Linear'],
        ['nstcomm'                  , 100],
        ['nstcalcenergy'            , 100],
        ['nstdhdl'                  , 100],
        ['comm-grps'                ,''],
        ['bd-fric'                  , 0],
        ['ld-seed'                  , 1993],
        ['emtol'                   , 100],
        ['emstep'                   , 0.01],
        ['niter'                    , 0],
        ['fcstep'                   , 0],
        ['nstcgsteep'               , 1000],
        ['nbfgscorr'                , 10],
        ['rtpi'                     , 0.05],
        ['nstxout'                  , 10000],
        ['nstvout'                  , 10000],
        ['nstfout'                  , 0],
        ['nstlog'                   , 1000],
        ['nstenergy'                , 100],
        ['nstxtcout'                , 100],
        ['xtc-precision'            , 1000],
        ['xtc-grps'                 ,''],
        ['energygrps'               , ''],
        ['nstlist'                  , 10],
        ['ns-type'                  , 'Grid'],
        ['pbc'                      , 'xyz'],
        ['periodic_molecules'       , 'no'],
        ['rlist'                    , 1.2],
        ['cutoff-scheme'            , 'verlet'],
        ['coulombtype'              , 'PME'],
        ['rcoulomb-switch'          , 0],
        ['rcoulomb'                 , 1.1],
        ['epsilon-r'                , 1],
        ['epsilon_rf'               , 1],
        ['vdw-modifier'             , 'Potential-switch'],
        ['rvdw-switch'              , 1],
        ['rvdw'                     , 1.1],
        ['DispCorr'                 , 'EnerPres'],
        ['table-extension'          , 1],
        ['energygrp_table'          ,''],
        ['fourierspacing'           , 0.12],
        ['fourier_nx'               , 0],
        ['fourier_ny'               , 0],
        ['fourier_nz'               , 0],
        ['pme_order'                , 4],
        ['ewald_rtol'               , 1e-05],
        ['ewald_geometry'           , '3d'],
        ['epsilon_surface'          , 0],
        ['optimize_fft'             , 'no'],
        ['implicit_solvent'         , 'No'],
#        ['gb_algorithm'             , 'Still'],
#        ['nstgbradii'               , 1],
#        ['rgbradii'                 , 2],
#        ['gb_epsilon_solvent'       , 80],
#        ['gb_saltconc'              , 0],
#        ['gb_obc_alpha'             , 1],
#        ['gb_obc_beta'              , 0.8],
#        ['gb_obc_gamma'             , 4.85],
#        ['sa_surface_tension'       , 2.092],
        ['tcoupl'                   , 'v-rescale'],
        ['tc-grps'                  , ['System']],
        ['tau-t'                    , [0.1]],
        ['ref-t'                    , [298]],
        ['Pcoupl'                   , 'Parrinello-Rahman'],
        ['Pcoupltype'               , 'Isotropic'],
        ['tau-p'                    , 5],
        ['compressibility'          , 4.6E-5],
        ['ref-p'                    , 1],
        ['refcoord_scaling'         , 'all'],
        ['andersen_seed'            , 815131],
        ['QMMM'                     , 'no'],
#        ['QMMM-grps'                ,''],
#        ['QMmethod'                 ,''],
#        ['QMMMscheme'               , 'normal'],
#        ['QMbasis'                  ,''],
#        ['QMcharge'                 ,''],
#        ['QMmult'                   ,''],
        ['annealing'                , ['no']],
        ['annealing_npoints'        , [2]],
        ['annealing_time'           , [0, 50]],
        ['annealing_temp'           , [0, 298]],
        ['gen-vel'                  , 'no'],
        ['gen-temp'                 , 298],
        ['gen-seed'                 , 173529],
        ['constraints'              , 'all-bonds'],
        ['constraint-algorithm'     , 'Lincs'],
#        ['continuation'             , 'yes'],
        ['Shake-SOR'                , 'no'],
        ['shake-tol'                , 1e-04],
        ['lincs-order'              , 4],
        ['lincs-iter'               , 1],
        ['lincs-warnangle'          , 30],
        ['morse'                    , 'no'],
        ['energygrp_excl'           ,''],
        ['nwall'                    , 0],
        ['wall_type'                , '9-3'],
        ['wall_r_linpot'            , -1],
        ['wall_atomtype'            ,''],
        ['wall_density'             ,''],
        ['wall_ewald_zfac'          , 3],
        ['pull'                    , 'no'],
        ['disre'                    , 'No'],
        ['disre-weighting'          , 'Equal'],
        ['disre-mixed'              , 'no'],
        ['disre-fc'                 , 1000],
        ['disre-tau'                , 0],
        ['nstdisreout'              , 100],
        ['orire'                    , 'no'],
        ['orire-fc'                 , 0],
        ['orire-tau'                , 0],
        ['orire-fitgrp'             ,''],
        ['nstorireout'              , 100],
        ['dihre'                    , 'No'],
        ['dihre-fc'                 , 1000],
        ['free-energy'              , 'no'],
        ['init-lambda'              , 0],
        ['delta-lambda'             , 0],
        ['sc-alpha'                 , 0.3],
        ['sc-power'                 , 1],
        ['sc-sigma'                 , 0.25],
        ['sc-coul'                  , 'yes'],
        ['couple-moltype'           ,''],
        ['couple-lambda0'           , 'vdw-q'],
        ['couple-lambda1'           , 'vdw-q'],
        ['couple-intramol'          , 'no'],
        ['acc-grps'                 ,''],
        ['accelerate'               ,''],
        ['freezegrps'               ,''],
        ['freezedim'                ,''],
        ['cos-acceleration'         , 0],
        ['deform'                   ,''],
        ['E-x'                      ,''],
        ['E-xt'                     ,''],
        ['E-y'                      ,''],
        ['E-yt'                     ,''],
        ['E-z'                      ,''],
        ['E-zt'                     ,''],
        ['user1-grps'               ,''],
        ['user2-grps'               ,''],
        ['userint1'                 , 0],
        ['userint2'                 , 0],
        ['userint3'                 , 0],
        ['userint4'                 , 0],
        ['userreal1'                , 0],
        ['userreal2'                , 0],
        ['userreal3'                , 0],
        ['userreal4'                , 0]
        ])

        if fn:
            self.read(fn)

    def __str__(self):
        line = ''
        for key, val in list(self.parameters.items()):
            if hasattr(val,"append"):
                s = ''
                for x in val:
                    s+=str(x)+' '
            else:
                s = str(val)
            line+="%-25s = %s\n" % (key, s)
        return line

    def __setitem__(self,item,value):
        if item not in self.parameters:
            raise MDPError("No such option %s" % item)

        self.parameters[item] = value

    def __getitem__(self, item):
        return self.parameters[item]


    def write(self, fp = None):

        if fp is None:
            fp = sys.stdout
        else:
            if not hasattr(fp,"write"):
                fp = open(fp,"w")
        print(self, file=fp)

    def read(self, filename):
        lines = open(filename).readlines()
        l = kickOutComments(lines,';')
        for line in l:
            entr = line.split('=')
            key = entr[0].strip()
            val = entr[1].strip().split()
            if key not in self.parameters:
                self.parameters[key] = val
#                print 'Warning! Ignoring entry \'%s\'' % key
            else:
                if len(val) == 0:
                    self[key] = ''
                elif len(val) == 1:
                    self[key] = val[0]
                else:
                    self[key] = val
        return self



def make_amber_residue_names(model):

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
                    if d < 2.5:
                        ss_bond = True
                        break
            if ss_bond:
                # terminal cys2 is ccyx
                rr = 'CYS2'
                res.set_resname(rr)
            else:
                res.set_resname('CYM')

        else:
            res.set_resname('CYN')
    lysl = model.fetch_residues('LYS')
    for res in lysl:
        at = res.fetch('HZ3')
        at2 = res.fetch('HZ2')
        if at or not at2:
            res.set_resname('LYP')
    # histidine
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

    aspl = model.fetch_residues('ASP')
    for res in aspl:
        bHD2 = False
        hd2 = res.fetch('HD2')
        if hd2:
            res.set_resname('ASH')

    glul = model.fetch_residues('GLU')
    for res in glul:
        bHD2 = False
        hd2 = res.fetch('HE2')
        if hd2:
            res.set_resname('GLH')
    for chain in model.chains:
        if chain.residues[0].is_protein_residue():
            first = chain.nterminus()
            last = chain.cterminus()
            first.set_resname('N'+first.resname) # rename e.g. ALA to NALA
            if last.resname == 'CYS2':
                last.set_resname('CCYX')   # rename e.g. ARG to CARG
            else:
                last.set_resname('C'+last.resname)   # rename e.g. ARG to CARG
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
                    print('pymacs_Warning_> No terminal oxygen atoms found in chain %s' % chain.id, file=sys.stderr)



def assign_ffamber99sb_params(m):
    m.get_symbol()
    m.rename_atoms()
    for c in m.chains:
        c.make_residue_tree()

    make_amber_residue_names( m)
    rtp = RTPParser('ffamber99sb.rtp')
    rtp.assign_params(m)
    bo = BondedParser('ffamber99sbbon.itp')
    nb = NBParser('ffamber99sbnb.itp')
    nb.assign_params( m )
    bo.assign_params( m )
    rtp.assign_dihedral_params( m, bo.directives )



def bond_energy(m):
    return _p.total_bond_energy(m.bond_list)
def angle_energy(m):
    return _p.total_angle_energy(m.angle_list)
def dihedral_energy(m):
    return _p.total_dihedral_energy(m.dihedral_list)
def improper_energy(m):
    return _p.total_improper_energy(m.improper_list)
def coul14_energy(m):
    return _p.coul14_energy( m.atoms )
def lj14_energy( m ):
    return _p.lj14_energy( m.atoms )
def nb_lj_energy( m ):
    return _p.nb_lj_energy( m.atoms )
def nb_coul_energy( m ):
    return _p.nb_coul_energy( m.atoms )
def nb_energy( m ):
    return _p.nb_energy( m.atoms )


def energy(m):

    bond_ene = bond_energy( m )
    angle_ene = angle_energy( m )
    dihedral_ene = dihedral_energy( m )
    improper_ene = improper_energy ( m )
    lj14_ene = lj14_energy( m )
    coul14_ene = coul14_energy( m )
    nb_ene = nb_energy( m )
##     print 'bonds = ', bond_ene
##     print 'angles = ',angle_ene
##     print 'dihedrals = ',dihedral_ene
##     print 'impropers = ',improper_ene
##     print 'nb = ',nb_ene
##     print 'lj14 = ',lj14_ene
##     print 'coul14 = ',coul14_ene

    return bond_ene + angle_ene + dihedral_ene + improper_ene + nb_ene + lj14_ene + coul14_ene
