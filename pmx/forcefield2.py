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
from .molecule import Molecule
from .odict import *
from .library import _aliases
from .ffparser import *
from  . import _pmx as _p

def TR ( s ):
    print("pmx.forcefield_> " + s)

#def cpp_parse_file(fn,cpp_defs=[],cpp_path=[os.environ.get('GMXDATA')+'/top'] ):
def cpp_parse_file(fn,cpp_defs=[],cpp_path=[os.environ.get('GMXLIB')], itp=False, ffpath=None ):

    defs = []
    incs = []
    for d in cpp_defs:
        defs.append('-D%s' % d)
    for i in cpp_path:
        incs.append('-I%s' % i)
    if itp:
        cmd1 = 'cpp -traditional %s %s %s ' % (' '.join(defs),' '.join(incs),fn)
        l1 = os.popen(cmd1,'r').readlines()
        if ffpath != None:
            ffname = ffpath+'/forcefield.itp'
            cmd2 = 'cpp -traditional %s %s %s ' % (' '.join(defs),' '.join(incs),ffname)
            l2 = os.popen(cmd2,'r').readlines()
            return(l1+l2)
        else:
            return(l1)
    else:
        cmd = 'cpp -traditional %s %s %s ' % (' '.join(defs),' '.join(incs),fn)
        return os.popen(cmd,'r').readlines()




class TopolBase:

    def __init__(self, filename, version = 'old'):
        self.filename = filename
        self.version = version
        self.is_itp = True
        if filename!=None:
            if os.path.splitext( filename )[1] == '.itp':
                self.is_itp = True
            else:
                self.is_itp = False
        self.has_atomtypes = False
        self.atomtypes = []
        self.atoms = []
        self.residues = []
        self.name = ''
        self.nrexcl = 0
        self.bonds = []
        self.constraints = []
        self.have_constraints = False
        self.pairs = []
        self.cmap  = []
        self.angles = []
        self.dihedrals = []
        self.virtual_sites2 = []
        self.virtual_sites3 = []
        self.virtual_sites4 = []
        self.has_vsites2 = False
        self.has_vsites3 = False
        self.has_vsites4 = False
        self.has_posre = False
        self.has_exclusions = False
        self.exclusions = []
        self.has_ii = False # intermolecular_interactions
        self.ii = {} # ii is a dictionary, which containts bonds, angles, dihedrals
        self.posre = []
        self.molecules = []
        self.system = ''
        self.qA = 0.
        self.qB = 0.
        if filename!=None:
            self.read()
    #===============================================================================
    # read functions

    def read( self ):
        lines = open(self.filename).readlines()
        lines = kickOutComments(lines,';')
        if not self.is_itp:
            self.read_header( lines )
        self.read_footer( lines )
        posre_sections = self.get_posre_section( lines )
        lines = kickOutComments(lines,'#')
        self.read_atomtypes(lines)
        self.read_moleculetype(lines)
        if self.name: # atoms, bonds, ... section
            self.read_atoms(lines)
            self.read_bonds(lines)
            self.read_constraints(lines)
            self.read_pairs(lines)
            self.read_angles(lines)
            self.read_dihedrals(lines)
            self.read_cmap(lines)
            self.read_vsites2(lines)
            self.read_vsites3(lines)
            self.read_vsites4(lines)
            self.read_exclusions(lines)
            if self.has_posre:
                self.read_posre(posre_sections)
            self.__make_residues()
        if not self.is_itp:
            self.read_system(lines)
            self.read_molecules(lines)

    def __atom_from_top_line(self, line):
        entr = line.split()
        idx = int(entr[0])
        atomtype = entr[1]
        resnr = int(entr[2])
        resname = entr[3]
        name = entr[4]
        cgnr = int(entr[5])
        q = float(entr[6])
        try:
            m = float(entr[7])
        except: # martini
            m = ""
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

    def __make_residues(self):
        cur_mol = None
        mol = None
        for atom in self.atoms:
            if atom.resnr == cur_mol:
                if mol:
                    mol.atoms.append( atom )
                else:
                    mol = Molecule()
                    cur_mol = atom.resnr
                    mol.resname = atom.resname
                    mol.id = cur_mol
                    mol.atoms.append( atom )
            else:
                if mol:
                    self.residues.append( mol )
                    mol = Molecule()
                    cur_mol = atom.resnr
                    mol.resname = atom.resname
                    mol.id = cur_mol
                    mol.atoms.append( atom )
                else:
                    mol = Molecule()
                    cur_mol = atom.resnr
                    mol.resname = atom.resname
                    mol.id = cur_mol
                    mol.atoms.append( atom )
        self.residues.append( mol )
        for r in self.residues:
            atom.molecule = r



    def read_system(self,lines):
        lst = readSection(lines,'[ system ]','[')
        self.system = lst[0].strip()

    def read_molecules(self,lines):
        lst = readSection(lines,'[ molecules ]','[')
        self.molecules = []
        for line in lst:
            entr = line.split()
            self.molecules.append([entr[0],int(entr[1])])

    def read_moleculetype(self, lines):
        l = readSection(lines,'[ moleculetype ]','[')
        if l:
            self.name, self.nrexcl =  l[0].split()[0], int(l[0].split()[1])

    def read_atomtypes( self, lines):
        l = readSection(lines,'[ atomtypes ]','[')
        lst = []
        try:
            lst = parseList('ssffsff',l) # gaff
        except:
            try:
                lst = parseList('sffsff',l) # another gaff or cgenff
            except:
                try:
                    lst = parseList('ssiffsff',l) # opls
                except:
                    print('Could not read atomtype format')
        if len(lst)>0:
            self.has_atomtypes = True
        self.atomtypes = lst

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
        try:
            idx = self.footer.index('[ system ]')
            self.footer = self.footer[:idx]
        except:
            pass

    def read_atoms(self,lines):
        lst = readSection(lines,'[ atoms ]','[')
        self.atoms = []
        for line in lst:
            a = self.__atom_from_top_line( line )
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
                self.bonds.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], idx[2], [idx[2],lA,kA],[idx[2],lB,kB]])

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
            elif len(entries) == 8 and entries[3]=='1':
                idx = [int(x) for x in entries[:4]]
                lA = float(entries[4])
                kA = float(entries[5])
                lB = float(entries[6])
                kB = float(entries[7])
                self.angles.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], \
                                    self.atoms[idx[2]-1], idx[3], [idx[3],lA,kA],[idx[3],lB,kB]])
            elif len(entries) == 8 and entries[3]=='5':
                idx = [int(x) for x in entries[:4]]
                lA1 = float(entries[4])
                kA1 = float(entries[5])
                lA2 = float(entries[6])
                kA2 = float(entries[7])
                self.angles.append([self.atoms[idx[0]-1], \
                        self.atoms[idx[1]-1],self.atoms[idx[2]-1],idx[3],[idx[3],lA1,kA1,lA2,kA2]])
            elif len(entries) == 12:
                idx = [int(x) for x in entries[:4]]
                lA1 = float(entries[4])
                kA1 = float(entries[5])
                lA2 = float(entries[6])
                kA2 = float(entries[7])
                lB1 = float(entries[8])
                kB1 = float(entries[9])
                lB2 = float(entries[10])
                kB2 = float(entries[11])
                self.angles.append([self.atoms[idx[0]-1], self.atoms[idx[1]-1], \
                                    self.atoms[idx[2]-1], idx[3],\
                                    [idx[3],lA1,kA1,lA2,kA2],[idx[3],lB1,kB1,lB2,kB2]])

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
#               foo = (self.atoms[idx[0]-1],self.atoms[idx[1]-1],self.atoms[idx[2]-1],self.atoms[idx[3]-1],func,rest)
#               print 'length %d' % len(foo)
    def read_cmap(self, lines):
        starts = []
        cmap = []
        for i, line in enumerate(lines):
            if line.strip().startswith('[ cmap ]'):
                starts.append(i)
        for s in starts:
            lst = readSection(lines[s:],'[ cmap ]','[')
            for line in lst:
                entr = line.split()
                idx = [int(x) for x in entr[:5]]

                func = int(entr[5])
                try:
                    rest = ' '.join(entr[6:])
                except:
                    rest = ''
                self.cmap.append([self.atoms[idx[0]-1],\
                                       self.atoms[idx[1]-1],\
                                       self.atoms[idx[2]-1],\
                                       self.atoms[idx[3]-1],\
                                       self.atoms[idx[4]-1],\
                                       func,rest])
#               foo = (self.atoms[idx[0]-1],self.atoms[idx[1]-1],self.atoms[idx[2]-1],self.atoms[idx[3]-1],func,rest)
#               print 'length %d' % len(foo)

    def read_exclusions(self, lines):
        starts = []
        for i, line in enumerate(lines):
            if line.strip().startswith('[ exclusions ]'):
                starts.append(i)
        if starts:
            self.has_exclusions = True
        for s in starts:
            lst = readSection(lines[s:], '[ exclusions ]', '[')
            for line in lst:
                entr = line.split()
                idx = [int(x) for x in entr[:]]
                toappend = [self.atoms[x-1] for x in idx]
                self.exclusions.append( toappend )

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

    def get_posre_section( self, lines ):
        starts = []
        bIfDef = False
        for i, line in enumerate(lines):
            if line.strip().startswith('#endif'):
                bIfDef = False
            if line.strip().startswith('#ifdef'):
                bIfDef = True
            if bIfDef:
                continue
            if line.strip().startswith('[ position_restraints ]'):
                starts.append(i)
        if starts:
            self.has_posre = True

        lstList = {}
        counter = 0
        for s in starts:
            lst = readSection(lines[s:],'[ position_restraints ]','[')
            lst = kickOutComments(lst,'#')
            lstList[counter] = lst
            counter+=1

        return( lstList )

    def read_posre(self, lstList):
        for lstKey in lstList:
            lst = lstList[lstKey]
            for line in lst:
                entr = line.split()
                idx = int(entr[0])

                func = int(entr[1])
                try:
                    rest = ' '.join(entr[2:])
                except:
                    rest = ''
                self.posre.append([self.atoms[idx-1],func,rest])


    #===============================================================================
    # write functions



    def write(self, outfile, stateBonded = 'AB', stateTypes = 'AB', stateQ = 'AB',
              scale_mass = False, dummy_qA = 'on', dummy_qB = 'on', target_qB = [],
              full_morphe = True):

        fp = open(outfile,'w')
        # determine if outfile is .itp
        is_out_itp = False
        if os.path.splitext( outfile )[1] == '.itp':
            is_out_itp = True

        if not is_out_itp:
            self.write_header(fp)
        if self.atoms:
            if self.has_atomtypes:
                self.write_atomtypes(fp)
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
            self.write_cmap(fp)
            if self.has_vsites2:
                self.write_vsites2(fp)
            if self.has_vsites3:
                self.write_vsites3(fp)
            if self.has_vsites4:
                self.write_vsites4(fp)
            if self.has_exclusions:
                self.write_exclusions(fp)
            if self.has_posre:
                self.write_posre(fp)
        if not is_out_itp:
            self.write_footer(fp)
            self.write_system(fp)
            self.write_molecules(fp)
        if self.atoms:
            if self.has_ii:
                self.write_ii(fp)
        fp.close()


    def write_header(self,fp):
        for line in self.header:
            print(line, file=fp)

    def write_footer(self,fp):
        try:
            for line in self.footer:
                print(line, file=fp)
        except:
            print("No footer in itp\n")

    def write_moleculetype(self, fp):
        print('[ moleculetype ]', file=fp)
        print('; Name        nrexcl', file=fp)
        print('%s  %d' % (self.name,self.nrexcl), file=fp)

    def write_atomtypes(self, fp):
        fp.write('[ atomtypes ]\n')
        for atype in self.atomtypes:
            for x in atype:
                fp.write(' %s' % x)
            fp.write('\n')
        fp.write('\n')

    def write_atoms(self, fp, charges = 'AB', atomtypes = 'AB', dummy_qA = 'on',\
                dummy_qB = 'on',  scale_mass=True, target_qB = [], full_morphe = True):

        self.qA = 0
        self.qB = 0
        for r in self.residues:
            if self.__is_perturbed_residue(r):
                try:
                    target_chargeB = target_qB.pop(0)
                except:
                    target_chargeB = 0
                TR( 'Making target charge %g for residue %s' % (round(target_chargeB,5), r.resname) )
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
                if round(qB_tot,5) != round(target_chargeB,5):
                    TR('State B has total charge of %g' % round(qB_tot,5))
                    TR('Applying charge correction to ensure integer charges')
                    latom = self.__last_perturbed_atom(r)
                    TR('Selecting atom %d-%s (%s) as perturbed atom with highest order' % (latom.id,latom.name, latom.resname))
                    newqB = latom.qqB-(qB_tot-target_chargeB)
                    TR('Changing chargeB of atom %s from %g to %g' % (latom.name, latom.qqB,newqB))
                    latom.qqB = newqB
                    qB_tot = sum([a.qqB for a in r.atoms])
                    TR('New total charge of B-state is %g' % round(qB_tot,5))
                else:
                    TR('No corrections applied to ensure integer charges')


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
                print('%6d %11s%7d%7s%7s%7d%11.6f%11.4f %11s%11.6f%11.4f' % \
                      (atom.id, atA, atom.resnr, atom.resname, atom.name, \
                       atom.cgnr, qqA, mA, atB, qqB, mB), file=fp)
                self.qA+=qqA
                self.qB+=qqB
            else:
                print('%6d %11s%7d%7s%7s%7d%11.6f%11.4f' % \
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
            if(len(p)==3):
                print('%6d %6d %6d' % (p[0].id, p[1].id, p[2]), file=fp)
            else:
                print('%6d %6d %6d %8s' % (p[0].id, p[1].id, p[2], p[3]), file=fp)

    def write_angles(self,fp, state='AB'):
        print('\n [ angles ]', file=fp)
        print(';  ai    aj    ak funct            c0            c1            c2            c3', file=fp)
        for ang in self.angles:
            if len(ang) == 4:
                print('%6d %6d %6d %6d' % (ang[0].id, ang[1].id, ang[2].id,ang[3]), file=fp)
            else:
                if state == 'A':
                    if ang[3]==1 :
                        print('%6d %6d %6d %6d %14.6f %14.6f' % (ang[0].id, ang[1].id, ang[2].id,ang[3],ang[4][0],ang[4][1]), file=fp)
                    elif ang[3]==5 :
                        if shape(ang[4])[0]==4:
                            print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % (ang[0].id, ang[1].id, ang[2].id,ang[3],ang[4][0],ang[4][1],ang[4][2],ang[4][3]), file=fp)
                        else:
                            print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % (ang[0].id, ang[1].id, ang[2].id,ang[3],ang[4][1],ang[4][2],ang[4][3],ang[4][4]), file=fp)
                    else :
                        print("Don't know how to print angletype %d" % ang[3])
                        exit()
                if state == 'AB':
#                   print ang[0].id, ang[1].id, ang[2].id
#                   print state
#                   print ang
                    #MS check type here, for charmm its different, Urey-Bradley
                    if ang[3]==1 :
                        # two possibilities: angle actually has a B-state or there is A state only
                        if(len(ang)>5): # B-state exists
                            print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                              (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1], \
                               ang[4][2], ang[5][1], ang[5][2], ang[0].name, ang[1].name, ang[2].name), file=fp)
                        else: # A-state only
                            print('%6d %6d %6d %6d %14.6f %14.6f' % (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][0], ang[4][1]), file=fp)

                    elif ang[3]==5:
                        # two possibilities: angle actually has a B-state or there is A state only
                        if(len(ang)>5): # B-state exists
                            print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                              (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1], \
                               ang[4][2], ang[4][3], ang[4][4], ang[5][1], \
                               ang[5][2], ang[5][3], ang[5][4], \
                               ang[0].name, ang[1].name, ang[2].name), file=fp)
                        else: # A-state only
                            print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1],ang[4][2], ang[4][3], ang[4][4] ), file=fp)
                    else :
                        print("Don't know how to print angletype %d" % ang[3])
                        exit()
                elif state == 'AA':
                    if ang[3]==1 :
                        print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1], \
                           ang[4][2], ang[4][1], ang[4][2], ang[0].name, ang[1].name, ang[2].name), file=fp)
                    elif ang[3]==5:
                        print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1], \
                           ang[4][2], ang[4][3], ang[4][4], ang[4][1], \
                           ang[4][2], ang[4][3], ang[4][4], \
                           ang[0].name, ang[1].name, ang[2].name), file=fp)
                    else :
                        print("Don't know how to print angletype %d" % ang[3])
                        exit()
                elif state == 'BB':
                    if ang[3]==1 :
                        print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[5][1], \
                           ang[5][2], ang[5][1], ang[5][2], ang[0].name, ang[1].name, ang[2].name), file=fp)
                    elif ang[3]==5:
                        print('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[5][1], \
                           ang[5][2], ang[5][3], ang[5][4], ang[5][1], \
                           ang[5][2], ang[5][3], ang[5][4], \
                           ang[0].name, ang[1].name, ang[2].name), file=fp)
                    else :
                        print("Don't know how to print angletype %d" % ang[3])
                        exit()

    def write_cmap(self, fp):
        print('\n [ cmap ]', file=fp)
        print(';  ai    aj    ak    al    am funct', file=fp)
        for d in self.cmap:
            print("%6d %6d %6d %6d %6d %4d" % ( d[0].id, d[1].id, d[2].id,d[3].id,d[4].id,d[5]), file=fp)

    def write_dihedrals(self, fp, state='AB'):
        print('\n [ dihedrals ]', file=fp)
        print(';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5', file=fp)
        for d in self.dihedrals:
            if len(d) == 5:
                print("%6d %6d %6d %6d %4d" % ( d[0].id, d[1].id, d[2].id, d[3].id, d[4]), file=fp)
            elif len(d) == 6:
                print("%6d %6d %6d %6d %4d %s" % ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], d[5]), file=fp)
            elif len(d) == 7:
                A, B = self.__check_case(d[:4])
                ast = d[5]
                bs = d[6]
                if ast == None or bs == None:
                    print(d[0].name, d[1].name, d[2].name, d[3].name, d[0].atomtype, d[1].atomtype, d[2].atomtype, d[3].atomtype, d[0].atomtypeB, d[1].atomtypeB, d[2].atomtypeB, d[3].atomtypeB)
                    print(d[0].type, d[1].type, d[2].type, d[3].type, d[0].typeB, d[1].typeB, d[2].typeB, d[3].typeB)
                if ast == 'NULL':
                    if d[4] == 3: # Ryckaert-Bellemans
                        ast = ' '.join(["%g" % x for x in [0,0,0,0,0,0]])
                    elif d[4] == 1 or d[4] == 4:
                        ast = ' '.join(["%g" % x for x in [0,0,0]])
                    elif d[4] == 9:
                        ast = ' '.join(["%g" % x for x in [0,0,0]])
                    elif d[4] == 2:
                        ast = ' '.join(["%g" % x for x in [0,0]])

                elif ast != 'NULL' and hasattr(ast,"append"):
                    ast = ' '.join(["%.10g" % x for x in d[5][1:]])
                if bs == 'NULL':
                    if d[4] == 3:
                        bs = ' '.join(["%g" % x for x in [0,0,0,0,0,0]])
                    elif d[4] == 1 or d[4] == 4:
                        bs = ' '.join(["%g" % x for x in [0,0,0]])
                    elif d[4] == 9:
                        bs = ' '.join(["%g" % x for x in [0,0,0]])
                    elif d[4] == 2:
                        bs = ' '.join(["%g" % x for x in [0,0]])

                elif bs !='NULL' and hasattr(bs,"append"):
                    bs = ' '.join(["%.10g" % x for x in d[6][1:]])
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

    def write_exclusions(self, fp):
        print('\n [ exclusions ]', file=fp)
        for excl in self.exclusions:
            towrite = ''
            for ex in excl:
                towrite += ' '+str(ex.id)
            print("%s" % towrite, file=fp)

    def write_vsites2(self, fp):
        print('\n [ virtual_sites2 ]', file=fp)
        print(';  ai    aj    ak  funct            c0            c1', file=fp)
        for vs in self.virtual_sites2:
            if len(vs) == 4:
                print("%6d %6d %6d %4d" % ( vs[0].id, vs[1].id, vs[2].id, vs[3]), file=fp)
            elif len(vs) == 5:
                print("%6d %6d %6d %4d %s" % ( vs[0].id, vs[1].id, vs[2].id, vs[3], vs[4]), file=fp)
            else:
                sys.stderr.write('EEK! Something went wrong while writing virtual_sites2!!!!\n')
                print(vs)
                sys.exit(1)


    def write_vsites3(self, fp):
        print('\n [ virtual_sites3 ]', file=fp)
        print(';  ai    aj    ak    al funct            c0            c1', file=fp)
        for vs in self.virtual_sites3:
            if len(vs) == 5:
                print("%6d %6d %6d %6d %4d" % ( vs[0].id, vs[1].id, vs[2].id, vs[3].id, vs[4]), file=fp)
            elif len(vs) == 6:
                print("%6d %6d %6d %6d %4d %s" % ( vs[0].id, vs[1].id, vs[2].id, vs[3].id, vs[4], vs[5]), file=fp)
            else:
                sys.stderr.write('EEK! Something went wrong while writing virtual_sites3!!!!\n')
                print(vs)
                sys.exit(1)

    def write_vsites4(self, fp):
        print('\n [ virtual_sites4 ]', file=fp)
        print(';  ai    aj    ak    al    am  funct            c0            c1          c2', file=fp)
        for vs in self.virtual_sites4:
            if len(vs) == 6:
                print("%6d %6d %6d %6d %6d %4d" % ( vs[0].id, vs[1].id, vs[2].id, vs[3].id, vs[4].id, vs[5]), file=fp)
            elif len(vs) == 7:
                print("%6d %6d %6d %6d %6d %4d %s" % ( vs[0].id, vs[1].id, vs[2].id, vs[3].id, vs[4].id, vs[5], vs[6]), file=fp)
            else:
                sys.stderr.write('EEK! Something went wrong while writing virtual_sites4!!!!\n')
                print(vs)
                sys.exit(1)

    def write_posre(self, fp):
        print('\n [ position_restraints ]', file=fp)
        print(';  ai    funct            c0            c1          c2', file=fp)
        for pr in self.posre:
            if len(pr) == 3:
                print("%6d %4d %s" % ( pr[0].id, pr[1], pr[2]), file=fp)
            else:
                sys.stderr.write('EEK! Something went wrong while writing position_restraints!!!!\n')
                print(pr)
                sys.exit(1)

    def write_ii(self, fp):
        fp.write('\n [ intermolecular_interactions ]\n')
        # bonds
        if 'bonds' in list(self.ii.keys()):
            fp.write(' [ bonds ]\n')
            for b in self.ii['bonds']:
                fp.write('%6d %6d %6d' % ( b[0].id, b[1].id, b[2] ))
                if len(b)>3:
                    for x in b[3]:
                        fp.write(' %14.6f' % x)
                fp.write('\n')
        # angles
        if 'angles' in list(self.ii.keys()):
            fp.write(' [ angles ]\n')
            for ang in self.ii['angles']:
                fp.write('%6d %6d %6d %6d' % ( ang[0].id, ang[1].id, ang[2].id, ang[3] ))
                if len(ang)>4:
                    for x in ang[4]:
                        fp.write(' %14.6f' % x)
                fp.write('\n')
        # dihedrals
        if 'dihedrals' in list(self.ii.keys()):
            fp.write(' [ dihedrals ]\n')
            for dih in self.ii['dihedrals']:
                fp.write('%6d %6d %6d %6d %6d' % ( dih[0].id, dih[1].id, dih[2].id, dih[3].id, dih[4] ))
                if len(dih)>5:
                    for x in dih[5]:
                        fp.write(' %14.6f' % x)
                fp.write('\n')

    def write_system(self,fp):
        print('[ system ]', file=fp)
        print(self.system, file=fp)

    def write_molecules(self,fp):
        print('[ molecules ]', file=fp)
        for mol, num in self.molecules:
            print("%s %d" % (mol,num), file=fp)


    #====================================================================================
    #    other functions

    def __check_case(self, atoms):
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


    def __atoms_morphe( self, atoms ):
        for atom in atoms:
            if atom.atomtypeB is not None and (atom.q!=atom.qB or atom.m != atom.mB or atom.atomtype != atom.atomtypeB): return True
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


#==================================================================================================
class ITPFile( TopolBase ):

    def __init__(self, filename ):
        TopolBase.__init__(self, filename)

#==================================================================================================
#==================================================================================================

class Topology( TopolBase ):

#    def __init__(self, filename, topfile = None, assign_types = True, cpp_path = [os.environ.get('GMXDATA')+'/top'], cpp_defs = [], version = 'old', ff = 'amber' ):
    def __init__(self, filename, topfile = None, assign_types = True, cpp_path = [os.environ.get('GMXLIB')], cpp_defs = [], version = 'old', ff = 'amber', ffpath=None ):
        TopolBase.__init__(self, filename, version)
        bItp = False
        if not topfile:
            topfile = filename
            bItp = True
        if assign_types:
            l = cpp_parse_file(topfile, cpp_defs = cpp_defs, cpp_path = cpp_path, itp = bItp, ffpath = ffpath)
            l = kickOutComments(l,'#')
            l = kickOutComments(l,';')
            self.BondedParams = BondedParser( l )
            self.NBParams = NBParser( l, version, ff = ff )
            self.assign_fftypes()


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


    def assign_fftypes(self):
        for atom in self.atoms:
            atom.type = self.NBParams.atomtypes[atom.atomtype]['bond_type']
            if atom.atomtypeB is not None:
                atom.typeB = self.NBParams.atomtypes[atom.atomtypeB]['bond_type']
            else:
                atom.typeB = atom.type

    def make_bond_params(self):
        for i, (at1,at2,func) in enumerate(self.bonds):
            param = self.BondedParams.get_bond_param(at1.type,at2.type)
            if param is None:
                print('Error! No bonded parameters found! (%s-%s)' % \
                      (at1.type, at2.type))
                sys.exit(1)
            self.bonds[i].append(param[1:])

    def make_angle_params(self):
        for i, (at1, at2, at3, func) in enumerate(self.angles):
            param = self.BondedParams.get_angle_param(at1.type, at2.type, at3.type)
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
                param = self.BondedParams.get_dihedral_param(at1.type, at2.type, \
                                                             at3.type, at4.type, \
                                                             func)
                if param is None:
                    print('Error! No dihedral parameters found! (%s-%s-%s-%s)' % \
                          (at1.type, at2.type, at3.type, at4.type))
                    print(func, dih)
                    sys.exit(1)
                del self.dihedrals[i][-1]
                self.dihedrals[i].append(param[1:])


#=================================================================================

class GAFFTopology( TopolBase ):

    def __init__(self, filename):
        TopolBase.__init__(self, filename)
        self.atomtypes = self.__read_atomtypes( filename )

    def __read_atomtypes( self, filename):
        l = open(filename).readlines()
        lst = readSection(l,'[ atomtypes ]','[')
        lst = parseList('ssffsff',lst)
        for line in lst:
            self.atomtypes[line[0]]=line[1:]

    def set_name( self, name ):
        self.name = name
        for atom in self.atoms:
            atom.name = name

#=================================================================================



class MDPError(Exception):
    def __init__(self, s):
        self.s = s
    def __str__(self):
        return repr(self.s)

#=================================================================================

class MDP:

    def __init__(self):
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
        ['nstcomm'                  , 1],
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
        ['coulombtype'              , 'PME'],
        ['rcoulomb-switch'          , 0],
        ['rcoulomb'                 , 1.2],
        ['epsilon-r'                , 1],
        ['epsilon_rf'               , 1],
        ['vdw-type'                 , 'switch'],
        ['rvdw-switch'              , 1],
        ['rvdw'                     , 1.1],
        ['DispCorr'                 , 'EnerPres'],
        ['table-extension'          , 1],
        ['energygrp_table'          ,''],
        ['fourierspacing'           , 0.14],
        ['fourier_nx'               , 0],
        ['fourier_ny'               , 0],
        ['fourier_nz'               , 0],
        ['pme_order'                , 4],
        ['ewald_rtol'               , 1e-05],
        ['ewald_geometry'           , '3d'],
        ['epsilon_surface'          , 0],
        ['optimize_fft'             , 'no'],
        ['implicit_solvent'         , 'No'],
        ['gb_algorithm'             , 'Still'],
        ['nstgbradii'               , 1],
        ['rgbradii'                 , 2],
        ['gb_epsilon_solvent'       , 80],
        ['gb_saltconc'              , 0],
        ['gb_obc_alpha'             , 1],
        ['gb_obc_beta'              , 0.8],
        ['gb_obc_gamma'             , 4.85],
        ['sa_surface_tension'       , 2.092],
        ['tcoupl'                   , 'v-rescale'],
        ['tc-grps'                  , ['Protein','non-protein']],
        ['tau-t'                    , [0.1, 0.1]],
        ['ref-t'                    , [298, 298]],
        ['Pcoupl'                   , 'Parrinello-Rahman'],
        ['Pcoupltype'               , 'Isotropic'],
        ['tau-p'                    , 1],
        ['compressibility'          , 4.6E-5],
        ['ref-p'                    , 1],
        ['refcoord_scaling'         , 'No'],
        ['andersen_seed'            , 815131],
        ['QMMM'                     , 'no'],
        ['QMMM-grps'                ,''],
        ['QMmethod'                 ,''],
        ['QMMMscheme'               , 'normal'],
        ['QMbasis'                  ,''],
        ['QMcharge'                 ,''],
        ['QMmult'                   ,''],
        ['SH'                       ,''],
        ['CASorbitals'              ,''],
        ['CASelectrons'             ,''],
        ['SAon'                     ,''],
        ['SAoff'                    ,''],
        ['SAsteps'                  ,''],
        ['MMChargeScaleFactor'      , 1],
        ['bOPT'                     ,''],
        ['bTS'                      ,''],
        ['annealing'                , ['no', 'no']],
        ['annealing_npoints'        , [2, 2]],
        ['annealing_time'           , [0, 50, 0, 50]],
        ['annealing_temp'           , [0, 298, 0, 298]],
        ['gen-vel'                  , 'no'],
        ['gen-temp'                 , 300],
        ['gen-seed'                 , 173529],
        ['constraints'              , 'all-bonds'],
        ['constraint-algorithm'     , 'Lincs'],
        ['continuation'             , 'yes'],
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
        ['free-energy'              , 'yes'],
        ['init-lambda'              , 0.5],
        ['delta-lambda'             , 0],
        ['sc-alpha'                 , 0.3],
        ['sc-power'                 , 1],
        ['sc-sigma'                 , 0.25],
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
                print('Warning! Ignoring entry \'%s\'' % key)
            else:
                if len(val) == 0:
                    self[key] = ''
                elif len(val) == 1:
                    self[key] = val[0]
                else:
                    self[key] = val
        return self


#=================================================================================

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
                    print('pmx_Warning_> No terminal oxygen atoms found in chain %s' % chain.id, file=sys.stderr)


#=================================================================================

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


#=================================================================================

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
    return bond_ene + angle_ene + dihedral_ene + improper_ene + nb_ene + lj14_ene + coul14_ene
