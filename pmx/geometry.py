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
This file contains the Rotation class.
A Rotation instance is built with two vectors.
Usage:
     >>> v1 = [1,2,3]
     >>> v2 = [2,3,4]
     >>> r = Rotation(v1,v2)  # create rotation object around v2-v1
     >>> v3 = [4,5,6]
     >>> v3 = r.apply(v3)    # rotate v3 around v2-v1



"""
from numpy import *
from .atom import Atom
from  . import _pmx as _p

class Rotation2:

    def __init__(self,v1,v2):
        """ creates a rotation object
        around the vector v2-v1"""

        self.v1 = array(v1)
        self.v2 = array(v2)
        tmp = array(v2)
        self.vec = tmp-self.v1
        self.m1 = None
        self.m2 = None
        x = 1./linalg.norm(self.vec)
        self.norm_vec = self.vec*x
        self.__rm1()
        self.__rm2()

    def __rm1(self):

        a = self.norm_vec
        self.m1 = matrix( [
            [ a[0]*a[0], a[0]*a[1], a[0]*a[2]],
            [ a[1]*a[0], a[1]*a[1], a[1]*a[2]],
            [ a[2]*a[0], a[2]*a[1], a[2]*a[2]]
            ] )

    def __rm2(self):

        a = self.norm_vec
        self.m2 = matrix( [
            [ 0.0, -a[2], a[1]],
            [ a[2], 0.0, -a[0]],
            [ -a[1], a[0], 0.0]
            ] )

    def apply(self,v, phi):

        vec = v - self.v2
        b = dot(self.m1,vec)
        d = dot(self.m2,vec)
        a = cos(phi) * vec
        c = -cos(phi) * b
        e = sin(phi) * d
        vec  = a + b + c + e
        v = self.v2 + vec
        # ?
        return [x for x in v.getA()[0]]

class Rotation:

    def __init__(self,v1,v2):
        """ creates a rotation object
        around the vector v2-v1"""

        self.v1 = array(v1)
#        self.v2 = array(v2)
        self.v2 = [v2[0], v2[1], v2[2]] #array(v2)
        tmp = array(v2)
        self.vec = tmp-self.v1
        self.m1 = None
        self.m2 = None
        x = 1./linalg.norm(self.vec)
        self.norm_vec = self.vec*x
        self.__rm1()
        self.__rm2()

    def __rm1(self):

        a = self.norm_vec
        self.m1 =  [
            [ a[0]*a[0], a[0]*a[1], a[0]*a[2]],
            [ a[1]*a[0], a[1]*a[1], a[1]*a[2]],
            [ a[2]*a[0], a[2]*a[1], a[2]*a[2]]
            ]

    def __rm2(self):

        a = self.norm_vec
        self.m2 =  [
            [ 0.0, -a[2], a[1]],
            [ a[2], 0.0, -a[0]],
            [ -a[1], a[0], 0.0]
            ]

    def apply(self,v, phi):
        return _p.apply_rotation( self, [v[0], v[1], v[2]], phi)




def vec_ang(v1,v2):
    x1 = linalg.norm(v1)
    x2 = linalg.norm(v2)
    return arccos(inner(v1,v2)/(x1*x2))

def bb_super(mol1,mol2, use_orig_mc_coords = True):
    """ superpose mol2 on mol1"""

    N1,CA1,C1 = mol1.fetchm(['N','CA','C'])
    N2,CA2,C2 = mol2.fetchm(['N','CA','C'])

    if( (mol1.resname=='GLY') or (mol2.resname=='GLY') or (mol2.resname[:2]=='G2') ):
        fit_atoms( [N1,CA1,C1], [N2,CA2,C2], mol2.atoms )
    else:
        N1,CA1,C1,CB1 = mol1.fetchm(['N','CA','C','CB'])
        N2,CA2,C2,CB2 = mol2.fetchm(['N','CA','C','CB'])
        fit_atoms( [N1,CA1,C1,CB1], [N2,CA2,C2,CB2], mol2.atoms )

#    fit_atoms( [N1,CA1,C1], [N2,CA2,C2], mol2.atoms )

    if use_orig_mc_coords:
        atom_set = ['N','CA','C','H','O','HA','HN']
        gly_atom_set = ['N','CA','C','H','O','HA1','HN']
        if mol1.resname == 'GLY':
            atoms1 = mol1.fetchm(gly_atom_set)
        else:
            atoms1 = mol1.fetchm(atom_set)
        if mol2.resname == 'GLY' or mol2.resname[:2] == 'G2':
            atoms2 = mol2.fetchm(gly_atom_set)
        else:
            atoms2 = mol2.fetchm(atom_set)
        assert len(atoms1) == len(atoms2), "%s -> %s" % ( '-'.join( [a.name for a in atoms1]),'-'.join( [a.name for a in atoms2]) )
        for atom1, atom2 in zip(atoms1, atoms2):
            atom2.x = atom1.x

def nuc_super(mol1,mol2,name1=None,name2=None):
    """ superpose mol2 on mol1"""

    if name1==None:
        name1 = mol1.resname[:2]
    if name2==None:
        name2 = mol2.resname[:2]

    if name1 in ['DT','DC','RC','RU']:
        fit1_atoms = ['C1\'', 'C6','N1','C2','C5','N3']
    else:
        fit1_atoms = ['C1\'', 'C8','N9','C4','N7','C5']

    if name2 in ['DT','DC','RC','RU']:
        fit2_atoms = ['C1\'', 'C6','N1','C2','C5','N3']
    else:
        fit2_atoms = ['C1\'', 'C8','N9','C4','N7','C5']

    atoms1 = mol1.fetchm(fit1_atoms)
    atoms2 = mol2.fetchm(fit2_atoms)

    fit_atoms( atoms1, atoms2, mol2.atoms)


##     N1,CA1,C1 = mol1.fetchm(['O4\'','C1\'','C2\''])
##     N2,CA2,C2 = mol2.fetchm(['O4\'','C1\'','C2\''])

##     fit_atoms( [N1,CA1,C1], [N2,CA2,C2], mol2.atoms )

##     rotN1 = None
##     rotC1 = None
##     rotN2 = None
##     rotC2 = None

##     if mol1.resname in ['DT','DC','RC','RU']:
##         rotN1, rotC1 = mol1.fetchm(['N1','C2'])
##     elif mol1.resname in ['DG','DA','RG','RA']:
##         rotN1, rotC1 = mol1.fetchm(['N9','C4'])
##     if mol2.resname in ['DTG','DTA','DCA','DCG',
##                         'DTC','DCT','RUC','RCU',
##                         'RUG','RUA','RCA','RCG']:
##         rotN2, rotC2 = mol2.fetchm(['N1','C2'])
##     elif mol2.resname in ['DGT','DAT','DAC','DGC',
##                           'DGA','DAG','RGA','RAG',
##                           'RGU','RAU','RAC','RGC']:
##         rotN2, rotC2 = mol2.fetchm(['N9','C4'])
##     if rotN2 is not None and rotC2 is not None:
##         print rotN1.name, rotC1.name, rotN2.name, rotC2.name
##         v1 = array(rotN1.x)-array(rotC1.x)
##         v2 = array(rotN2.x)-array(rotC2.x)
##         ang = vec_ang(v1,v2)
##         v3 = array(rotN1.x)-array(CA1.x)
## #        print 'ang2 = ' , ang, ang*180/pi
##         if linalg.det(matrix([v1,v3,v2])) < 0:
##  #           print 'DET =', linalg.det(matrix([v1,v3,v2]))
##             ang *= -1
##         r = Rotation(CA1.x, rotN1.x)
##         for atom in mol2.atoms:
##             atom.x = r.apply(atom.x,ang)
##         v1 = array(rotN1.x)-array(rotC1.x)
##         v2 = array(rotN2.x)-array(rotC2.x)
##         ang = vec_ang(v1,v2)
## #        print 'ang2 = ' , ang, ang*180/pi
## #    else:
## #        print 'No correction applied'
## #

def planarity(atom_list):

    coords = [a.x for a in atom_list]
    plan = _p.planarity(coords)
    return plan

def apply_fit_R( atoms, R):

    for atom in atoms:
        x_old = [x for x in atom.x]
        for r in range(3):
            atom.x[r] = 0
            for c in range(3):
                atom.x[r]+=R[r][c]*x_old[c]

def center_vector( v ):
    vout = _p.center_vec( v )
    return( vout )

def calc_fit_R( cs1, cs2, m ):
    R = _p.calc_fit_R(cs1, cs2, m)
    return(R)


def fit(model1, model2, atom_names = []):
    if atom_names:
        subset1 = model1.fetch_atoms( atom_names )
        subset2 = model2.fetch_atoms( atom_names )
        cs1 = [a.x for a in subset1]
        cs2 = [a.x for a in subset2]
    else:
        cs1 = model1.coords()
        cs2 = model2.coords()

    assert( len(cs1) == len(cs2) )
    m = [1. for x in cs1] # dummy array
    v = _p.center_vec( cs1 )
    v2 = _p.center_vec( cs2 )
    R = _p.calc_fit_R(cs1, cs2, m)
    model2.translate( [-v2[0], -v2[1], -v2[2] ] )
    apply_fit_R( model2.atoms, R)
    model2.translate( v )


def fit_by_ndx(ref, model, ndx1, ndx2):
    crd1 = [ref.atoms[i-1].x for i in ndx1]
    crd2 = [model.atoms[i-1].x for i in ndx2]

    assert( len(crd1) == len(crd2) )
    m = [1. for x in crd1] # dummy array
    v = _p.center_vec( crd1 )
    v2 = _p.center_vec( crd2 )
    R = _p.calc_fit_R(crd1, crd2, m)
    model.translate( [-v2[0], -v2[1], -v2[2] ] )
    apply_fit_R( model.atoms, R)
    model.translate( v )

def translate_by_ndx(struct, ndx):
    crd = [struct.atoms[i-1].x for i in ndx]
    m = [1. for x in crd]
    v = _p.center_vec( crd )
    struct.translate( [-v[0], -v[1], -v[2]] )
    return(v)

def fit_atoms( fit_atoms1, fit_atoms2, rot_atoms2 ):

    cs1 = [a.x for a in fit_atoms1]
    cs2 = [a.x for a in fit_atoms2]
    assert len(cs1) == len(cs2)
    m = [1. for x in cs1] # dummy array
    v = _p.center_vec( cs1 )
    v2 = _p.center_vec( cs2 )
    R = _p.calc_fit_R(cs1, cs2, m)
    for atom in rot_atoms2:
        atom.x[0]-=v2[0]
        atom.x[1]-=v2[1]
        atom.x[2]-=v2[2]
    apply_fit_R( rot_atoms2, R)
    for atom in rot_atoms2:
        atom.x[0]+=v[0]
        atom.x[1]+=v[1]
        atom.x[2]+=v[2]
