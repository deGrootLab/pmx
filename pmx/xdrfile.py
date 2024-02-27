#  -*- mode: python; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*-
#
#  $Id$
#
#  Copyright (c) Erik Lindahl, David van der Spoel 2003-2007.
#  Coordinate compression (c) by Frans van Hoesel.
#  Python wrapper (c) by Roland Schulz
#
#  IN contrast to the rest of Gromacs, XDRFILE is distributed under the
#  BSD license, so you can use it any way you wish, including closed source:
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#
#  Adapted by Daniel Seeliger for use in the pmx package (Aug 2015)
#
import numpy as np
from ctypes import *
import os.path

mTrr,mNumPy=1,2
auto_mode=0
out_mode=42


class Frame:
    #variables
    #x: rvec*natoms / numpy array if installed
    #box DIM*DIM
    #step
    #time
    #prec
    #lam: lambda

    def __init__(self,n,mode,x=None,box=None,units=None,v=None,f=None):
        #create vector for x
        self.natoms = n
        # x (coordinates)
        if mode==out_mode:
            scale = 1.0
            if units == 'A':
                scale=0.1
            self.x=((c_float*3)*n)()
            i=0
            for a in range(0,self.natoms):
                for dim in range(0,3):
                    self.x[a][dim] = scale*x[i]
                    i+=1
        elif mode&mNumPy and mode!=out_mode:
            self.x=empty((n,3),dtype=float32)
        else:
            self.x=((c_float*3)*n)()

        # dummy v and f for .trr
        self.v_size = c_size_t(0)
        self.f_size = c_size_t(0)
        if mode&mNumPy and mode!=out_mode:
            self.v=empty((n,3),dtype=float32)
            self.f=empty((n,3),dtype=float32)
        else:
            self.v=c_size_t(0)#((c_float*3)*n)()
            self.f=c_size_t(0)#((c_float*3)*n)()

        # box
        if box!=None:
            self.box = (c_float*3*3)()
            for r in range(0,3):
                for c in range(0,3):
                    self.box[r][c] = box[r][c]
        elif mode&mNumPy and mode!=out_mode:
            self.box = empty((3,3),float32)
        else:
            self.box = (c_float*3*3)()


    def update_box( self, box ):
        for i in range(3):
            for k in range(3):
                box[i][k] = self.box[i][k]

    def update_atoms( self, atom_sel ):
        for i, atom in enumerate(atom_sel.atoms):
            if atom_sel.unity == 'A':
                atom.x[0] = self.x[i][0]*10
                atom.x[1] = self.x[i][1]*10
                atom.x[2] = self.x[i][2]*10
            else:
                atom.x[0] = self.x[i][0]
                atom.x[1] = self.x[i][1]
                atom.x[2] = self.x[i][2]

    def update( self, atom_sel ):
        self.update_atoms(atom_sel )
        self.update_box( atom_sel.box )


    def get_natoms(self):
        return self.natoms
    def get_time(self):
        return self.time
    def get_step(self):
        return self.step
    def get_prec(self):
        return self.prec
        def get_box(self, mode='std'):
            if mode == 'std':
                box = [[0.,0.,0.],
                           [0.,0.,0.],
                           [0.,0.,0.]]
            elif mode == 'numpy':
                box = np.zeros((3,3))
            for i in range(3):
                for k in range(3):
                    box[i][k] = self.box[i][k]
            return box



    def __str__(self):

        s = '< xdrlib.Frame: natoms = %d | step = %d | time = %8.2f >' % (self.get_natoms(), self.get_step(), self.get_time())
        return s




class XDRFile:
    exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE, exdrINT, exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC, exdrNOMEM, exdrENDOFFILE, exdrNR = list(range(13))

    #
    def __init__(self,fn,mode="Auto",ft="Auto",atomNum=False):
        if mode=="NumPy":
            self.mode=mNumPy
            try:
                empty
            except NameError:
                raise IOError("NumPy selected but not correctly installed")
        elif mode=="Std":
            self.mode=0
        elif mode=="Auto":
            self.mode=auto_mode
        elif mode=='Out':
            self.mode=out_mode
        else:
            raise IOError("unsupported mode")

        if ft=="Auto":
            ft = os.path.splitext(fn)[1][1:]

        if self.mode!=out_mode:
            if ft=="trr":
                self.mode|=mTrr
            elif ft=="xtc":
                pass
            else:
                raise IOError("Only xtc and trr supported")

        #load libxdrfil
        try:
            p = os.path.join(os.path.dirname(__file__),'_xdrio.cpython-310-x86_64-linux-gnu.so')
            self.xdr=cdll.LoadLibrary(p)
        except:
            raise IOError("_xdrio.so can't be loaded")


        #open file
        if self.mode==out_mode:
            self.xd = self.xdr.xdrfile_open(fn,"w")
        else:
            self.xd = self.xdr.xdrfile_open(fn,"r")
        if not self.xd: raise IOError("Cannot open file: '%s'"%fn)

        #read natoms
        natoms=c_int()
        if self.mode==out_mode:
            if atomNum==False:
                raise Exception("To write an .xtc need to provide the number of atoms")
            self.natoms = atomNum
        else:
            if self.mode&mTrr:
                r=self.xdr.read_trr_natoms(fn,byref(natoms))
            else:
                r=self.xdr.read_xtc_natoms(fn,byref(natoms))
            if r!=self.exdrOK: raise IOError("Error reading: '%s'"%fn)
            self.natoms=natoms.value

        #for NumPy define argtypes - ndpointer is not automatically converted to POINTER(c_float)
        #alternative of ctypes.data_as(POINTER(c_float)) requires two version for numpy and c_float array
        if self.mode&mNumPy and self.mode!=out_mode:
            self.xdr.read_xtc.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float),
              ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),POINTER(c_float)]
            self.xdr.read_trr.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float),POINTER(c_float),
              ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),
              POINTER(c_float),POINTER(c_float)]

    def write_xtc_frame( self, step=0, time=0.0, prec=1000.0, lam=0.0, box=False, x=False, units='A', bTrr=False ):
        f = Frame(self.natoms,self.mode,box=box,x=x,units=units)
        step = c_int(step)
        time = c_float(time)
        prec = c_float(prec)
        lam = c_float(lam)
        if bTrr:
            result = self.xdr.write_trr(self.xd,self.natoms,step,time,lam,f.box,f.x,f.v,f.f)
        else:
            result = self.xdr.write_xtc(self.xd,self.natoms,step,time,f.box,f.x,prec)

    def __iter__(self):
        f = Frame(self.natoms,self.mode)
        #temporary c_type variables (frame variables are python type)
        step = c_int()
        time = c_float()
        prec = c_float()
        lam = c_float()
        if self.mode!=out_mode:
            while 1:
                #read next frame
                if not self.mode&mTrr:
                    result = self.xdr.read_xtc(self.xd,self.natoms,byref(step),byref(time),f.box,f.x,byref(prec))
                    f.prec=prec.value
                else:
                    result = self.xdr.read_trr(self.xd,self.natoms,byref(step),byref(time),byref(lam),f.box,f.x,None,None) #TODO: make v,f possible
                    f.lam=lam.value

                #check return value
                if result==self.exdrENDOFFILE: break
                if result==self.exdrINT and self.mode&mTrr:
                    break  #TODO: dirty hack. read_trr return exdrINT not exdrENDOFFILE
                if result!=self.exdrOK: raise IOError("Error reading xdr file")

                #convert c_type to python
                f.step=step.value
                f.time=time.value
                yield f
