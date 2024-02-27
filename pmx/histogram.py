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
This file contains the Histogram class. I wrote it before realising
that numpy has a much fancier way to generate histograms.
I recommand to use the numpy histogram.
"""
import sys,os
from numpy import *

class Histogram:
    """ class to store data as histogram """

    def __init__(self,begin,end,incr):
        """ initialize histogram (start, end, increment)"""
        self.values=arange(begin,end,incr)
        self.counter={}
        for i in range(len(self.values)):
            self.counter[self.values[i]]=0

    def add(self,x,weight=1.):
        """ add a value """
        if x>self.values[-1]:
#            print "Value out of Histogram range!"
            pass
        else:
            for i in range(len(self.values)-1):
                if x>=self.values[i] and \
                       x<self.values[i+1]:
                    self.counter[self.values[i]]+=weight
                    break

    def write(self,filename=None):
        """ write histogram data to stdout or file"""
        if not filename:
            outfile=sys.stdout
        else:
            outfile=file(filename,'w')
        for value in self.values:
            outfile.write("%10.3f %10.3f\n" % (value,self.counter[value]))

    def norm(self):
        """ scale histogram such that the integral is 1"""
        list=[]
        for value in self.values:
            list.append(self.counter[value])
        x=trapz(list,x=self.values)
        for k in list(self.counter.keys()):
            self.counter[k]=self.counter[k]/x

    def integ(self,  min_val = None,  max_val = None):
        """ integrate histogram """
        list = []
        x = []
        for value in self.values:
            if not min_val and not max_val:
                list.append(self.counter[value])
                x.append( value)
            else:
                if min_val and value > min_val:
                    if not max_val:
                        list.append(self.counter[value])
                        x.append( value)
                    elif value < max_val:
                        list.append(self.counter[value])
                        x.append( value)
                if max_val and value < max_val:
                    if not min_val:
                        list.append(self.counter[value])
                        x.append( value)
                    elif value > min_val:
                        list.append(self.counter[value])
                        x.append( value)

        return trapz(list,x=x)


    def mean(self):
        """ calculate the mean value"""
        tot=0.
        nval=0
        for i in range(len(self.values)):
            tot+=self.values[i]*self.counter[self.values[i]]
            nval+=self.counter[self.values[i]]
        self.meanval=tot/nval
        return self.meanval

    def variance(self):
        """ calculate the variance """
        m=self.mean()
        var=0.
        nval=0
        for i in range(len(self.values)):
            nval+=self.counter[self.values[i]]
            var+=self.counter[self.values[i]]*(self.values[i]-m)**2
        self.var=var/nval
        return self.var

    def stddev(self):
        """ calculate std deviation """
        self.stddeviation=sqrt(self.variance())
        return self.stddeviation


if __name__=='__main__':
    print('testing histogram')
    h = Histogram(0,10,1)
    import random
    for i in range(1000):
        n = random.randint(0,10)
        h.add(n)
    print('mean = ', h.mean())
    print('stdev = ', h.stddev())
    print('var = ', h.variance())
    print('integ = ', h.integ())
    print('norming histogram....')
    h.norm()
    print('new integral', h.integ())
