#!/usr/bin/env python
# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2017 by Daniel Seeliger
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


from pmx.parser import read_and_format
from pmx.estimators import Jarz, JarzGauss, Crooks, BAR, data2gauss, ks_norm_test
import sys
import os
import time
import re
import numpy as np
from scipy.integrate import simps
import pickle
import argparse
from .cli import check_unknown_cmd

# Constants
kb = 0.00831447215   # kJ/(K*mol)


# ==============================================================================
#                               FUNCTIONS
# ==============================================================================
def gauss_func(A, mean, dev, x):
    '''Given the parameters of a Gaussian and a range of the x-values, returns
    the y-values of the Gaussian function'''
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2.)/(2.0*(dev**2.))))
    return y


# -------------
# Files Parsing
# -------------
def _longest_dgdl_file( lst ):
    '''Takes a list of dgdl.xvg files and returns the id (starting from 0) of
    of the longest file.

    Parameters
    ----------
    lst : list
        list containing the paths to the dgdl.xvg files.

    Returns
    -------
    ind : int
        index of the longest file
    '''

    ind = 0
    maxlen = -9999
    for idx, f in enumerate(lst[0:]):
        try:
            fp = open(f)
            lines = fp.readlines()
            if len(lines)>maxlen:
                maxlen = len(lines)
                ind = idx
        except:
            pass
    return(ind)


def parse_dgdl_files(lst, lambda0=0, invert_values=False):
    '''Takes a list of dgdl.xvg files and returns the integrated work values

    Parameters
    ----------
    lst : list
        list containing the paths to the dgdl.xvg files.
    lambda0 : [0,1]
        whether the simulations started from lambda 0 or 1. Default is 0.
    invert_values : bool
        whether to invert the sign of the returned work value.

    Returns
    -------
    w : array
        array of work values.
    lst : list
        sorted list of input dgdl.avg files corresponding to the work values
        in w.
    '''

    # check lambda0 is either 0 or 1
    assert lambda0 in [0, 1]

    # identify file with the most entries
    imax = _longest_dgdl_file( lst )

    _check_dgdl(lst[imax], lambda0)
    first_w, ndata = integrate_dgdl(lst[imax], lambda0=lambda0,
                                    invert_values=invert_values)
    w_list = [first_w]
    for idx, f in enumerate(lst[0:]):
        if idx==imax: # already processed
            continue
        sys.stdout.write('\r    Reading %s' % f)
        sys.stdout.flush()

        w, _ = integrate_dgdl(f, ndata=ndata, lambda0=lambda0,
                              invert_values=invert_values)
        if w is not None:
            w_list.append(w)

    print('\n')

    return w_list


def integrate_dgdl(fn, ndata=-1, lambda0=0, invert_values=False):
    '''Integrates the data in a dgdl.xvg file.

    Parameters
    ----------
    fn : str
        the inpur dgdl.xvg file from Gromacs.
    ndata : int, optional
        number of datapoints in file. If -1, then ??? default is -1.
    lambda0 : [0,1]
        whether the simulations started from lambda 0 or 1. Default is 0.
    invert_values : bool
        whether to invert the sign of the returned work value.

    Returns
    -------
    integr : float
        result of the integration performed using Simpson's rule.
    ndata : int
        number of data points in the input file.
    '''

    # check lambda0 is either 0 or 1
    assert lambda0 in [0, 1]

    try:
        lines = open(fn).readlines()
    except:
        return None, None
    if not lines:
        return None, None

    # extract dgdl datapoints into r
    # TODO: we removed the check for file integrity. We could have an
    # optional files integrity check before calling this integration func

    lines = [l for l in lines if l[0] not in '#@&']
    try:
        r = [float(x.split()[1]) for x in lines]
    except:
        r = "incomplete_file"
        ndata = 1

    if ndata != -1 and len(r) != ndata:
        try:
            if "incomplete_file" in r:
                print(' !! Skipping %s (incomplete file, probably simulation crashed)\n' % fn)
            else:
                print(' !! Skipping %s ( read %d data points, should be %d )' % (fn, len(r), ndata))
        except:
            print(' !! Skipping %s ' % (fn))
        return None, None
    # convert time to lambda
    ndata = len(r)
    dlambda = 1./float(ndata)
    if lambda0 == 1:
        dlambda *= -1

    # arrays for the integration
    # --------------------------
    # array of lambda values
    x = [lambda0+i*dlambda for i, dgdl in enumerate(r)]
    # array of dgdl
    y = r

    if lambda0 == 1:
        x.reverse()
        y.reverse()

    if invert_values is True:
        integr = simps(y, x) * (-1)
        return integr, ndata
    else:
        integr = simps(y, x)
        return integr, ndata


def _check_dgdl(fn, lambda0):
    '''Prints some info about a dgdl.xvg file.'''
    l = open(fn).readlines()
    if not l:
        return None
    r = []
    for line in l:
        if line[0] not in '#@&':
            r.append([float(x) for x in line.split()])
    ndata = len(r)
    dlambda = 1./float(ndata)
    if lambda0 == 1:
        dlambda *= -1

    print('    # data points: %d' % ndata)
    print('    Length of trajectory: %8.3f ps' % r[-1][0])
    print('    Delta lambda: %8.5f' % dlambda)


def _dump_integ_file(outfn, f_lst, w_lst):
    with open(outfn, 'w') as f:
        for fn, w in zip(f_lst, w_lst):
            f.write('{dhdl} {work}\n'.format(dhdl=fn, work=w))


def _data_from_file(fn):
    data = read_and_format(fn, 'sf')
    return [a[1] for a in data]


# ------------------
# Plotting functions
# ------------------
def plot_work_dist(wf, wr, fname='Wdist.png', nbins=20, dG=None, dGerr=None,
                   units='kJ/mol', dpi=300, statesProvided='AB'):

    from matplotlib import pyplot as plt

    '''Plots forward and reverse work distributions. Optionally, it adds the
    estimate of the free energy change and its uncertainty on the plot.

    Parameters
    ----------
    wf : list
        list of forward work values.
    wr : list
        list of reverse work values.
    fname : str, optional
        filename of the saved image. Default is 'Wdist.png'.
    nbins : int, optional
        number of bins to use for the histogram. Default is 20.
    dG : float, optional
        free energy estimate.
    dGerr : float, optional
        uncertainty of the free energy estimate.
    units : str, optional
        the units of dG and dGerr. Default is 'kJ/mol'.
    dpi : int
        resolution of the saved image file.
    statesProvided: str
        work values for two states or only one

    Returns
    -------
    None

    '''

    def smooth(x, window_len=11, window='hanning'):

        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
            raise ValueError("Input vector needs to be bigger than "
                             "window size.")
        if window_len < 3:
            return x
        if window not in ['flat', 'hanning', 'hamming',
                          'bartlett', 'blackman']:
            raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                             "'bartlett', 'blackman'")
        s = np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
        # moving average
        if window == 'flat':
            w = np.ones(window_len, 'd')
        else:
            w = eval('np.' + window + '(window_len)')
        y = np.convolve(w/w.sum(), s, mode='same')
        return y[window_len-1:-window_len+1]

    plt.figure(figsize=(8, 6))
    plt.subplot(1, 2, 1)
    x1 = 0
    x2 = 0
    if 'A' in statesProvided:
        x1 = list(range(len(wf)))
    if 'B' in statesProvided:
        x2 = list(range(len(wr)))
    if x1 > x2:
        x = x1
    else:
        x = x2
    if 'A' in statesProvided:
        mf, devf, Af = data2gauss(wf)
    if 'B' in statesProvided:
        mb, devb, Ab = data2gauss(wr)

    if 'AB' in statesProvided:
        mini = min(wf+wr)
        maxi = max(wf+wr)
        sm1 = smooth(np.array(wf))
        sm2 = smooth(np.array(wr))
        plt.plot(x1, wf, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
        plt.plot(x1, sm1, 'g-', linewidth=3)
        plt.plot(x2, wr, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
        plt.plot(x2, sm2, 'b-', linewidth=3)
    elif 'A' in statesProvided:
        maxi = max(wf)
        mini = min(wf)
        sm1 = smooth(np.array(wf))
        plt.plot(x1, wf, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
        plt.plot(x1, sm1, 'g-', linewidth=3)
    elif 'B' in statesProvided:
        maxi = max(wr)
        mini = min(wr)
        sm2 = smooth(np.array(wr))
        plt.plot(x2, wr, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
        plt.plot(x2, sm2, 'b-', linewidth=3)

    plt.legend(shadow=True, fancybox=True, loc='upper center',
               prop={'size': 12})
    plt.ylabel(r'W [kJ/mol]', fontsize=20)
    plt.xlabel(r'# Snapshot', fontsize=20)
    plt.grid(lw=2)
    plt.xlim(0, x[-1]+1)
    xl = plt.gca()
    for val in list(xl.spines.values()):
        val.set_lw(2)
    plt.subplot(1, 2, 2)
    plt.hist(wf, bins=nbins, orientation='horizontal', facecolor='green',
             alpha=.75, density=True)
    plt.hist(wr, bins=nbins, orientation='horizontal', facecolor='blue',
             alpha=.75, density=True)

    x = np.arange(mini, maxi, .5)

    if 'AB' in statesProvided:
        y1 = gauss_func(Af, mf, devf, x)
        y2 = gauss_func(Ab, mb, devb, x)
        plt.plot(y1, x, 'g--', linewidth=2)
        plt.plot(y2, x, 'b--', linewidth=2)
        size = max([max(y1), max(y2)])
    elif 'A' in statesProvided:
        y1 = gauss_func(Af, mf, devf, x)
        plt.plot(y1, x, 'g--', linewidth=2)
        size = max(y1)
    elif 'B' in statesProvided:
        y2 = gauss_func(Ab, mb, devb, x)
        plt.plot(y2, x, 'b--', linewidth=2)
        size = max(y2)

    res_x = [dG, dG]
    res_y = [0, size*1.2]
    if dG is not None and dGerr is not None:
        plt.plot(res_y, res_x, 'k--', linewidth=2,
                 label=r'$\Delta$G = %.2f $\pm$ %.2f %s' % (dG, dGerr, units))
        plt.legend(shadow=True, fancybox=True, loc='upper center',
                   prop={'size': 12})
    elif dG is not None and dGerr is None:
        plt.plot(res_y, res_x, 'k--', linewidth=2,
                 label=r'$\Delta$G = %.2f %s' % (dG, units))
        plt.legend(shadow=True, fancybox=True, loc='upper center',
                   prop={'size': 12})
    else:
        plt.plot(res_y, res_x, 'k--', linewidth=2)

    plt.xticks([])
    plt.yticks([])
    xl = plt.gca()
    for val in list(xl.spines.values()):
        val.set_lw(2)
    plt.subplots_adjust(wspace=0.0, hspace=0.1)
    plt.savefig(fname, dpi=dpi)


# ---------------------
# Some helper functions
# ---------------------
def _tee(fp, s):
    print(s, file=fp)
    print(s)


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def time_stats(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return h, m, s


# ==============================================================================
#                      COMMAND LINE OPTIONS AND MAIN
# ==============================================================================
def parse_options():

    parser = argparse.ArgumentParser(description='Calculates free energies '
            'from fast growth thermodynamic integration simulations. '
            'Available methods for free energy estimation: '
            'Crooks Gaussian Intersection (CGI); '
            'Benett Acceptance Ratio (BAR); '
            'Jarzinski equality (JARZ);'
            'Jarzinski with a Gaussian approximatiuon (JARZ_Gauss).',)

    exclus = parser.add_mutually_exclusive_group()

    parser.add_argument('-fA',
                        metavar='dgdl',
                        dest='filesAB',
                        type=str,
                        help='dgdl.xvg files for the A->B simulations. Use '
                        'wildcard to select multiple xvg files: e.g. "-fA '
                        './forward_results/dgdl*.xvg"',
                        nargs='+')
    parser.add_argument('-fB',
                        metavar='dgdl',
                        dest='filesBA',
                        type=str,
                        help='dgdl.xvg files for the B->A simulations Use '
                        'wildcard to select multiple xvg files: e.g. "-fB '
                        './backward_results/dgdl*.xvg"',
                        nargs='+')
    parser.add_argument('-m',
                        metavar='method',
                        type=str.lower,
                        dest='methods',
                        help='Choose one or more estimators to use from the '
                        'available ones: CGI, BAR, JARZ. Default is all.',
                        default=['cgi', 'bar', 'jarz'],
                        nargs='+')
    parser.add_argument('-t',
                        metavar='temperature',
                        dest='temperature',
                        type=float,
                        help='Temperature in Kelvin. Default is 298.15.',
                        default=298.15)
    parser.add_argument('-o',
                        metavar='result file',
                        dest='outfn',
                        type=str,
                        help='Filename of output result file. Default is '
                        '"results.txt."',
                        default='results.txt')
    parser.add_argument('-b',
                        metavar='nboots',
                        dest='nboots',
                        type=int,
                        help='Number of bootstrap samples to use for the '
                        'bootstrap estimate of the standard errors. Default '
                        'is 0 (no bootstrap).',
                        default=100)
    parser.add_argument('-n',
                        metavar='nblocks',
                        dest='nblocks',
                        type=int,
                        help='Number of blocks to divide the data into for '
                        'an estimate of the standard error. You can use this '
                        'when multiple independent equilibrium simulations'
                        'have been run so to estimate the error from the '
                        'repeats. Default is 1 (i.e. no repeats). It assumes '
                        'the dhdl files for each repeat are read in order and '
                        'are contiguous, e.g. dhdl_0 to dhdl_9 is the first '
                        'repeat, dhdl_10 to dhdl_19 is the second one, etc.',
                        default=1)
    parser.add_argument('--integ_only',
                        dest='integ_only',
                        help='Whether to do integration only; the integrated '
                        'values are computed and saved, and the program '
                        'terminated. Default is False.',
                        default=False,
                        action='store_true')
    parser.add_argument('-iA',
                        metavar='work input',
                        dest='iA',
                        type=str,
                        help='Two-column dat file containing the list of input'
                        ' files and their respective integrated work values '
                        'for the forward (A->B) tranformation.')
    parser.add_argument('-iB',
                        metavar='work input',
                        dest='iB',
                        type=str,
                        help='Two-column dat file containing the list of input'
                        ' files and their respective integrated work values '
                        'for the reverse (B->A) tranformation.')
    parser.add_argument('-oA',
                        metavar='work output',
                        dest='oA',
                        type=str,
                        help='File where to save the list of input dhdl'
                        ' files and their respective integrated work values '
                        'for the forward (A->B) tranformation. Default is '
                        '"integA.dat"',
                        default='integA.dat')
    parser.add_argument('-oB',
                        metavar='work output',
                        dest='oB',
                        type=str,
                        help='File where to save the list of input dhdl'
                        ' files and their respective integrated work values '
                        'for the reverse (B->A) tranformation. Default is '
                        '"integB.dat"',
                        default='integB.dat')
    parser.add_argument('--reverseB',
                        dest='reverseB',
                        help='Whether to reverse the work values for the '
                        'backward (B->A) transformation. This is useful '
                        'when in Gromacs both forward and reverse simulations '
                        'were run from lambda zero to one.'
                        'Default is False.',
                        default=False,
                        action='store_true')
    # The following are mutually exclusive options
    exclus.add_argument('--skip',
                        metavar='',
                        dest='skip',
                        type=int,
                        help='Skip files, i.e. pick every nth work value. '
                        'Default is 1 (all); with 2, every other work value '
                        'is discarded, etc.',
                        default=1)
    exclus.add_argument('--slice',
                        metavar='',
                        dest='slice',
                        type=int,
                        help='Subset of trajectories to analyze.'
                        'Provide list slice, e.g. "10 50" will'
                        ' result in selecting dhdl_files[10:50].'
                        ' Default is all.',
                        default=None,
                        nargs=2)
    exclus.add_argument('--rand',
                        metavar='',
                        dest='rand',
                        type=int,
                        help='Take a random subset of trajectories. '
                        'Default is None (do not take random subset)',
                        default=None)
    exclus.add_argument('--index',
                        metavar='',
                        dest='index',
                        type=int,
                        help='Zero-based index of files to analyze (e.g.'
                        ' 0 10 20 50 60). It keeps '
                        'the dhdl.xvg files according to their position in the'
                        ' list, sorted according to the filenames. Default '
                        'is None (i.e. all dhdl are used).',
                        default=None,
                        nargs='+')
    parser.add_argument('--prec',
                        metavar='',
                        dest='precision',
                        type=int,
                        help='The decimal precision of the screen/file output.'
                        ' Default is 2.',
                        default=2)
    parser.add_argument('--units',
                        metavar='',
                        dest='units',
                        type=str.lower,
                        help='The units of the output. Choose from "kJ", '
                        '"kcal", "kT". Default is "kJ."',
                        default='kJ',
                        choices=['kj', 'kcal', 'kt'])
    parser.add_argument('--pickle',
                        dest='pickle',
                        help='Whether to save the free energy results from '
                        'the estimators in pickled files. Default is False.',
                        default=False,
                        action='store_true')
    parser.add_argument('--no_ks',
                        dest='do_ks_test',
                        help='Whether to do a Kolmogorov-Smirnov test '
                        'to check whether the Gaussian assumption for CGI '
                        'holds. Default is True; this flag turns it to False.',
                        default=True,
                        action='store_false')
    parser.add_argument('--work_plot',
                        metavar='',
                        dest='wplot',
                        type=str,
                        help='Name of image file showing the distribution of '
                        'work values. Default is "wplot.png". If you want to '
                        'avoid saving this plot, pass "none" to this flag. '
                        'If you choose to calculate the free energy with '
                        'multiple estimators, the dG values shown on the plot '
                        'will be chosen following the hierarchy '
                        'BAR > CGI > JARZ.',
                        default='wplot.png')
    parser.add_argument('--nbins',
                        metavar='',
                        dest='nbins',
                        type=int,
                        help='Number of histograms bins for the plot. '
                        'Default is 10.',
                        default=10)
    parser.add_argument('--dpi',
                        metavar='',
                        dest='dpi',
                        type=int,
                        help='Resolution of the plot. Default is 300.',
                        default=300)

    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    from pmx import __version__
    args.pmx_version = __version__

    return args


def main(args):
    """Run the main script.

    Parameters
    ----------
    args : argparse.Namespace
        The command line arguments
    """

    # start timing
    stime = time.time()

    # input arguments
    filesAB = []
    filesBA = []
    statesProvided = 'AB'
    out = open(args.outfn, 'w')
    if (args.iA is None) and (args.iB is None):
        if (args.filesAB is None) and (args.filesBA is None):
            exit('Need to provide dhdl.xvg files or integrated work values')
        elif args.filesAB is None:
            statesProvided = 'B'
            _tee(out, 'Only one directional Jarzynski estimator will be used')
            filesBA = natural_sort(args.filesBA)
        elif args.filesBA is None:
            statesProvided = 'A'
            _tee(out, 'Only one directional Jarzynski estimator will be used')
            filesAB = natural_sort(args.filesAB)
        else:
            filesAB = natural_sort(args.filesAB)
            filesBA = natural_sort(args.filesBA)
    elif (args.iA is None) or (args.iB is None):
        if args.iA is None:
            statesProvided = 'B'
        else:
            statesProvided = 'A'
        _tee(out, 'Only one directional Jarzynski estimator will be used')
    T = args.temperature
    skip = args.skip
    prec = args.precision
    methods = args.methods
    reverseB = args.reverseB
    integ_only = args.integ_only
    nboots = args.nboots
    nblocks = args.nblocks
    do_ks_test = args.do_ks_test

    # -------------------
    # Select output units
    # -------------------
    units = args.units
    if units.lower() == 'kj':
        # kJ is the input from GMX
        unit_fact = 1.
        units = 'kJ/mol'
    elif units == 'kcal':
        unit_fact = 1./4.184
        units = 'kcal/mol'
    elif units.lower() == 'kt':
        unit_fact = 1./(kb*T)
        units = 'kT'
    else:
        exit('No unit type \'%s\' available' % units)

    print("# analyze_crooks.py, pmx version = %s" % args.pmx_version, file=out)
    print("# pwd = %s" % os.getcwd(), file=out)
    print("# %s (%s)" % (time.asctime(), os.environ.get('USER')), file=out)
    print("# command = %s" % ' '.join(sys.argv), file=out)
    _tee(out, "\n")

    # ==========
    # Parse Data
    # ==========

    # If list of dgdl.xvg files are provided, parse dgdl
    if args.iA is None and args.iB is None:
        # If random selection is chosen, do this before reading files and
        # calculating the work values.
        if args.rand is not None:
            if 'A' in statesProvided:
                filesAB = np.random.choice(filesAB, size=args.rand, replace=False)
            if 'B' in statesProvided:
                filesBA = np.random.choice(filesBA, size=args.rand, replace=False)
            _tee(out, 'Selected random subset of %d trajectories.' % args.rand)

        # If slice values provided, select the files needed. Again before
        # reading files so speed up the process
        if args.slice is not None:
            first = args.slice[0]
            last = args.slice[1]
            _tee(out, ' First trajectory read: %s' % first )
            _tee(out, ' Last trajectory read: %s' % last )
            _tee(out, '')
            if 'A' in statesProvided:
                filesAB = filesAB[first:last]
            if 'B' in statesProvided:
                filesBA = filesBA[first:last]

        # If index values provided, select the files needed
        if args.index is not None:
            # Avoid index out of range error if "wrong" indices are provided
            if 'A' in statesProvided:
                filesAB = [filesAB[i] for i in args.index if i < len(filesAB)]
            if 'B' in statesProvided:
                filesBA = [filesBA[i] for i in args.index if i < len(filesBA)]
            # ...but warn if this happens
            if any(i > (len(filesAB) - 1) for i in args.index):
                print('\nWARNING: index out of range for some of your chosen '
                      '\nindices for the forward work values. This means you are'
                      '\ntrying to select input files that are not present.')
            if any(i > (len(filesBA) - 1) for i in args.index):
                print('\nWARNING: index out of range for some of your chosen'
                      '\nindices for the reverse work values. This means you are'
                      '\ntrying to select input files that are not present.')

        # when skipping start count from end: in this way the last frame is
        # always included, and what can change is the first one
        if 'A' in statesProvided:
            filesAB = list(reversed(filesAB[::-skip]))
        if 'B' in statesProvided:
            filesBA = list(reversed(filesBA[::-skip]))

        # --------------------
        # Now read in the data
        # --------------------
        print(' ========================================================')
        print('                   PROCESSING THE DATA')
        print(' ========================================================')
        res_ab = []
        res_ba = []
        if 'A' in statesProvided:
            print('  Forward Data')
            res_ab = parse_dgdl_files(filesAB, lambda0=0,invert_values=False)
            _dump_integ_file(args.oA, filesAB, res_ab)
        if 'B' in statesProvided:
            print('  Reverse Data')
            res_ba = parse_dgdl_files(filesBA, lambda0=1,invert_values=reverseB)
            _dump_integ_file(args.oB, filesBA, res_ba)

    # If work values are given as input instead, read those
    elif args.iA is not None or args.iB is not None:
        res_ab = []
        res_ba = []
        if 'A' in statesProvided:
            print('  Reading integrated values (A->B) from', args.iA)
            res_ab.extend(_data_from_file(args.iA))
        if 'B' in statesProvided:
            print('  Reading integrated values (B->A) from', args.iB)
            res_ba.extend(_data_from_file(args.iB))
        # If slice values provided, select the files needed.
        if args.slice is not None:
            first = args.slice[0]
            last = args.slice[1]
            _tee(out, '  Keeping work values from %i to %i.' % (first, last))
            _tee(out, '')
            if 'A' in statesProvided:
                res_ab = res_ab[first:last]
            if 'B' in statesProvided:
                res_ba = res_ba[first:last]
    else:
        exit('\nERROR: you need to provide either none of both sets of '
             'integrated work values.')

    # If asked to only do the integration of dhdl.xvg, exit
    if integ_only:
        print('\n    Integration done. Skipping analysis.')
        print('\n    ......done........\n')
        sys.exit(0)

    # ==============
    # Begin Analysis
    # ==============
    _tee(out, ' ========================================================')
    _tee(out, '                       ANALYSIS')
    _tee(out, ' ========================================================')
    if 'A' in statesProvided:
        _tee(out, '  Number of forward (0->1) trajectories: %d' % len(res_ab))
    if 'B' in statesProvided:
        _tee(out, '  Number of reverse (1->0) trajectories: %d' % len(res_ba))
    _tee(out, '  Temperature : %.2f K' % T)

    # ============================
    # Crooks Gaussian Intersection
    # ============================
    if 'cgi' in methods and 'AB' in statesProvided:
        _tee(out, '\n --------------------------------------------------------')
        _tee(out, '             Crooks Gaussian Intersection     ')
        _tee(out, ' --------------------------------------------------------')

        print('  Calculating Intersection...')
        cgi = Crooks(wf=res_ab, wr=res_ba, nboots=nboots, nblocks=nblocks)
        if args.pickle is True:
            pickle.dump(cgi, open("cgi_results.pkl", "wb"))

        _tee(out, '  CGI: Forward Gauss mean = {m:8.{p}f} {u} '
                  'std = {s:8.{p}f} {u}'.format(m=cgi.mf*unit_fact,
                                                s=cgi.devf*unit_fact,
                                                p=prec, u=units))
        _tee(out, '  CGI: Reverse Gauss mean = {m:8.{p}f} {u} '
                  'std = {s:8.{p}f} {u}'.format(m=cgi.mr*unit_fact,
                                                s=cgi.devr*unit_fact,
                                                p=prec, u=units))

        if cgi.inters_bool is False:
            _tee(out, '\n  Gaussians too close for intersection calculation')
            _tee(out, '   --> Taking difference of mean values')

        _tee(out, '  CGI: dG = {dg:8.{p}f} {u}'.format(dg=cgi.dg*unit_fact,
                                                       p=prec, u=units))
        _tee(out, '  CGI: Std Err (bootstrap:parametric) = {e:8.{p}f} {u}'.format(e=cgi.err_boot1*unit_fact,
                                                                                  p=prec, u=units))

        if nboots > 0:
            _tee(out, '  CGI: Std Err (bootstrap) = {e:8.{p}f} {u}'.format(e=cgi.err_boot2*unit_fact,
                                                                           p=prec, u=units))

        if nblocks > 1:
            _tee(out, '  CGI: Std Err (blocks) = {e:8.{p}f} {u}'.format(e=cgi.err_blocks*unit_fact,
                                                                        p=prec, u=units))

    # --------------
    # Normality test
    # --------------
    if do_ks_test and 'AB' in statesProvided:
        print('\n  Running KS-test...')
        q0, lam00, check0, bOk0 = ks_norm_test(res_ab)
        q1, lam01, check1, bOk1 = ks_norm_test(res_ba)

        _tee(out, '    Forward: gaussian quality = %3.2f' % q0)
        if bOk0:
            _tee(out, '             ---> KS-Test Ok')
        else:
            _tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f,'
                      ' lambda0 = %4.2f' % (q0, check0))
        _tee(out, '    Reverse: gaussian quality = %3.2f' % q1)
        if bOk1:
            _tee(out, '             ---> KS-Test Ok')
        else:
            _tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f,'
                      ' lambda0 = %4.2f' % (q1, check1))

    # ========================
    # Bennett Acceptance Ratio
    # ========================
    if 'bar' in methods and 'AB' in statesProvided:
        _tee(out, '\n --------------------------------------------------------')
        _tee(out, '             Bennett Acceptance Ratio     ')
        _tee(out, ' --------------------------------------------------------')

        print('  Running Nelder-Mead Simplex algorithm... ')

        bar = BAR(res_ab, res_ba, T=T, nboots=nboots, nblocks=nblocks)
        if args.pickle:
            pickle.dump(bar, open("bar_results.pkl", "wb"))

        _tee(out, '  BAR: dG = {dg:8.{p}f} {u}'.format(dg=bar.dg*unit_fact, p=prec, u=units))
        _tee(out, '  BAR: Std Err (analytical) = {e:8.{p}f} {u}'.format(e=bar.err*unit_fact, p=prec, u=units))

        if nboots > 0:
            _tee(out, '  BAR: Std Err (bootstrap)  = {e:8.{p}f} {u}'.format(e=bar.err_boot*unit_fact, p=prec, u=units))
        if nblocks > 1:
            _tee(out, '  BAR: Std Err (blocks)  = {e:8.{p}f} {u}'.format(e=bar.err_blocks*unit_fact, p=prec, u=units))

        _tee(out, '  BAR: Conv = %8.2f' % bar.conv)

        if nboots > 0:
            _tee(out, '  BAR: Conv Std Err (bootstrap) = %8.2f' % bar.conv_err_boot)

    # =========
    # Jarzynski
    # =========
    if 'jarz' in methods or 'AB' not in statesProvided:
        _tee(out, '\n --------------------------------------------------------')
        _tee(out, '             Jarzynski estimator     ')
        _tee(out, ' --------------------------------------------------------')

        jarz = Jarz(wf=res_ab, wr=res_ba, T=T, nboots=nboots, nblocks=nblocks, statesProvided=statesProvided)
        if args.pickle:
            pickle.dump(jarz, open("jarz_results.pkl", "wb"))

        if 'A' in statesProvided:
            _tee(out, '  JARZ: dG Forward = {dg:8.{p}f} {u}'.format(dg=jarz.dg_for*unit_fact,
                                                                p=prec, u=units))
        if 'B' in statesProvided:
            _tee(out, '  JARZ: dG Reverse = {dg:8.{p}f} {u}'.format(dg=jarz.dg_rev*unit_fact,
                                                                p=prec, u=units))
        if 'AB' in statesProvided:
            _tee(out, '  JARZ: dG Mean    = {dg:8.{p}f} {u}'.format(dg=jarz.dg_mean*unit_fact,
                                                                p=prec, u=units))
        if nboots > 0:
            if 'A' in statesProvided:
                _tee(out, '  JARZ: Std Err Forward (bootstrap) = {e:8.{p}f} {u}'.format(e=jarz.err_boot_for*unit_fact,
                                                                                    p=prec, u=units))
            if 'B' in statesProvided:
                _tee(out, '  JARZ: Std Err Reverse (bootstrap) = {e:8.{p}f} {u}'.format(e=jarz.err_boot_rev*unit_fact,
                                                                                    p=prec, u=units))

        if nblocks > 1:
            if 'A' in statesProvided:
                _tee(out, '  JARZ: Std Err Forward (blocks) = {e:8.{p}f} {u}'.format(e=jarz.err_blocks_for*unit_fact,
                                                                                 p=prec, u=units))
            if 'B' in statesProvided:
                _tee(out, '  JARZ: Std Err Reverse (blocks) = {e:8.{p}f} {u}'.format(e=jarz.err_blocks_rev*unit_fact,
                                                                                 p=prec, u=units))


        # -------------------------------------
        # Jarzynski with Gaussian approximation
        # -------------------------------------
        print('Running Jarzynski Gaussian approximation analysis...')
        jarzGauss = JarzGauss(wf=res_ab, wr=res_ba, T=T, nboots=nboots, nblocks=nblocks, statesProvided=statesProvided)
        if args.pickle:
            pickle.dump(jarzGauss, open("jarz_gauss_results.pkl", "wb"))

        if 'A' in statesProvided:
            _tee(out, '  JARZ_Gauss: dG Forward = {dg:8.{p}f} {u}'.format(dg=jarzGauss.dg_for*unit_fact,
                                                                    p=prec, u=units))
        if 'B' in statesProvided:
            _tee(out, '  JARZ_Gauss: dG Reverse = {dg:8.{p}f} {u}'.format(dg=jarzGauss.dg_rev*unit_fact,
                                                                    p=prec, u=units))
        if 'AB' in statesProvided:
            _tee(out, '  JARZ_Gauss: dG Mean    = {dg:8.{p}f} {u}'.format(dg=(jarzGauss.dg_for+jarzGauss.dg_rev)/2.0*unit_fact,
                                                                    p=prec, u=units))
        if 'A' in statesProvided:
            _tee(out, '  JARZ_Gauss: Std Err (analytical) Forward = {dg:8.{p}f} {u}'.format(dg=jarzGauss.err_for*unit_fact,
                                                                    p=prec, u=units))
        if 'B' in statesProvided:
            _tee(out, '  JARZ_Gauss: Std Err (analytical) Reverse = {dg:8.{p}f} {u}'.format(dg=jarzGauss.err_rev*unit_fact,
                                                                    p=prec, u=units))
        if nboots > 0:
            if 'A' in statesProvided:
                _tee(out, '  JARZ_Gauss: Std Err Forward (bootstrap) = {e:8.{p}f} {u}'.format(e=jarzGauss.err_boot_for*unit_fact, p=prec, u=units))
            if 'B' in statesProvided:
                _tee(out, '  JARZ_Gauss: Std Err Reverse (bootstrap) = {e:8.{p}f} {u}'.format(e=jarzGauss.err_boot_rev*unit_fact,p=prec, u=units))

        if nblocks > 1:
            if 'A' in statesProvided:
                _tee(out, '  JARZ_Gauss: Std Err Forward (blocks) = {e:8.{p}f} {u}'.format(e=jarzGauss.err_blocks_for*unit_fact,p=prec, u=units))
            if 'B' in statesProvided:
                _tee(out, '  JARZ_Gauss: Std Err Reverse (blocks) = {e:8.{p}f} {u}'.format(e=jarzGauss.err_blocks_rev*unit_fact,p=prec, u=units))

    _tee(out, ' ========================================================')

    # -----------------------
    # plot work distributions
    # -----------------------
    if args.wplot.lower() != 'none':
        print('\n   Plotting histograms......')
        # hierarchy of estimators: BAR > Crooks > Jarz
        if 'bar' in locals():
            show_dg = bar.dg * unit_fact
            # hierarchy of error estimates : blocks > boots > analytical
            if hasattr(bar, 'err_blocks'):
                show_err = bar.err_blocks * unit_fact
            elif hasattr(bar, 'err_boot') and not hasattr(bar, 'err_blocks'):
                show_err = bar.err_boot * unit_fact
            else:
                show_err = bar.err * unit_fact
            # plot
            plot_work_dist(fname=args.wplot, wf=res_ab, wr=res_ba, dG=show_dg,
                           dGerr=show_err, nbins=args.nbins, dpi=args.dpi,
                           units=units)
        elif 'bar' not in locals() and 'cgi' in locals():
            show_dg = cgi.dg * unit_fact
            # hierarchy of error estimates : blocks > boots
            if hasattr(cgi, 'err_blocks'):
                show_err = cgi.err_blocks * unit_fact
            elif hasattr(cgi, 'err_boot2') and not hasattr(cgi, 'err_blocks'):
                show_err = cgi.err_boot2 * unit_fact
            else:
                show_err = None
            # plot
            plot_work_dist(fname=args.wplot, wf=res_ab, wr=res_ba, dG=show_dg,
                           dGerr=show_err, nbins=args.nbins, dpi=args.dpi,
                           units=units)
        elif 'bar' not in locals() and 'cgi' not in locals() and 'jarz' in locals():
            # for the moment, show values only under specific circumstances
            if hasattr(jarz, 'dg_mean'):
                show_dg = jarz.dg_mean * unit_fact
            elif 'A' in statesProvided:
                show_dg = jarz.dg_for
            elif 'B' in statesProvided:
                show_dg = jarz.dg_rev
            else:
                show_dg = None
            show_err = None
            # plot
            plot_work_dist(fname=args.wplot, wf=res_ab, wr=res_ba, dG=show_dg,
                           dGerr=show_err, nbins=args.nbins, dpi=args.dpi,
                           units=units,statesProvided=statesProvided)

    print('\n   ......done...........\n')

    if args.pickle:
        print('   NOTE: units of results in pickled files are as in the\n'
              '   provided dhdl.xvg or integ.dat files. These are typically\n'
              '   in kJ/mol when using dhdl.xvg files from Gromacs.\n')
    # execution time
    etime = time.time()
    h, m, s = time_stats(etime-stime)
    print("   Execution time = %02d:%02d:%02d\n" % (h, m, s))


def entry_point():
    args = parse_options()
    main(args)


if __name__ == '__main__':
    entry_point()
