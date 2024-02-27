
import numpy as np
import sys
from scipy.optimize import fmin
from scipy.special import erf
import scipy.stats
from copy import deepcopy

# Constants
kb = 0.00831447215   # kJ/(K*mol)


class Jarz(object):
    '''Jarzynski estimator.

    Description...

    Parameters
    ----------
    wf : array_like
        array of forward work values.
    wr : array_like
        array of reverse work values.
    T : float or int
    nboots : int
        number of bootstrap samples to use for error estimation.

    Examples
    --------


    Attributes
    ----------

    '''

    def __init__(self, wf, wr, T, nboots=0, nblocks=1, statesProvided='AB'):
        if 'A' in statesProvided:
            self.wf = np.array(wf)
        if 'B' in statesProvided:
            self.wr = np.array(wr)
        self.T = float(T)
        self.nboots = nboots
        self.nblocks = nblocks

        # Calculate all Jarz properties available
        if 'A' in statesProvided:
            self.dg_for = self.calc_dg(w=self.wf, c=1.0, T=self.T)
        if 'B' in statesProvided:
            self.dg_rev = -1.0 * self.calc_dg(w=self.wr, c=-1.0, T=self.T)
        if 'AB' in statesProvided:
            self.dg_mean = (self.dg_for + self.dg_rev) * 0.5

        if nboots > 0:
            if 'A' in statesProvided:
                self.err_boot_for = self.calc_err_boot(w=self.wf, T=self.T,
                                                   c=1.0, nboots=nboots)
            if 'B' in statesProvided:
                self.err_boot_rev = self.calc_err_boot(w=self.wr, T=self.T,
                                                   c=-1.0, nboots=nboots)

        if nblocks > 1:
            if 'A' in statesProvided:
                self.err_blocks_for = self.calc_err_blocks(w=self.wf, c=1.0,
                                                       T=self.T,
                                                       nblocks=nblocks)
            if 'B' in statesProvided:
                self.err_blocks_rev = self.calc_err_blocks(w=self.wr, c=-1.0,
                                                       T=self.T,
                                                       nblocks=nblocks)

    @staticmethod
    def calc_dg(w, T, c):
        '''to be filled
        '''
        beta = 1./(kb*T)
        n = float(len(w))
        mexp = 0.0
        m = 0.0
        m2 = 0.0
        for i in w:
            mexp = mexp + np.exp(-beta*c*i)
            m = m + c*i
            m2 = m2 + i*i
        mexp = mexp/n
        m = m/n
        m2 = m2/n
        var = (m2-m*m)*(n/(n-1))
        # Jarzynski estimator
        dg = -kb*T*np.log(mexp)
        # Fluctuation-Dissipation estimator
        # FIXME: unused atm, remove or return?
        dg2 = m - beta*var/2.0

        return dg

    @staticmethod
    def calc_err_boot(w, T, c, nboots):
        '''Calculates the standard error via bootstrap. The work values are
        resampled randomly with replacement multiple (nboots) times,
        and the Jarzinski free energy recalculated for each bootstrap samples.
        The standard error of the estimate is returned as the standard d
        eviation of the bootstrapped free energies.

        Parameters
        ----------
        w : array_like
            work values.
        T : float
            temperature.
        c : [0,1]
            ???
        nboots: int
            number of bootstrap samples to use for the error estimate.

        Returns
        -------
        err : float
            standard error of the mean.
        '''
        dg_boots = []
        n = len(w)
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Std Err): iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()

            boot = np.random.choice(w, size=n, replace=True)
            dg_boot = -1.0 * Jarz.calc_dg(boot, T, c)
            dg_boots.append(dg_boot)
        sys.stdout.write('\n')
        err = np.std(dg_boots)
        return err

    @staticmethod
    def calc_err_blocks(w, T, c, nblocks):
        '''Calculates the standard error based on a number of blocks the
        work values are divided into. It is useful when you run independent
        equilibrium simulations, so that you can then use their respective
        work values to compute the standard error based on the repeats.

        Parameters
        ----------
        w : array_like
            array of work values.
        T : float
            temperature.
        c : [0,1]
            ???
        nblocks: int
            number of blocks to divide the data into. This can be for
            instance the number of independent equilibrium simulations
            you ran.
        '''

        dg_blocks = []
        # loosely split the arrays
        w_split = np.array_split(w, nblocks)

        # calculate all dg
        for w_block in w_split:
            dg_block = -1.0 * Jarz.calc_dg(w_block, T, c)
            dg_blocks.append(dg_block)

        # get std err
        err_blocks = scipy.stats.sem(dg_blocks, ddof=1)

        return err_blocks


class JarzGauss:
    '''Jarzynski estimator using a Gaussian approximation. [6]_
    Both the forward and reverse estimates.
    The standard error is calculated using the analytical expression derived
    by Hummer (2001). [6]_ When ``nboots``>0, the error is also calculated
    by bootstrap. When  ``nblocks``>1, it is calculated by
    separating the work values into blocks.
    Parameters
    ----------
    wf : array_like
        array of forward work values.
    wr : array_like
        array of reverse work values.
    T : float, optional
        temperature in Kelvin. Default is 298.15 K.
    nboots : int, optional
        how many bootstrap samples to draw for estimating the standard error.
        Default is zero (do not estimate the error).
    nblocks : int, optional
        how many blocks to divide the input work values into for the estimation
        of the standard error. Default is one (do not estimate the error).
    statesProvided: str, optional
        two directions or one
    Examples
    --------
    >>> estimate = JarzGauss(wf, wr, T=300, nboots=1000, nblocks=10)
    >>> dg_forward = estimate.dg_for
    >>> dg_reverse = estimate.dg_rev
    >>> dg_forward_err1 = estimate.err_for
    >>> dg_forward_err2 = estimate.err_boot_for
    >>> dg_forward_err3 = estimate.err_blocks_for
    Attributes
    ----------
    dg_for : float
        the forward free energy estimate.
    dg_rev : float
        the reverse free energy estimate.
    err_for : float
        analytical estimate of the standard error of the forward free energy
        estimate.
    err_rev : float
        analytical estimate of the standard error of the reverse free energy
        estimate.
    err_boot_for : float
        standard error of the forward free energy estimate calculated
        via bootstrap.
    err_boot_rev : float
        standard error of the reverse free energy estimate calculated
        via bootstrap.
    err_blocks_for : float
        standard error of the forward free energy estimate calculated by
        separating the input work values into groups/blocks.
    err_blocks_rev : float
        standard error of the reverse free energy estimate calculated by
        separating the input work values into groups/blocks.
    '''

    def __init__(self, wf, wr, T=298.15, nboots=0, nblocks=1, statesProvided='AB'):
        if 'A' in statesProvided:
            self.wf = np.array(wf)
        if 'B' in statesProvided:
            self.wr = np.array(wr)
        self.T = float(T)
        self.nboots = nboots
        self.nblocks = nblocks

        # Calculate all Jarz properties available
        if 'A' in statesProvided:
            self.dg_for = self.calc_dg(w=self.wf, T=self.T, bReverse=False)
            self.err_for = self.calc_err(w=self.wf, T=self.T, bReverse=False)
        if 'B' in statesProvided:
            self.dg_rev = self.calc_dg(w=self.wr, T=self.T, bReverse=True)
            self.err_rev = self.calc_err(w=self.wr, T=self.T, bReverse=True)

        if nboots > 0:
            if 'A' in statesProvided:
                self.err_boot_for = self.calc_err_boot(w=self.wf, T=self.T,
                                                       nboots=self.nboots,
                                                       bReverse=False)
            if 'B' in statesProvided:
                self.err_boot_rev = self.calc_err_boot(w=self.wr, T=self.T,
                                                       nboots=self.nboots,
                                                       bReverse=True)

        if nblocks > 1:
            if 'A' in statesProvided:
                self.err_blocks_for = self.calc_err_blocks(w=self.wf,
                                                           T=self.T,
                                                           nblocks=self.nblocks,
                                                           bReverse=False)
            if 'B' in statesProvided:
                self.err_blocks_rev = self.calc_err_blocks(w=self.wr,
                                                           T=self.T,
                                                           nblocks=self.nblocks,
                                                           bReverse=True)

    @staticmethod
    def calc_dg(w, T, bReverse=False):
        '''Calculates the free energy difference using the Jarzynski estimator
        with a Gaussian approximation. [6]_
        Parameters
        ----------
        w : array_like
            array of work values.
        T : float
            temperature in Kelvin.
        bReverse : bool, optional
            whether the work values provided are for the reverse transition.
            Default if False. If they are for the reverse transition, set it to
            True.
        Returns
        -------
        dg : float
            estimate of the free energy difference.
        '''
        beta = 1./(kb*T)
        if bReverse is False:
            c = 1.0
        elif bReverse is True:
            c = -1.0

        dg = np.mean(c*w) - (beta * np.var(c*w, ddof=1)) * 0.5
        return c * dg

    @staticmethod
    def calc_err(w, T, bReverse=False):
        '''Calculates the standard error via an analytic expression.
        The expression is derived by Hummer, 2001, JChemPhys. [6]_
        Parameters
        ----------
        w : array_like
            work values.
        T : float
            temperature.
        bReverse : bool, optional
            whether the work values provided are for the reverse transition.
            Default if False. If they are for the reverse transition, set it to
            True.
        Returns
        -------
        err : float
            standard deviation of the estimator.
        '''

        beta = 1./(kb*T)
        w_var = np.var(w, ddof=1)
        n = float(len(w))
        dg_var = w_var/n + np.power(beta*w_var, 2) / (2.0 * (n-1.0))
        dg_stderr = np.sqrt(dg_var)

        return dg_stderr

    @staticmethod
    def calc_err_boot(w, T, nboots, bReverse=False):
        '''Calculates the standard error via bootstrap. The work values are
        resampled randomly with replacement multiple (nboots) times,
        and the Gaussian approximation for Jarzinski free energy
        is recalculated for each bootstrap samples.
        The standard error of the estimate is returned as the standard d
        eviation of the bootstrapped free energies.
        Parameters
        ----------
        w : array_like
            work values.
        T : float
            temperature.
        nboots: int
            number of bootstrap samples to use for the error estimate.
        bReverse : bool, optional
            whether the work values provided are for the reverse transition.
            Default if False. If they are for the reverse transition, set it to
            True.
        Returns
        -------
        err : float
            standard error of the mean.
        '''
        dg_boots = []
        n = len(w)
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Std Err): iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()

            boot = np.random.choice(w, size=n, replace=True)
            dg_boot = JarzGauss.calc_dg(boot, T, bReverse)
            dg_boots.append(dg_boot)
        sys.stdout.write('\n')
        err = np.std(dg_boots)
        return err

    @staticmethod
    def calc_err_blocks(w, T, nblocks, bReverse=False):
        '''Calculates the standard error based on a number of blocks the
        work values are divided into. It is useful when you run independent
        equilibrium simulations, so that you can then use their respective
        work values to compute the standard error based on the repeats.
        Parameters
        ----------
        w : array_like
            array of work values.
        T : float
            temperature.
        nblocks: int
            number of blocks to divide the data into. This can be for
            instance the number of independent equilibrium simulations
            you ran.
        bReverse : bool, optional
            whether the work values provided are for the reverse transition.
            Default if False. If they are for the reverse transition, set it to
            True.
        Returns
        -------
        sderr : float
            the standard error of the estimate.
        '''

        dg_blocks = []
        # loosely split the arrays
        w_split = np.array_split(w, nblocks)

        # calculate all dg
        for w_block in w_split:
            dg_block = JarzGauss.calc_dg(w_block, T, bReverse)
            dg_blocks.append(dg_block)

        # get std err
        err_blocks = scipy.stats.sem(dg_blocks, ddof=1)

        return err_blocks


class Crooks(object):
    '''Crooks Gaussian Intersection (CGI) estimator. The forward and reverse work
    values are fitted to Gaussian functions and their intersection is taken
    as the free energy estimate. In some cases, when the two Gaussians are very
    close to each other, the intersection cannot be taken and the average of
    the two Gaussian means is taken as the free energy estimate insted. The
    standard error is by default calculated via bootstrap using 1000 samples.

    Parameters
    ----------
    wf : array_like
        array of forward work values.
    wr : array_like
        array of reverse work values.

    Examples
    --------
    >>> Crooks(wf, wr)
    >>> cgi_dg = Crooks.dg
    >>> cgi_err = Crooks.err_boot

    Attributes
    ----------
    dg : float
        the free energy estimate.
    err_boot : float
        standard error of the free energy estimate calculated via bootstrap.
    inters_bool : bool
        whether the interection could be taken. If False, the free energy
        estimate is the average of the two Gaussian means.
    Af : float
        height of the forward Gaussian.
    mf : float
        mean of the forward Gaussian.
    devf : float
        standard deviation of the forward Gaussian.
    Ar : float
        height of the reverse Gaussian.
    mr : float
        mean of the reverse Gaussian.
    devr : float
        standard deviation of the reverse Gaussian.
    '''

    def __init__(self, wf, wr, nboots=0, nblocks=1):

        # inputs
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.nboots = nboots
        self.nblocks = nblocks
        # params of the gaussians
        self.mf, self.devf, self.Af = data2gauss(wf)
        self.mr, self.devr, self.Ar = data2gauss(wr)

        # Calculate Crooks properties
        self.dg, self.inters_bool = self.calc_dg(wf=self.wf, wr=self.wr)

        self.err_boot1 = self.calc_err_boot1(m1=self.mf, s1=self.devf,
                                             n1=len(wf), m2=self.mr,
                                             s2=self.devr, n2=len(wr),
                                             nboots=1000)
        if nboots > 0:
            self.err_boot2 = self.calc_err_boot2(wf=self.wf, wr=self.wr,
                                                 nboots=nboots)

        if nblocks > 1:
            self.err_blocks = self.calc_err_blocks(self.wf, self.wr, nblocks)

    @staticmethod
    def calc_dg(wf, wr):
        '''Calculates the free energy difference using the Crooks Gaussian
        Intersection method. It finds the intersection of two Gaussian
        functions. If the intersection cannot be computed, the average of
        the two Gaussian locations is returned.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.

        Returns
        -------
        float
            location of the intersection.
        bool
            whether the intersection could be calculated. If the intersection
            was calculated as expected a True value is returned.
            If the Gaussians are too close to each other, the intersection
            cannot be calculated and a False value is returned; in this case,
            the first float value retured is the average of the Gaussian means.
        '''

        m1, s1, A1 = data2gauss(wf)
        m2, s2, A2 = data2gauss(wr)

        p1 = m1/s1**2-m2/s2**2
        p2 = np.sqrt(1/(s1**2*s2**2)*(m1-m2)**2+2*(1/s1**2-1/s2**2)*np.log(s2/s1))
        p3 = 1/s1**2-1/s2**2
        x1 = (p1+p2)/p3
        x2 = (p1-p2)/p3

        # determine which solution to take
        if x1 > m1 and x1 < m2 or x1 > m2 and x1 < m1:
            dg = x1
            return dg, True
        elif x2 > m1 and x2 < m2 or x2 > m2 and x2 < m1:
            dg = x2
            return dg, True
        else:
            # we do not take the intersection but the average of the means
            dg = (m1 + m2) * 0.5
            return dg, False

    # Possible change of behaviour compared to the original script:
    # here it is not determined in advanced whether to take the intersection
    # or the mean, but for each bootstrap sample if the intersecion cannot
    # be taken, then the mean is used automatically.
    @staticmethod
    def calc_err_boot1(m1, s1, n1, m2, s2, n2, nboots=1000):
        '''Calculates the standard error of the Crooks Gaussian Intersection
        via parametric bootstrap. Given the parameters of the forward and
        reverse Gaussian distributions, multiple (nboots) bootstrap samples
        are built by random sampling from these two Gaussian distributions.
        The CGI free energy is then calculated for each bootstrap sample
        (forward and reverse Gaussians). The standard error of the estimate
        is returned as the standard deviation of the bootstrapped free
        energies.

        Parameters
        ----------
        m1 : float
            mean of the forward Gaussian.
        s1 : float
            standard deviation of the forward Gaussian.
        n1 : int
            number of work values to which the first Gaussian was fit.
        m2 : float
            mean of the reverse Gaussian.
        s2 : float
            standard deviation of the reverse Gaussian.
        n2 : int
            number of work values to which the second Gaussian was fit.
        nboots: int
            number of bootstrap samples to use for the error estimate.
            Parametric bootstrap is used where work values are resampled from
            two Gaussians.

        Returns
        -------
        float
            standard error of the mean.
        '''

        dg_boots = []
        for k in range(nboots):
            bootA = np.random.normal(loc=m1, scale=s1, size=n1)
            bootB = np.random.normal(loc=m2, scale=s2, size=n2)

            dg_boot, _ = Crooks.calc_dg(bootA, bootB)
            dg_boots.append(dg_boot)
        err = np.std(dg_boots)
        return err

    @staticmethod
    def calc_err_boot2(wf, wr, nboots):
        '''Calculates the standard error of the Crooks Gaussian Intersection
        via non-parametric bootstrap. The work values are resampled randomly
        with replacement multiple (nboots) times, and the CGI free energy
        recalculated for each bootstrap samples. The standard error of
        the estimate is returned as the standard deviation of the bootstrapped
        free energies.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        nboots: int
            number of bootstrap samples to use for the error estimate.

        Returns
        -------
        err : float
            standard error of the mean.
        '''

        nf = len(wf)
        nr = len(wr)

        dg_boots = []
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Std Err): iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()

            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)

            dg_boot, _ = Crooks.calc_dg(bootA, bootB)
            dg_boots.append(dg_boot)

        sys.stdout.write('\n')
        err = np.std(dg_boots)
        return err

    @staticmethod
    def calc_err_blocks(wf, wr, nblocks):
        '''Calculates the standard error based on a number of blocks the
        work values are divided into. It is useful when you run independent
        equilibrium simulations, so that you can then use their respective
        work values to compute the standard error based on the repeats.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        nblocks: int
            number of blocks to divide the data into. This can be for
            instance the number of independent equilibrium simulations
            you ran.
        '''

        dg_blocks = []
        # loosely split the arrays
        wf_split = np.array_split(wf, nblocks)
        wr_split = np.array_split(wr, nblocks)

        # calculate all dg
        for wf_block, wr_block in zip(wf_split, wr_split):
            dg_block, _ = Crooks.calc_dg(wf_block, wr_block)
            dg_blocks.append(dg_block)

        # get std err
        err_blocks = scipy.stats.sem(dg_blocks, ddof=1)

        return err_blocks


class BAR(object):
    '''Bennett acceptance ratio (BAR).

    Description...

    Parameters
    ----------

    Examples
    --------
    '''

    def __init__(self, wf, wr, T, nboots=0, nblocks=1):
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.T = float(T)
        self.nboots = nboots
        self.nblocks = nblocks

        self.nf = len(wf)
        self.nr = len(wr)
        self.beta = 1./(kb*self.T)
        self.M = kb * self.T * np.log(float(self.nf) / float(self.nr))

        # Calculate all BAR properties available
        self.dg = self.calc_dg(self.wf, self.wr, self.T)
        self.err = self.calc_err(self.dg, self.wf, self.wr, self.T)
        if nboots > 0:
            self.err_boot = self.calc_err_boot(self.wf, self.wr, nboots,
                                               self.T)
        self.conv = self.calc_conv(self.dg, self.wf, self.wr, self.T)
        if nboots > 0:
            self.conv_err_boot = self.calc_conv_err_boot(self.dg, self.wf,
                                                         self.wr, nboots,
                                                         self.T)
        if nblocks > 1:
            self.err_blocks = self.calc_err_blocks(self.wf, self.wr, nblocks,
                                                   self.T)

    @staticmethod
    def calc_dg(wf, wr, T):
        '''Estimates and returns the free energy difference.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature

        Returns
        ----------
        dg : float
            the BAR free energy estimate.
        '''

        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1./(kb*T)
        M = kb * T * np.log(nf/nr)

        def func(x, wf, wr):
            sf = 0
            for v in wf:
                sf += 1./(1+np.exp(beta*(M+v-x)))

            sr = 0
            for v in wr:
                sr += 1./(1+np.exp(-beta*(M+v-x)))

            r = sf-sr
            return r**2

        avA = np.average(wf)
        avB = np.average(wr)
        x0 = (avA+avB)/2.
        dg = fmin(func, x0=x0, args=(wf, wr), disp=0)

        return float(dg)

    @staticmethod
    def calc_err(dg, wf, wr, T):
        '''Calculates the analytical error estimate.

        Parameters
        ----------
        dg : float
            the BAR free energy estimate
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        '''

        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1./(kb*T)
        M = kb * T * np.log(nf/nr)

        err = 0
        for v in wf:
            err += 1./(2+2*np.cosh(beta * (M+v-dg)))
        for v in wr:
            err += 1./(2+2*np.cosh(beta * (M+v-dg)))
        N = nf + nr
        err /= float(N)
        tot = 1/(beta**2*N)*(1./err-(N/nf + N/nr))

        err = float(np.sqrt(tot))
        return err

    @staticmethod
    def calc_err_boot(wf, wr, nboots, T):
        '''Calculates the error by bootstrapping.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        nboots: int
            number of bootstrap samples.

        '''

        nf = len(wf)
        nr = len(wr)
        dg_boots = []
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Std Err): iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()

            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)
            dg_boot = BAR.calc_dg(bootA, bootB, T)
            dg_boots.append(dg_boot)

        sys.stdout.write('\n')
        err_boot = np.std(dg_boots)

        return err_boot

    @staticmethod
    def calc_err_blocks(wf, wr, nblocks, T):
        '''Calculates the standard error based on a number of blocks the
        work values are divided into. It is useful when you run independent
        equilibrium simulations, so that you can then use their respective
        work values to compute the standard error based on the repeats.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        nblocks: int
            number of blocks to divide the data into. This can be for
            instance the number of independent equilibrium simulations
            you ran.
        '''

        dg_blocks = []
        # loosely split the arrays
        wf_split = np.array_split(wf, nblocks)
        wr_split = np.array_split(wr, nblocks)

        # calculate all dg
        for wf_block, wr_block in zip(wf_split, wr_split):
            dg_block = BAR.calc_dg(wf_block, wr_block, T)
            dg_blocks.append(dg_block)

        # get std err
        err_blocks = scipy.stats.sem(dg_blocks, ddof=1)

        return err_blocks

    @staticmethod
    def calc_conv(dg, wf, wr, T):
        '''Evaluates BAR convergence as described in Hahn & Then, Phys Rev E
        (2010), 81, 041117. Returns a value between -1 and 1: the closer this
        value to zero the better the BAR convergence.

        Parameters
        ----------
        dg : float
            the BAR free energy estimate
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature

        '''

        wf = np.array(wf)
        wr = np.array(wr)

        beta = 1./(kb*T)
        nf = len(wf)
        nr = len(wr)
        N = float(nf + nr)

        ratio_alpha = float(nf)/N
        ratio_beta = float(nr)/N
        bf = 1.0/(ratio_beta + ratio_alpha * np.exp(beta*(wf-dg)))
        tf = 1.0/(ratio_alpha + ratio_beta * np.exp(beta*(-wr+dg)))
        Ua = (np.mean(tf) + np.mean(bf))/2.0
        Ua2 = (ratio_alpha * np.mean(np.power(tf, 2)) +
               ratio_beta * np.mean(np.power(bf, 2)))
        conv = (Ua-Ua2)/Ua
        return conv

    @staticmethod
    def calc_conv_err_boot(dg, wf, wr, nboots, T):
        nf = len(wf)
        nr = len(wr)
        conv_boots = []
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Conv): '
                             'iteration %s/%s' % (k+1, nboots))
            sys.stdout.flush()

            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)
            conv_boot = BAR.calc_conv(dg, bootA, bootB, T)
            conv_boots.append(conv_boot)

        sys.stdout.write('\n')
        err = np.std(conv_boots)
        return err


# ==============================================================================
#                               FUNCTIONS
# ==============================================================================
def ks_norm_test(data, alpha=0.05, refks=None):
    '''Performs a Kolmogorov-Smirnov test of normality.

    Parameters
    ----------
    data : array_like
        a one-dimensional array of values. This is the distribution tested
        for normality.
    alpha : float
        significance level of the statistics. Default if 0.05.
    refks : ???
        ???

    Returns
    -------
    Q : float
    lam0 : float
    check : float
    bOk : bool
    '''

    def ksref():
        f = 1
        potent = 10000
        lamb = np.arange(0.25, 2.5, 0.001)
        q = np.zeros(len(lamb), float)
        res = []
        for k in range(-potent, potent):
            q = q + f*np.exp(-2.0*(k**2)*(lamb**2))
            f = -f
        for i in range(len(lamb)):
            res.append((lamb[i], q[i]))
        return res

    def ksfunc(lamb):
        f = 1
        potent = 10000
        q = 0
        for k in range(-potent, potent):
            q = q + f*np.exp(-2.0*(k**2)*(lamb**2))
            f *= -1
        return q

    def edf(dg_data):
        edf_ = []
        ndata = []
        data = deepcopy(dg_data)
        data.sort()
        N = float(len(data))
        cnt = 0
        for item in data:
            cnt += 1
            edf_.append(cnt/N)
            ndata.append(item)
        ndata = np.array(ndata)
        edf_ = np.array(edf_)
        return ndata, edf_

    def cdf(dg_data):
        data = deepcopy(dg_data)
        data.sort()
        mean = np.average(data)
        sig = np.std(data)
        cdf = 0.5*(1+erf((data-mean)/float(sig*np.sqrt(2))))
        return cdf

    N = len(data)
    nd, ed = edf(data)
    cd = cdf(data)
    siglev = 1-alpha
    dval = []
    for i, val in enumerate(ed):
        d = abs(val-cd[i])
        dval.append(d)
        if i:
            d = abs(ed[i-1]-cd[i])
            dval.append(d)
    dmax = max(dval)
    check = np.sqrt(N)*dmax
    if not refks:
        refks = ksref()
    lst = [x for x in refks if x[1] > siglev]
    lam0 = lst[0][0]
    if check >= lam0:
        bOk = False
    else:
        bOk = True

    q = ksfunc(check)
    return (1-q), lam0, check, bOk


def data2gauss(data):
    '''Takes a one dimensional array and fits a Gaussian.

    Returns
    -------
    float
        mean of the distribution.
    float
        standard deviation of the distribution.
    float
        height of the curve's peak.
    '''
    m = np.average(data)
    dev = np.std(data)
    A = 1./(dev*np.sqrt(2*np.pi))
    return m, dev, A
