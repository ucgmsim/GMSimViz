"""
Acceleration amplification models.

@date 24 June 2016
@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz

Implemented Models
==============================
cb_amp (version = "2008"):
    Based on Campbell and Bozorgnia 2008 - added 24 June 2016
cb_amp (version = "2014"):
    Based on Campbell and Bozorgnia 2014 - added 22 September 2016

Usage
==============================
from siteamp_models import cb_amp (or *)
cb_amp(variables, ...)
"""

# math functions faster than numpy for non-vector data
from math import ceil, exp, log

import numpy as np
# divide by zero = inf
np.seterr(divide = 'ignore')

def cb_amp(dt, n, vref, vsite, vpga, pga, version = '2008', \
        flowcap = 0.0, fmin = 0.01, fmidbot = 0.2, fmid = 1.0, \
        fhigh = 3.33, fhightop = 10.0, fmax = 15.0):
    # cb08 constants
    n_per = 22
    scon_c = 1.88
    scon_n = 1.18
    per = np.array([0.00, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, \
            0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, \
            1.00, 1.50, 2.00, 3.00, 4.00, 5.00, 7.50, 10.0])
    m_freq = 1.0 / per
    f1_src = np.hstack(([1000.0], m_freq[1:]))
    if version == '2008':
        c10 = np.array([1.058, 1.058, 1.102, 1.174, 1.272, 1.438, 1.604, \
                1.928, 2.194, 2.351, 2.460, 2.587, 2.544, 2.133, \
                1.571, 0.406, -0.456, -0.82, -0.82, -0.82, -0.82, -0.82])
    elif version == '2014':
        c10 = np.array([1.090, 1.094, 1.149, 1.290, 1.449, 1.535, 1.615, \
                1.877, 2.069, 2.205, 2.306, 2.398, 2.355, 1.995, \
                1.447, 0.330, -0.514, -0.848, -0.793, -0.748, -0.664, -0.576])
    else:
        raise Exception('BAD CB AMP version specified.')
    k1 = np.array([865.0, 865.0, 865.0, 908.0, 1054.0, 1086.0, 1032.0, \
            878.0, 748.0, 654.0, 587.0, 503.0, 457.0, 410.0, \
            400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0])
    k2 = np.array([-1.186, -1.186, -1.219, -1.273, -1.346, -1.471, -1.624, \
            -1.931, -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, \
            -1.955, -1.025, -0.299, 0.0, 0.0, 0.0, 0.0, 0.0])

    # f_site function domains
    # TODO: normal functions for speed, will require testing
    fs_low = lambda T, vs30, a1100 : c10[T] * log(vs30 / k1[T]) \
            + k2[T] * log((a1100 + scon_c \
            * exp(scon_n * log(vs30 / k1[T]))) / (a1100 + scon_c))
    fs_mid = lambda T, vs30, a1100 = None : \
            (c10[T] + k2[T] * scon_n) * log(vs30 / k1[T])
    fs_high = lambda T, vs30 = None, a1100 = None : \
            (c10[T] + k2[T] * scon_n) * log(1100.0 / k1[T])
    fs_auto = lambda T, vs30 : \
            fs_low if vs30 < k1[T] else fs_mid if vs30 < 1100.0 else fs_high
    fs1100 = fs_high(0)

    # default amplification is 1.0 (keeping values the same)
    ampf = np.ones(n / 2, np.float)

    fs_vpga = fs_auto(0, vpga)(0, vpga, pga)
    a1100 = pga * exp(fs1100 - fs_vpga)

    # calculate factor for each period
    it = (exp(fs_auto(T, vsite)(T, vsite, a1100) \
            - fs_auto(T, vref)(T, vref, a1100)) \
            for T in xrange(n_per))
    ampf0 = np.fromiter(it, np.float, count = n_per)

    try:
        # T is the first occurance of a value <= flowcap
        # throws IndexError if no results (the second [0])
        T = np.nonzero((m_freq[:-1] <= flowcap))[0][0]
        # T cannot be the last value because of the following logic
        ampf0[T + 1:] = ampf0[T]
    except IndexError:
        pass
    # frequencies of fourier transform
    ftfreq = np.arange(1, n / 2) * (1.0 / (n * dt))

    # calculate ampv based on period group
    j = n_per - 1
    f0 = f1_src[j]
    a0 = ampf0[j]
    f1 = f0
    a1 = a0
    dadf = 0.0
    for i, ftf in enumerate(ftfreq):
        if ftf > f1:
            f0 = f1
            a0 = a1
            if j - 1:
                j -= 1
                a1 = ampf0[j]
                f1 = f1_src[j]
                dadf = (a1 - a0) / log(f1 / f0)
            else:
                dadf = 0.0
        ampv = a0 + dadf * log(ftf / f0)

        # scale amplification factor by frequency
        if ftf < fmin:
            continue
        if ftf < fmidbot:
            ampf[i + 1] = 1.0 + log(ftf / fmin) \
                    * (ampv - 1.0) / log(fmidbot / fmin)
        elif ftf < fhightop:
            ampf[i + 1] = ampv
        elif ftf < fmax:
            ampf[i + 1] = ampv + log(ftf / fhightop) \
                    * (1.0 - ampv) / log(fmax / fhightop)

    return ampf
