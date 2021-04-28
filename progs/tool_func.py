#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: tool_func.py
"""
Created on Thu Apr 23 17:39:40 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Some tool functions. The comment will be added when I am free.

"""


import numpy as np
# scikit-learn bootstrap
from sklearn.utils import resample

from astropy.table import Table
from astropy import units as u

# My progs
from my_progs.catalog.read_icrf import read_icrf3, read_icrf2
from my_progs.catalog.pos_diff import radio_cat_diff_calc

from my_progs.vsh.vsh_fit import rotgli_fit_4_table


# -----------------------------  FUNCTIONS -----------------------------
def rot_calc(omega, err):
    """Calculate the module of a vector
    """

    w = np.sqrt(np.sum(omega**2))
    w_err = np.sqrt(np.sum(err**2))

    return w, w_err


def bootstrap_sample(tab, sam_size=500, with_replace=True):
    """Randomly select items of some number.
    """

    N = len(tab)
    ind = resample(np.arange(N), replace=with_replace, n_samples=sam_size)
    new_tab = tab[ind]

    return new_tab


def sample_clean(pos_oft, rho0=10, print_log=False):
    """Outlier elimination for VSH fitting.
    """

    # Remove the outlier (consider the normalized separation)
    N0 = len(pos_oft)
    X0 = np.sqrt(np.log(N0) * 2)

    mask = ((pos_oft["nor_sep"] <= X0)
            & (pos_oft["ang_sep"] < rho0))

    # Table of a clean sample
    pos_oft_cln = pos_oft[mask]
    N1 = len(pos_oft_cln)

    if print_log:
        print("For a sample of %d sources, "
              "the number of the outlier is smaller than 1 when X >= %.2f." % (N0, X0))
        print("After elimination, there are %d sources in the clean sample." % N1)
        print("The outlier rate is %.2f%%.\n" % ((N0 - N1) / N0 * 100))

    return pos_oft_cln


def vsh_fit_for_pos(pos_oft, print_log=False):
    """Estimate the VSH coefficient.
    """

    output = rotgli_fit_4_table(pos_oft, verbose=False)

    # Keep rotation only and mas -> uas
    pmt = output["pmt1"][3:] * 1.e3
    sig = output["sig1"][3:] * 1.e3

    # Print resultss
    if print_log:
        print("Estimates (%6d sources)\n"
              "----------------------------------------------\n"
              "               Rotation [uas]                 \n"
              "       x             y             z          \n"
              "----------------------------------------------\n"
              "  %+4.0f +/- %3.0f  %+4.0f +/- %3.0f  %+4.0f +/- %3.0f\n"
              "----------------------------------------------\n" %
              (len(pos_oft),
               pmt[0], sig[0], pmt[1], sig[1], pmt[2], sig[2]))

    return pmt, sig


def calc_orient(pos_oft):
    """Calculate orientation angles based on positional difference
    """

    N0 = len(pos_oft)
    pos_oft_cln = sample_clean(pos_oft, rho0=10)
    N1 = len(pos_oft_cln)
    pmt, sig = vsh_fit_for_pos(pos_oft_cln)

    return N0, N1, pmt, sig


def orientation_estimate(pos_oft, opt, clean=False):
    """Estimate the orientation offsets
    """

    if clean:
        pos_oft = sample_clean(pos_oft, rho0=10)

    pmti, sigi = vsh_fit_for_pos(pos_oft)
    opti = np.hstack((pmti, sigi))
    opt = np.vstack((opt, opti))

    return opt


def orient_angle_sampling(pos_oft, iter_time=1000, sam_size=2000, with_replace=False):
    """Orientation angle sampling.
    """

    # Array to store data in the form of [r1, r2, r3, sig1, sig2, sig3]
    opt = np.empty(dtype=np.float, shape=(6,))  # For all sources
    opt1 = np.empty(dtype=np.float, shape=(6,))  # For a clean sample

    for i in range(iter_time):
        print(">>>>{:4d}th iteration:".format(i+1), end="")
        pos_oft_sub = bootstrap_sample(pos_oft, sam_size, with_replace)
        # Sample size is around 60% of the common sources, 3410 * 0.6 = 2046

        # All sources
        opt = orientation_estimate(pos_oft_sub, opt)

        # Removing Outliers
        opt1 = orientation_estimate(pos_oft_sub, opt1, True)
        print(" done!")

    return opt, opt1


def save_data(data, fname):
    tab = Table(data, names=["R1", "R2", "R3", "R1_err", "R2_err", "R3_err"])
    tab.write(fname, overwrite=True)


def vsh_fit_for_pm(apm_table):
    """Estimate VSH coefficients from apparent proper motion.
    """

    output = rotgli_fit_4_table(apm_table, verbose=False)

    # Keep rotation only
    pmt = output["pmt1"][3:]
    sig = output["sig1"][3:]

    # Calculate the total rotation
    w, w_err = rot_calc(pmt, sig)

    # Concatenate the results
    pmt = np.hstack((pmt, [w]))
    sig = np.hstack((sig, [w_err]))

    return pmt, sig


def vsh_fit_for_pm2(apm_table):
    """Estimate VSH coefficients from apparent proper motion.

    Only rotation vector is estimated.
    """

    output = rotgli_fit_4_table(apm_table, fit_type="T", verbose=False)

    # Keep rotation only
    pmt = output["pmt1"]
    sig = output["sig1"]

    # Calculate the total rotation
    w, w_err = rot_calc(pmt, sig)

    # Concatenate the results
    pmt = np.hstack((pmt, [w]))
    sig = np.hstack((sig, [w_err]))

    return pmt, sig

# --------------------------------- END --------------------------------
