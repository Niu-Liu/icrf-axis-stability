#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: tool_func.py
"""
Created on Thu Apr 23 17:39:40 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Some tool functions

"""

from astropy.table import Table
from astropy import units as u
# scikit-learn bootstrap
from sklearn.utils import resample
import numpy as np
from numpy.random import choice

# My progs
from my_progs.catalog.read_icrf import read_icrf3, read_icrf2
from my_progs.catalog.pos_diff import radio_cat_diff_calc
from my_progs.catalog.vsh_deg1_cor import vsh_deg01_fitting


# -----------------------------  FUNCTIONS -----------------------------
def sample_extraction(pos_oft, sam_size=500, with_replace=True):
    """Randomly select items of some number.
    """

    N = len(pos_oft)
    ind = resample(np.arange(N), replace=with_replace, n_samples=sam_size)

    pos_oft_sub = pos_oft[ind]

    return pos_oft_sub


def sample_clean(pos_oft, rho0=10, print_log=False):
    """Outlier elimination for VSH fitting.
    """

    # Remove the outlier (consider the normalized separation)
    N0 = len(pos_oft)
    # X0 = np.sqrt(np.log(N0) * 2)
    X0 = 5

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


def vsh_fitting(pos_oft, print_log=False):
    """Estimate the VSH coefficient.
    """

    # Transform astropy.table.Column into np.array
    dra = np.array(pos_oft["dra"])
    ddec = np.array(pos_oft["ddec"])
    dra_err = np.array(pos_oft["dra_err"])
    ddec_err = np.array(pos_oft["ddec_err"])
    ra_rad = np.array(pos_oft["ra"].to(u.radian))
    dec_rad = np.array(pos_oft["dec"].to(u.radian))
    dra_ddec_cov = np.array(pos_oft["dra_ddec_cov"])

    # Transformation parameters
    pmt, sig, corr = vsh_deg01_fitting(
        dra, ddec, ra_rad, dec_rad, dra_err, ddec_err,
        cov=dra_ddec_cov, elim_flag="None", fit_type="rotation")

    # mas -> uas
    pmt = pmt * 1.e3
    sig = sig * 1.e3

    # Print results
    if print_log:
        print("Estimates (%6d sources)\n"
              "----------------------------------------------\n"
              "               Rotation [uas]                 \n"
              "       x             y             z          \n"
              "----------------------------------------------\n"
              "  %+4.0f +/- %3.0f  %+4.0f +/- %3.0f  %+4.0f +/- %3.0f\n"
              "----------------------------------------------\n" %
              (dra.size,
               pmt[0], sig[0], pmt[1], sig[1], pmt[2], sig[2]))

    return pmt, sig


def orientation_estimate(pos_oft, opt):
    """Estimate the orientation offsets
    """

    pmti, sigi = vsh_fitting(pos_oft)

    opti = np.hstack((pmti, sigi))

    opt = np.vstack((opt, opti))

    return opt


def orientation_estimate_cln(pos_oft, opt):
    """Estimate the orientation offsets after removing outliers
    """

    pos_oft_cln = sample_clean(pos_oft, rho0=10)

    pmti, sigi = vsh_fitting(pos_oft_cln)

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
        pos_oft_sub = sample_extraction(pos_oft, sam_size, with_replace)
        # Sample size is around 60% of the common sources, 3410 * 0.6 = 2046

        opt = orientation_estimate(pos_oft_sub, opt)
        opt1 = orientation_estimate_cln(pos_oft_sub, opt1)
        print(" done!")

    return opt, opt1


def save_data(data, fname):
    tab = Table(data, names=["R1", "R2", "R3", "R1_err", "R2_err", "R3_err"])
    tab.write(fname, overwrite=True)


def calc_orient_4_yearly_crf(pos_oft):
    """
    """

    N0 = len(pos_oft)
    pos_oft_cln = sample_clean(pos_oft, rho0=10)
    N1 = len(pos_oft_cln)
    pmt, sig = vsh_fitting(pos_oft_cln)

    return N0, N1, pmt, sig


def sample_clean1(pos_oft, rho0=10, print_log=False):
    """Outlier elimination for VSH fitting.
    """

    # Remove the outlier (consider the normalized separation)
    N0 = len(pos_oft)
    # X0 = np.sqrt(np.log(N0) * 2)
    X0 = 5

    ang_sep = np.sqrt(pos_oft["dra"]**2 + pos_oft["ddec"]**2)
    ang_sep_err = np.sqrt(pos_oft["dra_err"]**2 + pos_oft["ddec_err"]**2)
    nor_sep = ang_sep / ang_sep_err

    mask = ((nor_sep <= X0)
            & (ang_sep < rho0))

    # Table of a clean sample
    pos_oft_cln = pos_oft[mask]
    N1 = len(pos_oft_cln)

    if print_log:
        print("For a sample of %d sources, "
              "the number of the outlier is smaller than 1 when X >= %.2f." % (N0, X0))
        print("After elimination, there are %d sources in the clean sample." % N1)
        print("The outlier rate is %.2f%%.\n" % ((N0 - N1) / N0 * 100))

    return pos_oft_cln


def vsh_fitting1(pos_oft, print_log=False):
    """Estimate the VSH coefficient.
    """

    # Transform astropy.table.Column into np.array
    dra = np.array(pos_oft["dra"])
    ddec = np.array(pos_oft["ddec"])
    dra_err = np.array(pos_oft["dra_err"])
    ddec_err = np.array(pos_oft["ddec_err"])
    ra_rad = np.deg2rad(np.array(pos_oft["ra"]))
    dec_rad = np.deg2rad(np.array(pos_oft["dec"]))

    # Transformation parameters
    pmt, sig, corr = vsh_deg01_fitting(
        dra, ddec, ra_rad, dec_rad, dra_err, ddec_err,
        cov=None, elim_flag="None", fit_type="rotation")

    # mas -> uas
    pmt = pmt * 1.e3
    sig = sig * 1.e3

    return pmt, sig


def calc_orient_4_yearly_crf1(pos_oft):
    """
    """

    N0 = len(pos_oft)
    pos_oft_cln = sample_clean(pos_oft, rho0=10)
    N1 = len(pos_oft_cln)
    pmt, sig = vsh_fitting(pos_oft_cln)

    return N0, N1, pmt, sig
# --------------------------------- END --------------------------------
