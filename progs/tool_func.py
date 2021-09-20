#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: tool_func.py
"""
Created on Thu Apr 23 17:39:40 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Some tool functions. The comment will be added when I am free.

"""


from my_progs.vsh.vsh_fit import rotgli_fit_4_table
from my_progs.catalog.pos_diff import radio_cat_diff_calc
from my_progs.catalog.read_icrf import read_icrf3, read_icrf2
from astropy.stats import sigma_clip, mad_std
from astropy import units as u
from astropy.table import Table
from sklearn.utils import resample
import numpy as np
# np.random.seed(28)


# -----------------------------  FUNCTIONS -----------------------------
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
    X0 = 10000000
    rho0 = 1000000

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
    pmt = output["pmt1"] * 1.e3
    sig = output["sig1"] * 1.e3

    # Calculate the total rotation
    w, w_err = output["R"] * 1.e3, output["R_err"] * 1.e3
    g, g_err = output["G"] * 1.e3, output["G_err"] * 1.e3

    # Concatenate the results
    pmt = np.hstack((pmt, [g, w]))
    sig = np.hstack((sig, [g_err, w_err]))

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
               pmt[3], sig[3], pmt[4], sig[4], pmt[5], sig[5]))

    return pmt, sig


def calc_orient(pos_oft):
    """Calculate orientation angles based on positional difference
    """

    N0 = len(pos_oft)
    pos_oft_cln = sample_clean(pos_oft, rho0=10)
    N1 = len(pos_oft_cln)
    pmt, sig = vsh_fit_for_pos(pos_oft_cln)

    return N0, N1, pmt, sig


def calc_orient_new(pos_oft):
    """Calculate orientation angles based on positional difference
    """

    N0 = len(pos_oft)
    pmt, sig = vsh_fit_for_pos(pos_oft)

    return N0, pmt, sig


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

    # Array to store data in the form of
    # [g1, g2, g3, r1, r2, r3, g, r,
    # sig_r1, sig_r2, sig_r3, sig_g1, sig_g2, sig_g3, sig_g, sig_r]
    opt = np.empty(dtype=np.float, shape=(16,))  # For all sources
    opt1 = np.empty(dtype=np.float, shape=(16,))  # For a clean sample

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
    tab = Table(data, names=["G1", "G2", "G3", "R1", "R2", "R3",
                             "G", "R",
                             "G1_err", "G2_err", "G3_err",
                             "R1_err", "R2_err", "R3_err",
                             "G_err", "R_err"])
    tab.write(fname, overwrite=True)


def vsh_fit_for_pm(apm_table):
    """Estimate VSH coefficients from apparent proper motion.
    """

    output = rotgli_fit_4_table(apm_table, verbose=False)

    # Keep rotation only
    pmt = output["pmt1"]
    sig = output["sig1"]

    # Calculate the total rotation
    w, w_err = output["R"], output["R_err"]
    g, g_err = output["G"], output["G_err"]

    # Concatenate the results
    pmt = np.hstack((pmt[3:], [w], pmt[:3], [g]))
    sig = np.hstack((sig[3:], [w_err], sig[:3], [g_err]))

    return pmt, sig, output


def vsh_fit_for_pm2(apm_table):
    """Estimate VSH coefficients from apparent proper motion.

    Only rotation vector is estimated.
    """

    output = rotgli_fit_4_table(apm_table, fit_type="T", verbose=False)

    # Keep rotation only
    pmt = output["pmt1"]
    sig = output["sig1"]

    # Calculate the total rotation
    w, w_err = output["R"], output["R_err"]

    # Concatenate the results
    pmt = np.hstack((pmt, [w]))
    sig = np.hstack((sig, [w_err]))

    return pmt, sig, output


def calc_mean_std(y):
    """Esimate robustly the mean and standard deviation
    """

    filtered_data = sigma_clip(y, sigma=3, maxiters=1, stdfunc=mad_std)
    ymean, ystd = np.mean(filtered_data), np.std(filtered_data)

    return ymean, ystd


def random_walk(epoch, t_scale=5, sigma_var=2):
    """
    """

    dt = epoch[1:] - epoch[:-1]
    dt = np.concatenate(([0], dt))

    # Positional offset
    dra = np.zeros(len(epoch))
    ddec = np.zeros(len(epoch))

    dra[0] = np.random.random() - 0.5
    ddec[0] = np.random.random() - 0.5

    for i in range(len(epoch)):
        # Exponential factor
        exp_fac_i = np.exp(-dt[i]/t_scale)

        # Gaussian factor
        sigma_i = sigma_var * np.sqrt(1-np.exp(-2*dt[i]/t_scale))
        g_ra_i = (np.random.random_sample()-0.5) * sigma_i
        g_dec_i = (np.random.random_sample()-0.5) * sigma_i

        dra[i] = exp_fac_i * dra[i-1] + g_ra_i
        ddec[i] = exp_fac_i * ddec[i-1] + g_dec_i

    return dra, ddec

# --------------------------------- END --------------------------------
