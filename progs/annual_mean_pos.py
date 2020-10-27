#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: annual_mean_pos.py
"""
Created on Mon May  4 16:31:15 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, setdiff
import numpy as np
from read_ts import get_ts


# -----------------------------  FUNCTIONS -----------------------------
def calc_wmean(x, err):
    """Calculate the weighted mean.
    """

    if len(x) == 1:
        wmean,  wmerr = x[0], err[0]
    else:
        wmean = np.dot(x, err**-2) / np.dot(1/err, 1/err)
        wmerr = np.sum(1/err) / np.sum(err**-2)

    return wmean, wmerr


# Source list
sou = Table.read("../data/sou_list.txt", format="ascii")

# Output
f_out = open("../data/yearly-mean-dpos.txt", "w")

print("iers_name, year, dra, dra_err, ddec, ddec_err, num_ses", file=f_out)

N = len(sou)
i = 0

for soui in sou["iers_name"]:

    print("Processing source {:s} [{:d}/{:d}]".format(soui, i+1, N))
    ts = get_ts(soui, "../data/ts")

    # Eliminate data points that have ratio of offset to formal larger than 10 or offset > 10 mas
    X_ra = ts["dra"] / ts["ra_err"]
    X_dec = ts["ddec"] / ts["dec_err"]

    mask = ((X_ra <= 10) & (ts["dra"] < 10) & (X_dec <= 10) & (ts["ddec"] < 10))
    ts1 = ts[mask]

    # Caculate the annual mean position
    for year in range(1979, 2020):
        ts2 = ts1[((ts1["mjy"] >= year) & (ts1["mjy"] < year+1))]

        if len(ts2):
            dra0, dra0_err = calc_wmean(ts1["dra"], ts1["ra_err"])
            ddc0, ddc0_err = calc_wmean(ts1["ddec"], ts1["dec_err"])
            ra, dec = ts1["ra"][0], ts1["dec"][0]

            print("{:8s}, {:.0f}, {:+7.3f},  {:7.3f},  {:+7.3f},  {:7.3f},  "
                  "{:7.3f},  {:+7.3f}, {:3d}".format(
                      soui, year, dra0, dra0_err, ddc0, ddc0_err, ra, dec, len(ts2)), file=f_out)

    i += 1

    #### TO-DO ####
    # Allan Variance
    #### TO-DO ####

f_out.close()
# --------------------------------- END --------------------------------
