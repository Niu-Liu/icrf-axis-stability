#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: calc_annual_mean_pos.py
"""
Created on Mon May  4 16:31:15 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table
import numpy as np
from my_progs.vlbi.ts_func import get_ts


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
sou = Table.read("../data/ts-sou.list", format="ascii")

# Output
f_out = open("../data/yearly-mean-position-from-ts.txt", "w")

print("iers_name, year, num_ses, ra, ra_err, dec, dec_err, ra_dec_corr", file=f_out)

N = len(sou)
i = 0

# Annually average the source positions
for soui in sou["iers_name"]:

    print("Processing source {:s} [{:4d}/{:4d}]".format(soui, i+1, N))
    ts = get_ts(soui)

    # Caculate the annual mean position
    for year in range(1979, 2021):
        ts2 = ts[((ts["jyear"] >= year) & (ts["jyear"] < year+1))]

        if len(ts2):
            ra0, ra0_err = calc_wmean(ts2["ra"], ts2["ra_err"])
            dc0, dc0_err = calc_wmean(ts2["dec"], ts2["dec_err"])
            ra_dec_cor = np.median(ts2["ra_dec_corr"])

            print("{:8s}, {:.0f}, {:3d}, "
                  "{:15.10f}, {:6.3f}, {:+15.10f}, {:6.3f}, {:6.3f}".format(
                      soui, year, len(ts2),
                      ra0, ra0_err, dc0, dc0_err, ra_dec_cor), file=f_out)
    i += 1

f_out.close()

# Generate yearly CRF
pos_tab = Table.read("../data/yearly-mean-position-from-ts.txt", format="ascii")
for year in np.arange(1979, 2021):
    mask = (pos_tab["year"] == year)
    sub_tab = pos_tab[mask]

    if len(sub_tab):
        sub_tab.write("../data/yearly-ts/{:d}.dat".format(year),
                      format="ascii.csv", overwrite=True)

# --------------------------------- END --------------------------------
