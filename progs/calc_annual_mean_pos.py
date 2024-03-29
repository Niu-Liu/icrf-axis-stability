#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: calc_annual_mean_pos.py
"""
Created on Mon May  4 16:31:15 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table
import numpy as np

from myprogs.vlbi.ts_func import get_ts


# -----------------------------  FUNCTIONS -----------------------------
def calc_wmean(x, err):
    """Calculate the weighted mean.
    """

    if len(x) == 1:
        wmean,  wmerr = x[0], err[0]
    else:
        wmean = np.dot(x, err**-2) / np.sum(err**-2)
        wmerr = np.sum(1/err) / np.sum(err**-2)

    return wmean, wmerr


# ts_type = "nju_glo_20step"
ts_type = "nju_glo_10step"
# ts_type = "nju_glo_8step"
# ts_type = "nju_glo_4step"

if ts_type == "opa_ind":
    sou_list = "../data/sou-ts.list"
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou"
    pm_file = "../data/ts_pm_fit_3sigma.dat"
    mean_pos_file = "../data/yearly-mean-position-from-ts.txt"
    crf_dir = "../data/yearly-ts"
elif ts_type == "opa_glo":
    sou_list = "../data/sou-ts-glo.list"
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou-from-glo"
    mean_pos_file = "../data/yearly-mean-position-from-ts-glo.txt"
    crf_dir = "../data/yearly-ts-glo"
elif ts_type == "nju_glo_4step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-4step"
    sou_list = "{}/sou_list_4_ts.txt".format(ts_dir)
    mean_pos_file = "../data/yearly-mean-position-from-ts-nju.txt"
    crf_dir = "../data/yearly-ts-nju"
elif ts_type == "nju_glo_8step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-8step"
    sou_list = "{}/sou_list_4_ts.txt".format(ts_dir)
    mean_pos_file = "../data/yearly-mean-position-from-ts-nju-8step.txt"
    crf_dir = "../data/yearly-ts-nju-8step"
elif ts_type == "nju_glo_10step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-10step"
    sou_list = "{}/sou_list_4_ts.txt".format(ts_dir)
    mean_pos_file = "../data/yearly-mean-position-from-ts-nju-10step.txt"
    crf_dir = "../data/yearly-ts-nju-10step"
elif ts_type == "nju_glo_20step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-20step"
    sou_list = "{}/sou_list_4_ts.txt".format(ts_dir)
    mean_pos_file = "../data/yearly-mean-position-from-ts-nju-20step.txt"
    crf_dir = "../data/yearly-ts-nju-20step"

# Source list
sou = Table.read(sou_list, format="ascii")

# Output
f_out = open(mean_pos_file, "w")

print("iers_name, year, num_ses, ra, ra_err, dec, dec_err, ra_dec_corr", file=f_out)

N = len(sou)
i = 0

# Annually average the source positions
for soui in sou["iers_name"]:

    print("Processing source {:s} [{:4d}/{:4d}]".format(soui, i+1, N))

    ts = get_ts(soui, ts_dir)

    # Caculate the annual mean position
    for year in range(1979, 2021):
        ts2 = ts[((ts["jyear"] >= year) & (ts["jyear"] < year+1))]

        if len(ts2) >= 1:
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
pos_tab = Table.read(mean_pos_file, format="ascii")
for year in np.arange(1979, 2021):
    mask = (pos_tab["year"] == year)
    sub_tab = pos_tab[mask]

    if len(sub_tab):
        sub_tab.write("{:s}/{:d}.dat".format(crf_dir, year),
                      format="ascii.csv", overwrite=True)

# --------------------------------- END --------------------------------
