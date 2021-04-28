#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: apm_fit.py
"""
Created on Sat Apr 18 13:58:59 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from astropy.table import Table

# My progs
from my_progs.vlbi.ts_func import get_ts
from linear_fit import linfit2d

# -----------------------------  FUNCTIONS -----------------------------


# --------------------------------- MAINS --------------------------------
ts_stats = Table.read("../logs/ts_stat_20210315.log", format="ascii")

# First, we removed all sources with num_ses < 5
mask = (ts_stats["num_ses"] >= 5)
ts_stats = ts_stats[mask]

# Store results
f_sta = open("../logs/ts_pm_fit_4_test.dat", "w")
print("iers_name, "
      "num_cln, num_outl, outl_ra, outl_dec, "
      "pm_ra, pm_dec, pmra_err, pmdec_err, pmra_pmdec_corr, "
      "pm, pm_err, pa, pa_err", file=f_sta)

N0 = len(ts_stats)
i = 0

for soui, rai, deci in ts_stats["iers_name", "mean_ra", "mean_dec"]:
    print("Processing source {:s} [{:4d}/{:4d}]".format(soui, i+1, N0))

    dataline = []

    ts = get_ts(soui)
    fac = np.cos(np.deg2rad(deci))
    dra = (ts["ra"] - rai) * 3.6e6 * fac  # Degree -> mas
    ddec = (ts["dec"] - deci) * 3.6e6  # Degree -> mas
    ts.add_columns([dra, ddec], names=["dra", "ddec"])
    ts["ra_err"] = ts["ra_err"] * fac

    # Eliminate data points that have ratio of offset to formal larger than 10 or offset > 10 mas
    X_ra = dra / ts["ra_err"]
    X_dec = ddec / ts["dec_err"]

    # R.A.
    mask1 = ((np.fabs(X_ra) <= 10) & (np.fabs(ts["dra"]) < 10))
    out_ra = len(ts) - len(ts[mask1])

    # declination
    mask2 = ((np.fabs(X_dec) <= 10) & (np.fabs(ts["ddec"]) < 10))
    out_dec = len(ts) - len(ts[mask2])

    ts1 = ts[mask1 & mask2]
    num_cln = len(ts1)
    out_tot = len(ts) - num_cln

    if num_cln >= 1000:
        # Add source name to data line
        dataline.append("{:8s}".format(soui))
        dataline.append(",{:4d},{:4d},{:4d},{:4d}".format(
            num_cln, out_tot, out_ra, out_dec))

        # Apparent proper motion fit
        # Estimate pm_ra and pm_dec
        res1 = linfit2d(ts1["jyear"],
                        ts1["dra"],
                        ts1["ddec"],
                        x_err=ts1["ra_err"],
                        y_err=ts1["dec_err"],
                        xy_cor=ts1["ra_dec_corr"],
                        fit_type="sep")

        pm_ra, pmra_err = res1["x1"], res1["x1_err"]
        pm_dec, pmdec_err = res1["y1"], res1["y1_err"]
        pmra_pmdec_corr = res1["x1y1_cor"]

        dataline.append(",{:+10.5f},{:+10.5f},{:8.5f},{:8.5f},{:8.5f}".format(
            pm_ra, pm_dec, pmra_err, pmdec_err, pmra_pmdec_corr))

        # Estimate pm and pa
        res2 = linfit2d(ts1["jyear"],
                        ts1["dra"],
                        ts1["ddec"],
                        x_err=ts1["ra_err"],
                        y_err=ts1["dec_err"],
                        xy_cor=ts1["ra_dec_corr"],
                        fit_type="tot")

        pm, pm_err = res2["pm"], res2["pm_err"]
        pa, pa_err = res2["pa"], res2["pa_err"]

        dataline.append(",{:+10.5f},{:8.5f},{:+8.3f},{:8.3f}".format(
            pm, pm_err, pa, pa_err))

    print("".join(dataline), file=f_sta)

    i += 1

f_sta.close()
# --------------------------------- END --------------------------------