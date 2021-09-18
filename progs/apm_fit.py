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


# --------------------------------- MAINS --------------------------------
ts_type = "nju_glo_8step"
max_sigma = 3

if ts_type == "opa_ind":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou"
    stats_file = "../logs/ts_stat_20210315"
    pm_file = "../data/ts_pm_fit_{}sigma.dat".format(max_sigma)
elif ts_type == "opa_glo":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou-from-glo"
    stats_file = "../logs/ts_stat_glo_20210315.log"
    pm_file = "../data/ts_glo_pm_fit_{}sigma.dat".format(max_sigma)
elif ts_type == "nju_glo_4step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-4step"
    stats_file = "../logs/ts_stat_nju_20210531.log"
    pm_file = "../data/ts_nju_pm_fit_{}sigma.dat".format(max_sigma)
elif ts_type == "nju_glo_8step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-8step"
    stats_file = "../logs/ts_stat_nju_20210605.log"
    pm_file = "../data/ts_nju_pm_fit_{}sigma-8step.dat".format(max_sigma)
elif ts_type == "nju_glo_10step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-10step"
    stats_file = "../logs/ts_stat_nju_20210608.log"
    pm_file = "../data/ts_nju_pm_fit_{}sigma-10step.dat".format(max_sigma)
elif ts_type == "nju_glo_20step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-20step"
    stats_file = "../logs/ts_stat_nju_20210609.log"
    pm_file = "../data/ts_nju_pm_fit_{}sigma-20step.dat".format(max_sigma)


with open(pm_file, "w") as f_sta:

    print("iers_name, num_cln, num_outl, outl_ra, outl_dec, ra, dec, "
          "pmra, pmdec, pmra_err, pmdec_err, pmra_pmdec_cor", file=f_sta)

    ts_stats = Table.read(stats_file, format="ascii")

    # First, we removed all sources with num_ses < 5
    mask = (ts_stats["num_ses"] >= 5)
    ts_stats = ts_stats[mask]

    N0 = len(ts_stats)

    for i, (soui, rai, deci) in enumerate(ts_stats["iers_name", "mean_ra", "mean_dec"]):
        print("Processing source {:s} [{:4d}/{:4d}]".format(soui, i+1, N0))

        dataline = []

        # Add source name to data line
        dataline.append("{:8s}".format(soui))

        ts = get_ts(soui, ts_dir)

        # Only keep data point before 2021
        ts = ts[ts["jyear"] < 2021]
        fac = np.cos(np.deg2rad(deci))
        dra = (ts["ra"] - rai) * 3.6e6 * fac  # Degree -> mas
        ddec = (ts["dec"] - deci) * 3.6e6  # Degree -> mas
        ts.add_columns([dra, ddec], names=["dra", "ddec"])
        ts["ra_err"] = ts["ra_err"] * fac

        # Eliminate data points that have ratio of offset to formal larger than 10 or offset > 10 mas
        X_ra = dra / ts["ra_err"]
        X_dec = ddec / ts["dec_err"]

        lim = max_sigma

        # R.A.
        mask1 = (np.fabs(X_ra) <= lim)
        out_ra = len(ts) - len(ts[mask1])

        # declination
        mask2 = (np.fabs(X_dec) <= lim)
        out_dec = len(ts) - len(ts[mask2])

        ts1 = ts[mask1 & mask2]
        num_cln = len(ts1)
        out_tot = len(ts) - num_cln

        dataline.append(",{:4d},{:4d},{:4d},{:4d}".format(num_cln, out_tot, out_ra, out_dec))
        dataline.append(",{:8.3f},{:8.3f}".format(rai, deci))

        if num_cln < 5:
            dataline.append(",,,,,")
        else:
            # Apparent proper motion fit
            # Estimate pm_ra and pm_dec
            res1 = linfit2d(ts1["jyear"]-np.median(ts1["jyear"]),
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

        print("".join(dataline), file=f_sta)


# --------------------------------- END --------------------------------
