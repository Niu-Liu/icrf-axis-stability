#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: ts_stat.py
"""
Created on Sat Apr 18 13:58:59 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import cos, sin, tan
from astropy.table import Table
import statsmodels.api as sm

# My progs
from myprogs.StatsFunc.simple_func import wgt_mean_calc
from myprogs.vlbi.ts_func import get_ts


# --------------------------------- MAINS --------------------------------

# ts_type = "nju_glo_4step"
# ts_type = "nju_glo_8step"
ts_type = "nju_glo_10step"
# ts_type = "nju_glo_20step"

if ts_type == "opa_ind":
    sou_list = "../data/sou-ts.list"
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou"
    stats_file = "../logs/ts_stat_20210315"
elif ts_type == "opa_glo":
    sou_list = "../data/sou-ts-glo.list"
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou-from-glo"
    stats_file = "../logs/ts_stat_glo_20210315.log"
elif ts_type == "nju_glo_4step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-4step"
    sou_list = "{}/sou_list_4_ts.txt".format(ts_dir)
    stats_file = "../logs/ts_stat_nju_20210531.log"
elif ts_type == "nju_glo_8step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-8step"
    sou_list = "{}/sou_list_4_ts.txt".format(ts_dir)
    stats_file = "../logs/ts_stat_nju_20210605.log"
elif ts_type == "nju_glo_10step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-10step"
    sou_list = "{}/sou_list_4_ts.txt".format(ts_dir)
    stats_file = "../logs/ts_stat_nju_20210608.log"
elif ts_type == "nju_glo_20step":
    ts_dir = "/Users/Neo/Astronomy/data/vlbi/nju/series-20step"
    sou_list = "{}/sou_list_4_ts.txt".format(ts_dir)
    stats_file = "../logs/ts_stat_nju_20210609.log"

# Output
with open(stats_file, "w") as f_sta:

    print("iers_name, num_ses, obs_len, mean_epoch, beg_epoch, end_epoch, "
          "mean_ra, sigma_ra, mean_dec, sigma_dec, med_err_ra, med_err_dec",
          file=f_sta)

    sou = Table.read(sou_list, format="ascii")
    N0 = len(sou)

    for i, soui in enumerate(sou["iers_name"]):

        print("Processing source {:s} [{:5d}/{:5d}]".format(soui, i+1, N0))
        ts = get_ts(soui, ts_dir)

        data_line = []

        # Number of sessions
        N = len(ts)

        # length of observing time span
        tmax, tmin = ts["jyear"].max(), ts["jyear"].min()
        dti = tmax - tmin
        mean_epo = ts["jyear"].mean()

        # Add to data line
        data_line.append("{:8s},{:4d},{:5.1f},{:8.2f},{:8.2f},{:8.2f}".format(
            soui, N, dti, mean_epo, tmin, tmax))

        # # Mean position
        # # Apply a global rotation
        # ra_rad = np.deg2rad(ts["ra"])
        # dec_rad = np.deg2rad(ts["dec"])
        #
        # r1 = 36e-6 / 3600
        # r2 = -52e-6 / 3600
        # r3 = -10e-6 / 3600
        # dra = (-r1*cos(ra_rad) - r2*sin(ra_rad)) * tan(dec_rad) + r3
        # ddec = r1 * sin(ra_rad) - r2 * cos(ra_rad)
        #
        # # RA
        # mean_ra, std_ra = wgt_mean_calc(ts["ra"]-dra, ts["ra_err"])
        #
        # # declination
        # mean_dec, std_dec = wgt_mean_calc(ts["dec"]-ddec, ts["dec_err"])

        # RA
        mean_ra, std_ra = wgt_mean_calc(ts["ra"], ts["ra_err"])

        # declination
        mean_dec, std_dec = wgt_mean_calc(ts["dec"], ts["dec_err"])

        # Degree -> mas
        std_ra = std_ra * 3.6e6
        std_dec = std_dec * 3.6e6
        data_line.append(",{:15.10f},{:6.3f},{:+15.10f},{:6.3f}".format(
            mean_ra, std_ra, mean_dec, std_dec))

        # Median formal uncertainty
        med_err_ra = np.median(ts["ra_err"])
        med_err_dec = np.median(ts["dec_err"])
        data_line.append(",{:6.3f},{:6.3f}".format(med_err_ra, med_err_dec))

        print("".join(data_line), file=f_sta)


# --------------------------------- END --------------------------------
