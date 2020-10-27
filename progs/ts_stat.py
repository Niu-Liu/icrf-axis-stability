#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: ts_stat.py
"""
Created on Sat Apr 18 13:58:59 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from read_ts import get_ts
import numpy as np
from astropy.table import Table
import statsmodels.api as sm

# My progs
from my_progs.stat_func.rms_calc import rms_calc


# -----------------------------  FUNCTIONS -----------------------------
def linear_fit(mjy, x, x_err):
    """A simple linear regression function
    """

    mjy = np.array(mjy-2015)
    mjy = sm.add_constant(mjy)

    mod_wls = sm.WLS(x, mjy, weights=x_err**-2)
    res_wls = mod_wls.fit()
    x0, dx_dt = res_wls.params
    x0_err, dx_dt_err = res_wls.bse

    return x0, x0_err, dx_dt, dx_dt_err


sou = Table.read("../data/sou_list.txt", format="ascii")

# A special check on the ICRF3 defining sources
def_sou = Table.read("../data/icrf3sx-def-sou.txt", format="ascii")

# Store results
f_sta = open("../logs/ts_stat.log", "w")
print("iers_name, num_ses, obs_len, sou_type, "
      "num_ses1, obs_len1, wmean_ra, wrms_ra, slp_ra, slp_ra_err, ra, "
      "num_ses2, obs_len2, wmean_dc, wrms_dc, slp_dc, slp_dc_err, dec", file=f_sta)

N0 = len(sou)
i = 0

for soui in sou["iers_name"]:

    print("Processing source {:s} [{:d}/{:d}]".format(soui, i+1, N0))
    ts = get_ts(soui, "../data/ts")

    dataline = []

    # Defining source group
    if soui in def_sou["iers_name"]:
        type = "D"
    else:
        type = " "

    # Number of sessions
    N = len(ts)

    # length of observing time span
    dti = ts["mjy"].max() - ts["mjy"].min()

    # Add to data line
    dataline.append("{:8s},{:4d},{:5.1f},{:s}".format(soui, N, dti, type))

    # Eliminate data points that have ratio of offset to formal larger than 10 or offset > 10 mas
    X_ra = ts["dra"] / ts["ra_err"]
    X_dec = ts["ddec"] / ts["dec_err"]

    # R.A.
    mask1 = ((X_ra <= 10) & (ts["dra"] < 10))
    ts1 = ts[mask1]

    N1 = len(ts1)

    if N1 < 2:
        dataline.append(",,,,,,,")
    else:
        # len of observing window
        dti1 = ts1["mjy"].max() - ts1["mjy"].min()

        # WRMS, WMean (mas)
        wmean_ra, wrms_ra, std_ra = rms_calc(ts1["dra"], ts1["ra_err"])
        ra = ts1["ra"][0]

        if N1 >= 5 and dti1 > 10:
            # Slope, Bias (mas/yr, mas)
            ra0, ra0_err, slp_ra, slp_ra_err = linear_fit(ts1["mjy"], ts1["dra"], ts1["ra_err"])
            dataline.append(",{:4d},{:5.1f},{:+6.3f},{:+6.3f},{:+6.3f},{:5.3f},{:5.3f}".format(
                N1, dti1, wmean_ra, wrms_ra, slp_ra, slp_ra_err, ra))
        else:
            dataline.append(",{:4d},{:5.1f},{:+6.3f},{:+6.3f},,{:5.3f}".format(
                N1, dti1, wmean_ra, wrms_ra, ra))

    # declination
    mask2 = ((X_dec <= 10) & (ts["ddec"] < 10))
    ts2 = ts[mask2]

    N2 = len(ts2)

    if N2 < 2:
        dataline.append(",,,,,,,")
    else:
        # len of observing window
        dti2 = ts2["mjy"].max() - ts2["mjy"].min()

        # WRMS, WMean (mas)
        wmean_dec, wrms_dec, std_dec = rms_calc(ts2["ddec"], ts2["dec_err"])
        dec = ts2["dec"][0]

        if N2 >= 5 and dti2 > 10:
            # Slope, Bias (mas/yr, mas)
            dec0, dec0_err, slp_dec, slp_dec_err = linear_fit(
                ts2["mjy"], ts2["ddec"], ts2["dec_err"])
            dataline.append(",{:4d},{:5.1f},{:+6.3f},{:+6.3f},{:+6.3f},{:5.3f},{:5.3f}".format(
                N2, dti2, wmean_dec, wrms_dec, slp_dec, slp_dec_err, dec))
        else:
            dataline.append(",{:4d},{:5.1f},{:+6.3f},{:+6.3f},,,{:5.3f}".format(
                N2, dti2, wmean_dec, wrms_dec, dec))

    #### TO-DO ####
    # Allan Variance
    #### TO-DO ####

    print("".join(dataline), file=f_sta)

    i += 1

f_sta.close()
# --------------------------------- END --------------------------------
