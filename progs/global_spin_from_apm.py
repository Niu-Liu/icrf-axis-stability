#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: global_spin_from_apm.py
"""
Created on Tue Apr 20 23:19:28 2021

@author: Neo(niu.liu@nju.edu.cn)

Estimate the global spin of the CRF based on the apparent proper motion of
radio sources.
"""

import numpy as np
from astropy.table import Table, join
from tool_func import vsh_fit_for_pm, bootstrap_sample

np.random.seed(28)


# -----------------------------  MAIN -----------------------------
# ts_type = "opa_glo"
# ts_type = "opa_ind"
# ts_type = "nju_glo_4step"
# ts_type = "nju_glo_8step"
ts_type = "nju_glo_10step"
# ts_type = "nju_glo_20step"

if ts_type == "opa_ind":
    pm_file = "../logs/ts_pm_fit_3sigma.dat"
    spin_file = "../logs/spin_fit.txt"
elif ts_type == "opa_glo":
    g
    pm_file = "../logs/ts_glo_pm_fit_3sigma.dat"
    spin_file = "../logs/spin_fit_from_glo_ts.txt"
elif ts_type == "nju_glo_4step":
    pm_file = "../data/ts_nju_pm_fit_10sigma.dat"
    spin_file = "../logs/spin_fit_from_nju_ts.txt"
elif ts_type == "nju_glo_8step":
    pm_file = "../data/ts_nju_pm_fit_10sigma-8step.dat"
    spin_file = "../logs/spin_fit_from_glo_ts-8step.txt"
elif ts_type == "nju_glo_10step":
    pm_file = "../data/ts_nju_pm_fit_3sigma-10step.dat"
    spin_file = "../logs/spin_fit_from_glo_ts-10step.txt"
elif ts_type == "nju_glo_20step":
    pm_file = "../data/ts_nju_pm_fit_10sigma-20step.dat"
    spin_file = "../logs/spin_fit_from_glo_ts-20step.txt"

# APM data
apm_tab = Table.read(pm_file, format="ascii.csv")

# Only ICRF3 defining sources
icrf3_def = Table.read("../data/icrf3sx-def-sou.txt", format="ascii")
apm_tab = join(icrf3_def, apm_tab, keys="iers_name")

# convert mas/yr into muas/yr
apm_tab["pmra"] = apm_tab["pmra"] * 1e3
apm_tab["pmra_err"] = apm_tab["pmra_err"] * 1e3
apm_tab["pmdec"] = apm_tab["pmdec"] * 1e3
apm_tab["pmdec_err"] = apm_tab["pmdec_err"] * 1e3

# Remove sources without apparant proper motion estimate.
mask = apm_tab["num_cln"] >= 5
apm_tab = apm_tab[mask]

# Output file
with open(spin_file, "w") as f_out:

    # Print header information
    print("num_sou, "
          "rx_mean, ry_mean, rz_mean, r_mean, gx_mean, gy_mean, gz_mean, g_mean, "
          "rx_std, ry_std, rz_std, r_std, gx_std, gy_std, gz_std, g_std, "
          "rx_q1, ry_q1, rz_q1, r_q1, gx_q1, gy_q1, gz_q1, g_q1, "
          "rx_q2, ry_q2, rz_q2, r_q2, gx_q2, gy_q2, gz_q2, g_q2, "
          "rx_q3, ry_q3, rz_q3, r_q3, gx_q3, gy_q3, gz_q3, g_q3", file=f_out)

    # Number of sources to be extracted from the whole sample
    num_sou = np.arange(10, 290, 10)

    # For each value of num_sou, do the bootstrap resampling and estimate the
    # global spin. Repeat this procedure 100 times.
    num_iter = 100

    for i, num_soui in enumerate(num_sou):

        print("Iteration : [{:4d}/{:4d}]".format(i+1, len(num_sou)))

        pmt = np.zeros((num_iter, 8))
        sig = np.zeros((num_iter, 8))

        for j in range(num_iter):
            # Resampling
            new_tab = bootstrap_sample(apm_tab, sam_size=num_soui, with_replace=False)

            # Do the LSQ fitting
            pmt[j, :], sig[j, :], output = vsh_fit_for_pm(new_tab)

        pmt_mean = np.mean(pmt, axis=0)
        pmt_std = np.std(pmt, axis=0)
        pmt_q1 = np.percentile(pmt, 16, axis=0)  # Lower limit for 1-sigma
        pmt_q2 = np.percentile(pmt, 50, axis=0)
        pmt_q3 = np.percentile(pmt, 84, axis=0)  # Upper limit for 1-sigma

        print("{:d},{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f},"
              "{:8.3f},{:8.3f},{:8.3f},{:8.3f}".format(
                  num_soui, *pmt_mean, *pmt_std, *pmt_q1, *pmt_q2, *pmt_q3), file=f_out)


# --------------------------------- END --------------------------------
