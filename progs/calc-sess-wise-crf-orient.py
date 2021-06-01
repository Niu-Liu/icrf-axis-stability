#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: calc-sess-wise-crf-orient.py
"""
Created on Sat May  8 20:34:45 2021

@author: Neo(niu.liu@nju.edu.cn)
"""

import numpy as np

from astropy.table import Table, join

from my_progs.vlbi.ts_func import get_ts
from my_progs.catalog.read_icrf import read_icrf3
from my_progs.catalog.pos_diff import radio_cat_diff_calc
from tool_func import sample_clean, vsh_fit_for_pos

# ICRF3 catalog
icrf3sx = read_icrf3(wv="sx")
icrf3def = icrf3sx[icrf3sx["type"] == "D"]
icrf3def = Table([icrf3def["iers_name"]])

# Source coordinate time series in one file
sou_ts = get_ts("sou-ti-concat", "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou")
# sou_ts = get_ts("sou-ti-concat", "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou-from-glo")

# Group by the session name
sou_ts_g = sou_ts.group_by("db_name")
num_sess = len(sou_ts_g.groups)

f_out_all = open("../data/sess-wise-orient-all.txt", "w")
header = ["sess_name, mjd, num_sou, "
          "gx, gy, gz, rx, ry, rz, g, r, "
          "gx_err, gy_err, gz_err, rx_err, ry_err, rz_err, g_err, r_err, "
          "jyear"][0]
print(header, file=f_out_all)
f_out_def = open("../data/sess-wise-orient-def.txt", "w")
print(header, file=f_out_def)

# Line format
lfmt = "{:9s},{:9.3f},{:4d}," + "{:6.0f}," * 16 + "{:f}"

print("There are {:5d} sessions in total.".format(num_sess))
i = 0


for group in sou_ts_g.groups:

    sess_name = group["db_name"][0]
    mjd = group["mjd"][0]
    jyear = group["jyear"][0]

    print("Processing session {:10s} [{:4d}/{:4d}]".format(sess_name, i+1, num_sess))
    i += 1

    # Number of source should be greater than 5
    if len(group) < 6:
        continue

    # Offset wrt
    pos_oft = radio_cat_diff_calc(group, icrf3sx, sou_name="iers_name")
    pos_oft_cln = sample_clean(pos_oft, rho0=10)

    if len(pos_oft_cln) >= 6:
        N_all = len(pos_oft_cln)
        pmt, sig = vsh_fit_for_pos(pos_oft_cln)

        print(lfmt.format(sess_name, mjd, N_all, *pmt, *sig, jyear),
              file=f_out_all)

    pos_oft_def = join(pos_oft, icrf3def, keys="iers_name")

    if len(pos_oft_def) < 6:
        continue

    pos_oft_def_cln = sample_clean(pos_oft_def, rho0=10)

    if len(pos_oft_def_cln) >= 6:
        N_def = len(pos_oft_def_cln)
        pmt, sig = vsh_fit_for_pos(pos_oft_def_cln)

        print(lfmt.format(sess_name, mjd, N_def, *pmt, *sig, jyear),
              file=f_out_def)

# # Save the data
f_out_all.close()
f_out_def.close()
