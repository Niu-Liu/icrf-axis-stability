#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: plot-long-ts.py
"""
Created on Fri Mar 27 16:52:20 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from astropy.table import Table
#
from plot_ts import plot_ts


# -----------------------------  FUNCTIONS -----------------------------
# sou = Table.read("../data/sou_list_longhist.txt", format="ascii")
# sou = Table.read("../data/ts-sou.list", format="ascii")
#
# for soui in sou["iers_name"]:
#     plot_ts(soui, "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou",
#             "/Users/Neo/Astronomy/data/vlbi/opa/ts-plot")
#     print("Plot the ts for {}: done!".format(soui))

sou = Table.read("../data/sou-ts-glo.list", format="ascii")
for soui in sou["iers_name"]:
    plot_ts(soui, "/Users/Neo/Astronomy/data/vlbi/opa/ts-sou-from-glo",
            "/Users/Neo/Astronomy/data/vlbi/opa/ts-plot-from-glo")
    print("Plot the ts for {}: done!".format(soui))
# --------------------------------- END --------------------------------
