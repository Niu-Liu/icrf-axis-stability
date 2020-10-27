#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: plot-long-ts.py
"""
Created on Fri Mar 27 16:52:20 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from astropy.table import Table
import numpy as np
#
from plot_ts import plot_ts


# -----------------------------  FUNCTIONS -----------------------------
sou = Table.read("../data/sou_list_longhist.txt", format="ascii")

for soui in sou["iers_name"]:
    plot_ts(soui)
    print("Plot the ts for {}: done!".format(soui))
# --------------------------------- END --------------------------------
