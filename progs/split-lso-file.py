#!/usr/bin/env python
# coding: utf-8

import numpy as np

from my_progs.vlbi.read_lso import read_lso, write_ts

# read .lso file
opa2019r = read_lso("/Users/Neo/Astronomy/data/vlbi/opa/ts_sou/opa2019r.lso")

# Group by name
tabg = opa2019r.group_by("iers_name")
soulist = tabg.groups.keys

print("There are {:d} sources with coordinate time series available".format(len(soulist)))

# Write the data file
with open("../logs/split-lso-file.log", "w") as flog:
    for i, sou in enumerate(soulist["iers_name"]):
        tabsub = tabg.groups[i]
        tabsub.sort("epoch")
        write_ts(tabsub, datadir="../data/ts")

        print("Retrieve coordinate ts of {:s}: done!".format(sou), file=flog)

soulist.write("../data/sou_list.txt", format="ascii", overwrite=True)

