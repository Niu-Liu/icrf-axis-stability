#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: generate-sou-list.py
"""
Created on Wed 02 Jun 2021 07:21:35 PM CST

@author: Neo(niu.liu@nju.edu.cn)
"""

import os
import sys
import numpy as np

from astropy.table import join, setdiff, Table, vstack


# -----------------------------  FUNCTIONS -----------------------------
def read_icrf3(icrf3_file="/data/catalogs/icrf/icrf3sx.txt"):
    """Read the ICRF3 catalog.

    Parameter
    ---------
    icrf3_file : string
        file name and path of the ICRF3 catalog

    Return
    ------
    icrf3 : an astropy.Table object
        data in the catalog
    """

    icrf3 = Table.read(icrf3_file,
                       format="ascii.fixed_width", data_start=16,
                       names=["icrf_name", "iers_name", "type"],
                       col_starts=[0, 25, 35],
                       col_ends=[20, 32, 35])

    return icrf3


def read_sou(sou_file):
    """Read radio source positions

    Parameters
    ----------
    sou_file : string
        the full path of .sou file

    Return
    ------
    t_sou : astropy.table object
        |
        -- sou_name : str
            IVS or IERS source name
    """

    if not os.path.isfile(sou_file):
        sys.exit()

    t_sou = Table.read(sou_file, format="ascii.fixed_width_no_header",
                       names=("sou_name"),
                       col_starts=(10),
                       col_ends=(18))
    return t_sou


def write_nnrs(soun, list_file):
    """Write the No_Net_Rotation constraint source list.

    Parameters
    ----------
    soun: array, string
        IVS/IERS designation name of radio source
    list_file: Output file
    """

    N = 8  # 8 source pre line

    fopt = open(list_file, "w")

    for i, souni in enumerate(soun):
        if not i % N:
            if i:
                print("\\\n     ", end="", file=fopt)
            else:
                print("     ", end="", file=fopt)

        print("%-9s" % souni, end="", file=fopt)

    fopt.close()


# -----------------------------  MAINS -----------------------------
if len(sys.argv) != 2:
    print("Usage: {} <num_step>  <sou_file>".format(sys.argv[0]))
    print("Example: {} 20 opa2021a.sou".format(sys.argv[0]))
    print("num_step : positive integer.")
    print("sou_file : .sou file outputed from Calc/Solve.")
    print("Note you may wish to edit the line18 for the path to the icrf3sx.txt file,")
    print("which can be downloaded at https://hpiers.obspm.fr/icrs-pc/newwww/icrf/icrf3sx.txt.")

# Number of steps to generate the global time series
nb_step = int(sys.argv[1])
sou_file = sys.argv[2]

# ICRF3 S/X catalog
icrf3sx = read_icrf3()

# Create a sub-table for the ICRF3 defining sources
icrf3_def = Table([icrf3sx[icrf3sx["type"] == "D"]["sou_name"]])

# Keep only source name
icrf3sx.keep_columns(["sou_name"])

# Radio source catalog from the global solution
glo_sou = read_sou(sou_file)
glo_sou.keep_columns(["sou_name"])

# Cross-match with the ICRF3 S/X-band catalog
# ICRF3 defining sources
def_sou = join(glo_sou, icrf3_def)
oth_sou = setdiff(glo_sou, def_sou)

print("There are {} sources in the global solutions, "
      "{} being the ICRF3 defining sources while {} being others.".format(
          len(glo_sou), len(def_sou), len(oth_sou)))

# Divide the ICRF3 defining and non-defining sources into N subsets
print("I will divide sources into {} subsets.".format(nb_step))

for i in range(nb_step):
    def0 = def_sou[i::nb_step]
    oth0 = oth_sou[nb_step-i::nb_step]

    # Combine the defning and non-defning source list
    cat0 = vstack((def0, oth0))

    print("Writing the source list #{}...".format(i+1))
    file0 = "list{}".format(i)
    print("output file:", file0)
    write_nnrs(cat0["sou_name"], file0)
    print("The number of sources in the subset is {}.".format(len(cat0)))

print("Done!")
# --------------------------------- END --------------------------------
