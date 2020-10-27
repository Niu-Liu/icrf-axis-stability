#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: download_ts.py
"""
Created on Fri Jan 17 12:38:50 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Download the coordinate time series from OPA
"""

from astropy.table import Table
import astropy.units as u
import numpy as np
import os
import sys
from urllib.parse import quote

# Link to OPA data center
dataurl = "http://ivsopar.obspm.fr/radiosources/series/"


# -----------------------------  FUNCTIONS -----------------------------
def find_iers_name():
    """Find the

    Parameters
    ----------

    Returns
    -------
    """


def download_ts(souname, outdir="../data/"):
    """Download coordinate time series from OPA.

    Parameters
    ----------
    souname: string
        source name of IERS designation
    outdir: string
        directory to store the data file

    Returns
    -------
    None
    """

    ##### TO-DO  #####
    if souname not in ierstab:
        souname, status = find_iers_name(souname)

    if not status:
        print("Could not find the IERS name for", souname)
        sys.exit()

    ##### TO-DO  #####

    ##### TO-DO  #####
    # Save log somewhere
    # default name by download moment
    ##### TO-DO  #####

    ##### TO-DO  #####
    # Check download log for more information
    ##### TO-DO  #####

    os.system("wget -c -N -q -P {} {}{}.txt".format(
        outdir, dataurl, quote(souname)))

    print("Download coordinate time series "
          "for source {} successfully.".format(souname))


if __name__ == "__main__":
    download_ts(sys.argv[1])
# --------------------------------- END --------------------------------
