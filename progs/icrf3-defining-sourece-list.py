#!/usr/bin/env python
# coding: utf-8

# This notebook is used to generate the defining source lists of ICRF3 at S/X-, K-, and X/Ka-band.

from my_progs.catalog.read_icrf import read_icrf3


# S/X band
sx = read_icrf3(wv="sx")

mask = sx["type"] == "D"
sxdef = sx[mask]

sxdef.write("../data/icrf3sx-def-sou.txt", format="ascii",
           include_names=["iers_name"], overwrite=True)

# K band
k = read_icrf3(wv="k")

mask = k["type"] == "D"
kdef = k[mask]

kdef.write("../data/icrf3k-def-sou.txt", format="ascii",
           include_names=["iers_name"], overwrite=True)


# X/Ka band
xka = read_icrf3(wv="xka")

mask = xka["type"] == "D"
xkadef = xka[mask]

xkadef.write("../data/icrf3xka-def-sou.txt", format="ascii",
           include_names=["iers_name"], overwrite=True)

