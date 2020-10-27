#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: icrf-orient-check.py
"""
Created on Thu Apr 23 17:49:29 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import join, Table
# from astropy.table import Column, Table
from astropy import units as u
import numpy as np

# My progs
from my_progs.vlbi.read_sou import read_crf
from my_progs.catalog.read_icrf import read_icrf3, read_icrf2
from my_progs.catalog.pos_diff import radio_cat_diff_calc
from tool_func import orient_angle_sampling, save_data

# -----------------------------  FUNCTIONS -----------------------------
icrf2 = read_icrf2()
icrf3sx = read_icrf3(wv="sx")
icrf3k = read_icrf3(wv="k")
icrf3xka = read_icrf3(wv="xka")
opa2019a = read_crf("../data/opa2019a.crf", drop_few_obs=True)

oft_2_sx = radio_cat_diff_calc(icrf2, icrf3sx, sou_name="iers_name")

oft_k_sx = radio_cat_diff_calc(icrf3k, icrf3sx, sou_name="iers_name")

oft_ka_sx = radio_cat_diff_calc(icrf3xka, icrf3sx, sou_name="iers_name")

oft_opa_sx = radio_cat_diff_calc(opa2019a, icrf3sx, sou_name="iers_name")

# Do sampling of the orientation angle
# Sample size is around 66.7% of the common sources
# ------------------------------------
# 1) sampling without replacement

# ICRF2 vs. ICRF3 SX
# 3410 * 2 / 3 = 2273.3
# print("ICRF2 - ICRF3 SX:")
# opt2, opt21 = orient_angle_sampling(oft_2_sx, sam_size=2275)
# save_data(opt2, "../logs/icrf2_icrf3sx_orient.fits")
# save_data(opt21, "../logs/icrf2_icrf3sx_orient_cln.fits")

# ICRF3 K vs. ICRF3 SX
# 793 * 2 / 3 = 528.6
# print("\n\nICRF3 K - ICRF3 SX:")
# optk, optk1 = orient_angle_sampling(oft_k_sx, sam_size=530)
# save_data(optk, "../logs/icrf3k_icrf3sx_orient.fits")
# save_data(optk1, "../logs/icrf3k_icrf3sx_orient_cln.fits")

# ICRF3 XKa vs. ICRF3 SX
# 638 * 2 / 3 = 425.3
# print("\n\nICRF3 XKa - ICRF3 SX:")
# optka, optka1 = orient_angle_sampling(oft_ka_sx, sam_size=425)
# save_data(optka, "../logs/icrf3xka_icrf3sx_orient.fits")
# save_data(optka1, "../logs/icrf3xka_icrf3sx_orient_cln.fits")

# opa2019a vs. ICRF3 SX
# 4380 * 2 / 3 = 2920
# print("\n\nopa2019a - ICRF3 SX:")
# opt, opt1 = orient_angle_sampling(oft_opa_sx, sam_size=2920)
# save_data(opt, "../logs/opa2019a_icrf3sx_orient.fits")
# save_data(opt1, "../logs/opa2019a_icrf3sx_orient_cln.fits")

# ------------------------------------
# 2) sampling with replacement

# ICRF2 vs. ICRF3 SX
# 3410 * 2 / 3 = 2273.3
print("ICRF2 - ICRF3 SX:")
opt2, opt21 = orient_angle_sampling(oft_2_sx, sam_size=2275, with_replace=True)
save_data(opt2, "../logs/icrf2_icrf3sx_orient_norep.fits")
save_data(opt21, "../logs/icrf2_icrf3sx_orient_cln_norep.fits")

# ICRF3 K vs. ICRF3 SX
# 793 * 2 / 3 = 528.6
print("\n\nICRF3 K - ICRF3 SX:")
optk, optk1 = orient_angle_sampling(oft_k_sx, sam_size=530, with_replace=True)
save_data(optk, "../logs/icrf3k_icrf3sx_orient_norep.fits")
save_data(optk1, "../logs/icrf3k_icrf3sx_orient_cln_norep.fits")

# ICRF3 XKa vs. ICRF3 SX
# 638 * 2 / 3 = 425.3
print("\n\nICRF3 XKa - ICRF3 SX:")
optka, optka1 = orient_angle_sampling(oft_ka_sx, sam_size=425, with_replace=True)
save_data(optka, "../logs/icrf3xka_icrf3sx_orient_norep.fits")
save_data(optka1, "../logs/icrf3xka_icrf3sx_orient_cln_norep.fits")

# opa2019a vs. ICRF3 SX
# 4380 * 2 / 3 = 2920
print("\n\nopa2019a - ICRF3 SX:")
opt, opt1 = orient_angle_sampling(oft_opa_sx, sam_size=2920, with_replace=True)
save_data(opt, "../logs/opa2019a_icrf3sx_orient_norep.fits")
save_data(opt1, "../logs/opa2019a_icrf3sx_orient_cln_norep.fits")
# --------------------------------- END --------------------------------
