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
from my_progs.vlbi.sou_func import read_crf
from my_progs.catalog.read_icrf import read_icrf3, read_icrf2
from my_progs.catalog.read_gaia import read_edr3_icrf_sou
from my_progs.catalog.pos_diff import radio_cat_diff_calc
from tool_func import orient_angle_sampling, save_data

# -----------------------------  FUNCTIONS -----------------------------
np.random.seed(2)

icrf2 = read_icrf2()
icrf3sx = read_icrf3(wv="sx")
icrf3k = read_icrf3(wv="k")
icrf3xka = read_icrf3(wv="xka")
gedr3 = read_edr3_icrf_sou()

# Various solutions
opa2019a = read_crf("../data/ivs-solutions/opa2019a.crf",
                    drop_few_obs=True, analy_cen="opa")
asi2020a = read_crf("../data/ivs-solutions/asi2020a.crf",
                    drop_few_obs=True, analy_cen="asi")
aus2020b = read_crf("../data/ivs-solutions/aus2020b.crf",
                    drop_few_obs=True, analy_cen="aus")
usn2019c = read_crf("../data/ivs-solutions/usn2019c.crf",
                    drop_few_obs=True, analy_cen="usn")

oft_2_sx = radio_cat_diff_calc(icrf2, icrf3sx, sou_name="iers_name")
oft_k_sx = radio_cat_diff_calc(icrf3k, icrf3sx, sou_name="iers_name")
oft_ka_sx = radio_cat_diff_calc(icrf3xka, icrf3sx, sou_name="iers_name")
oft_g_sx = radio_cat_diff_calc(gedr3, icrf3sx, sou_name="iers_name")

oft_opa_sx = radio_cat_diff_calc(opa2019a, icrf3sx, sou_name="iers_name")
oft_asi_sx = radio_cat_diff_calc(asi2020a, icrf3sx, sou_name="iers_name")
oft_aus_sx = radio_cat_diff_calc(aus2020b, icrf3sx, sou_name="iers_name")
oft_usn_sx = radio_cat_diff_calc(usn2019c, icrf3sx, sou_name="iers_name")

# Do sampling of the orientation angle
# Sample size is around 66.7% of the common sources
# ------------------------------------
# 1) sampling without replacement

# ICRF2 vs. ICRF3 SX
# 3410 * 2 / 3 = 2273.3
print("ICRF2 - ICRF3 SX:")
opt2, opt21 = orient_angle_sampling(oft_2_sx, sam_size=2300)
save_data(opt2, "../logs/icrf2_icrf3sx_orient.fits")
save_data(opt21, "../logs/icrf2_icrf3sx_orient_cln.fits")

# ICRF3 K vs. ICRF3 SX
# 793 * 2 / 3 = 528.6
print("\n\nICRF3 K - ICRF3 SX:")
optk, optk1 = orient_angle_sampling(oft_k_sx, sam_size=550)
save_data(optk, "../logs/icrf3k_icrf3sx_orient.fits")
save_data(optk1, "../logs/icrf3k_icrf3sx_orient_cln.fits")

# ICRF3 XKa vs. ICRF3 SX
# 638 * 2 / 3 = 425.3
print("\n\nICRF3 XKa - ICRF3 SX:")
optka, optka1 = orient_angle_sampling(oft_ka_sx, sam_size=450)
save_data(optka, "../logs/icrf3xka_icrf3sx_orient.fits")
save_data(optka1, "../logs/icrf3xka_icrf3sx_orient_cln.fits")

# Gaia EDR3 vs. ICRF3 SX
# 3142 * 2 / 3 = 2095
print("\n\nGaia EDR3 - ICRF3 SX:")
optka, optka1 = orient_angle_sampling(oft_g_sx, sam_size=2100)
save_data(optka, "../logs/gedr3_icrf3sx_orient.fits")
save_data(optka1, "../logs/gedr3_icrf3sx_orient_cln.fits")

# opa2019a vs. ICRF3 SX
# 4380 * 2 / 3 = 2920
print("\n\nopa2019a - ICRF3 SX:")
opt, opt1 = orient_angle_sampling(oft_opa_sx, sam_size=3000)
save_data(opt, "../logs/opa2019a_icrf3sx_orient.fits")
save_data(opt1, "../logs/opa2019a_icrf3sx_orient_cln.fits")

# asi2020a vs. ICRF3 SX
# 4296 * 2 / 3 = 2864
print("\n\nasi2020a - ICRF3 SX:")
opt, opt1 = orient_angle_sampling(oft_asi_sx, sam_size=3000)
save_data(opt, "../logs/asi2020a_icrf3sx_orient.fits")
save_data(opt1, "../logs/asi2020a_icrf3sx_orient_cln.fits")

# aus2020b vs. ICRF3 SX
# 4456 * 2 / 3 = 2971
print("\n\n us2020b - ICRF3 SX:")
opt, opt1 = orient_angle_sampling(oft_aus_sx, sam_size=3000)
save_data(opt, "../logs/aus2020b_icrf3sx_orient.fits")
save_data(opt1, "../logs/aus2020b_icrf3sx_orient_cln.fits")

# usn2019c vs. ICRF3 SX
# 2263 * 2 / 3 = 1509
print("\n\nusn2019c - ICRF3 SX:")
opt, opt1 = orient_angle_sampling(oft_usn_sx, sam_size=1500)
save_data(opt, "../logs/usn2019c_icrf3sx_orient.fits")
save_data(opt1, "../logs/usn2019c_icrf3sx_orient_cln.fits")

# ------------------------------------
# 2) sampling with replacement
# ICRF2 vs. ICRF3 SX
# 3410 * 2 / 3 = 2273.3
# print("ICRF2 - ICRF3 SX:")
# opt2, opt21 = orient_angle_sampling(oft_2_sx, sam_size=2275, with_replace=True)
# save_data(opt2, "../logs/icrf2_icrf3sx_orient_norep.fits")
# save_data(opt21, "../logs/icrf2_icrf3sx_orient_cln_norep.fits")

# ICRF3 K vs. ICRF3 SX
# 793 * 2 / 3 = 528.6
# print("\n\nICRF3 K - ICRF3 SX:")
# optk, optk1 = orient_angle_sampling(oft_k_sx, sam_size=530, with_replace=True)
# save_data(optk, "../logs/icrf3k_icrf3sx_orient_norep.fits")
# save_data(optk1, "../logs/icrf3k_icrf3sx_orient_cln_norep.fits")

# ICRF3 XKa vs. ICRF3 SX
# 638 * 2 / 3 = 425.3
# print("\n\nICRF3 XKa - ICRF3 SX:")
# optka, optka1 = orient_angle_sampling(oft_ka_sx, sam_size=425, with_replace=True)
# save_data(optka, "../logs/icrf3xka_icrf3sx_orient_norep.fits")
# save_data(optka1, "../logs/icrf3xka_icrf3sx_orient_cln_norep.fits")

# opa2019a vs. ICRF3 SX
# 4380 * 2 / 3 = 2920
# print("\n\nopa2019a - ICRF3 SX:")
# opt, opt1 = orient_angle_sampling(oft_opa_sx, sam_size=2920, with_replace=True)
# save_data(opt, "../logs/opa2019a_icrf3sx_orient_norep.fits")
# save_data(opt1, "../logs/opa2019a_icrf3sx_orient_cln_norep.fits")
# --------------------------------- END --------------------------------
