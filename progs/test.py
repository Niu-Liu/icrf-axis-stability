# Various solutions
opa2019a = read_crf("../data/ivs-solutions/opa2019a.crf",
                    drop_few_obs=True, analy_cen="opa")
asi2020a = read_crf("../data/ivs-solutions/asi2020a.crf",
                    drop_few_obs=True, analy_cen="asi")
aus2020b = read_crf("../data/ivs-solutions/aus2020b.crf",
                    drop_few_obs=True, analy_cen="aus")
usn2019c = read_crf("../data/ivs-solutions/usn2019c.crf",
                    drop_few_obs=True, analy_cen="usn")

oft_opa_sx = radio_cat_diff_calc(opa2019a, icrf3sx, sou_name="iers_name")
oft_opa_sx_def = join(oft_opa_sx, def_list, keys="iers_name")

oft_asi_sx = radio_cat_diff_calc(asi2020a, icrf3sx, sou_name="iers_name")
oft_asi_sx_def = join(oft_asi_sx, def_list, keys="iers_name")

oft_aus_sx = radio_cat_diff_calc(aus2020b, icrf3sx, sou_name="iers_name")
oft_aus_sx_def = join(oft_aus_sx, def_list, keys="iers_name")

oft_usn_sx = radio_cat_diff_calc(usn2019c, icrf3sx, sou_name="iers_name")
oft_usn_sx_def = join(oft_usn_sx, def_list, keys="iers_name")

    # asi2020a
    [
        np.mean(optasi1["R1"]),
        np.std(optasi1["R1"]),
        np.mean(optasi1["R2"]),
        np.std(optasi1["R2"]),
        np.mean(optasi1["R3"]),
        np.std(optasi1["R3"])
    ],
    # aus2020b
    [
        np.mean(optaus1["R1"]),
        np.std(optaus1["R1"]),
        np.mean(optaus1["R2"]),
        np.std(optaus1["R2"]),
        np.mean(optaus1["R3"]),
        np.std(optaus1["R3"])
    ],
    # usn2019c
    [
        np.mean(optusn1["R1"]),
        np.std(optusn1["R1"]),
        np.mean(optusn1["R2"]),
        np.std(optusn1["R2"]),
        np.mean(optusn1["R3"]),
        np.std(optusn1["R3"])
