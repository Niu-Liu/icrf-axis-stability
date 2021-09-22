fig, ax = plt.subplots()

ax.errorbar(output1_4["G_ra"], output1_4["G_dec"],
            xerr=output1_4["G_ra_err"],
            yerr=output1_4["G_dec_err"],
            capsize=3, elinewidth=0.5,
            fmt="v", color="r", label=" 4-step")
ax.errorbar(output1_8["G_ra"], output1_8["G_dec"],
            xerr=output1_8["G_ra_err"],
            yerr=output1_8["G_dec_err"],
            capsize=3, elinewidth=0.5,
            fmt="x", color="g", label=" 8-step")
ax.errorbar(output1_10["G_ra"], output1_10["G_dec"],
            xerr=output1_10["G_ra_err"],
            yerr=output1_10["G_dec_err"],
            capsize=3, elinewidth=0.5,
            fmt="^", color="b", label="10-step")
ax.errorbar(output1_20["G_ra"], output1_20["G_dec"],
            xerr=output1_20["G_ra_err"],
            yerr=output1_20["G_dec_err"],
            capsize=3, elinewidth=0.5,
            fmt="s", color="m", label="20-step")


ax.plot(265, -30, "ko", label="GC")

ax.legend()

ax.set_title("Direction of glide axis", fontsize=15)
ax.set_xlabel("RA (degree)", fontsize=15)
ax.set_ylabel("Decl. (degree)", fontsize=15)

plt.tight_layout()
