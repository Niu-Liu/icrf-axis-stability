fig, (ax0, ax1, ax2, ax3) = plt.subplots(figsize=(8, 8),
                                         nrows=4,
                                         sharex=True,
                                         sharey=True)

color = "grey"

ax0.hist(apm_tab4["pmdec"],
         bins=bin_array,
         color=color,
         fill=False,
         label="All")
ax0.hist(apm_def4["pmdec"], bins=bin_array, color=color, label="Defining")

ax1.hist(apm_tab8["pmdec"],
         bins=bin_array,
         color=color,
         fill=False,
         label="All")
ax1.hist(apm_def8["pmdec"], bins=bin_array, color=color, label="Defining")

ax2.hist(apm_tab10["pmdec"],
         bins=bin_array,
         color=color,
         fill=False,
         label="All")
ax2.hist(apm_def10["pmdec"], bins=bin_array, color=color, label="Defining")

ax3.hist(apm_tab20["pmdec"],
         bins=bin_array,
         color=color,
         fill=False,
         label="All")
ax3.hist(apm_def20["pmdec"], bins=bin_array, color=color, label="Defining")

ax3.set_xlabel("$\\mu_{\\alpha*}$ (mas/yr)", fontsize=15)

ax0.set_ylabel(" 4-Step", color=color, fontsize=15)
ax1.set_ylabel(" 8-Step", color=color, fontsize=15)
ax2.set_ylabel("10-Step", color=color, fontsize=15)
ax3.set_ylabel("20-Step", color=color, fontsize=15)

ax3.axis([-100, 100, 0, 400])
ax0.legend(loc="upper left")

plt.tight_layout()
