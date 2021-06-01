# Difference in RA and decl.
fig, (ax0, ax1) = plt.subplots(figsize=(10, 6), ncols=2, sharex=True, sharey=True)

# lim = 10
# x = np.arange(-lim, lim, 0.1)
# ax0.plot(x, x, "r--", lw=0.5)
# ax1.plot(x, x, "r--", lw=0.5)

ax0.plot(com_sou_00_01['used_obs_0'], com_sou_00_01['used_obs_1'], "b.", ms=2)
ax1.plot(com_sou_00_01['used_ses_0'], com_sou_00_01['used_ses_1'], "b.", ms=2)

# ax0.axis("square")
# ax1.axis("square")
# ax1.axis([-lim, lim, -lim, lim])

# ax0.set_xlabel('$\Delta\\alpha*$ (2003-c, mas)', fontsize=15)
# ax0.set_ylabel('$\Delta\\alpha*$ (2004-c, mas)', fontsize=15)
# ax1.set_xlabel('$\Delta\\delta$ (2003-c, mas)', fontsize=15)
ax1.set_ylabel('$\Delta\\delta$ (2004-c, mas)', fontsize=15)

fig.tight_layout()
