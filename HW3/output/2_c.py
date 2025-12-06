fig, axs = plt.subplots(2, 2, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN))

axs[0, 0].imshow(u_sol_tv, cmap="gray", norm=gray_norm)
axs[0, 0].set_title("TV Inpainting")

axs[0, 1].imshow(u_sol_tv_noiseless, cmap="gray", norm=gray_norm)
axs[0, 1].set_title("TV Inpainting (Noiseless)")

axs[1, 0].imshow(u_sol_harmonic, cmap="gray", norm=gray_norm)
axs[1, 0].set_title("Harmonic Inpainting")

axs[1, 1].imshow(u_sol_harmonic_noiseless, cmap="gray", norm=gray_norm)
axs[1, 1].set_title("Harmonic Inpainting (Noiseless)")

for ax in axs.flatten():
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "2_c.png", dpi=EXPORT_DPI, bbox_inches='tight')
