# Compute naive solution
x_naive = np.linalg.inv(A) @ y_noisy

# Perform SVD
j = np.arange(1, n + 1)
U, S, Vh = np.linalg.svd(A)

# Prepare Picard
expansion_coefficients = np.abs(U.T @ y_noisy)
reconstruction_coefficients = expansion_coefficients / S
# Perform Picard plot
plt.figure(figsize=(PAPER_WIDTH_IN, 2.3))
plt.title("Picard Plot for Abel Problem with $x(s) = s^2$")
plt.plot(j, S, marker="o", ls="-", ms=5, label="Singular,\n$μ_k$")
plt.plot(
    j,
    expansion_coefficients,
    marker="s",
    ls="-",
    ms=5,
    label="Expansion,\n$|(y^δ, y_k)|$",
)
plt.plot(
    j,
    reconstruction_coefficients,
    marker="^",
    ls="-",
    ms=5,
    label="Reconstruction,\n$|(y^δ, y_k)|/μ_k$",
)
plt.xlim(0, 50)
plt.ylim(1e-3, 1e1)
plt.xlabel("Index, $k$")
plt.ylabel("Magnitude")
plt.semilogy()
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "2_b_picard.png", dpi=600, bbox_inches="tight")
display(plt.gcf())

plt.gcf().set_size_inches(PAPER_WIDTH_IN, 4)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "2_b_picard_large.png", dpi=600, bbox_inches="tight")

# Reconstruction plots
fig, (ax_real, ax_error) = plt.subplots(
    1, 2, figsize=(PAPER_WIDTH_IN, 2.0), sharex=True
)

# ks = [5, 10, 25, 50, 100, "naive"]
ks = [5, 10, 20, 50, 100]
plot_kws_config = {
    50: {"ls": "--"},
    100: {"ls": ":"},
}
for k in ks:
    if k == "naive":
        x_reconstructed = x_naive
        show_error = False
    else:
        # To TSVD
        U_truncated = U[:, :k]
        S_truncated = np.diag(S)[:k, :k]
        V_truncated = Vh.T[:, :k]
        A_inv_truncated = V_truncated @ np.linalg.inv(S_truncated) @ U_truncated.T
        x_reconstructed = A_inv_truncated @ y_noisy
        show_error = True

    plot_kws = plot_kws_config.get(k, {})
    plot_kws.setdefault("lw", 1)
    plot_kws.setdefault("ls", "-")
    real_plot = ax_real.plot(s, x_reconstructed, label=f"$x_{{{k}}}(s)$", **plot_kws)

    if show_error:
        residual = x_true - x_reconstructed
        relative_error = (x_true - x_reconstructed) / np.linalg.norm(x_true)
        ax_error.plot(s, relative_error, color=real_plot[0].get_color(), **plot_kws)

ax_real.plot(s, x_true, label="$x(s)$", lw=1, ls="--", color="k")
ax_real.set_xlabel("Spatial Coordinate, $s$")
ax_error.set_xlabel("Spatial Coordinate, $s$")
ax_real.set_title("Truncated SVD Reconstructions")
ax_error.set_title("Relative Reconstruction Errors")

# Legend stuff
# Collect all handles and labels from both axes
handles_real, labels_real = ax_real.get_legend_handles_labels()
handles_error, labels_error = ax_error.get_legend_handles_labels()

# Place combined legend on ax_error, to the right
ax_error.legend(
    handles_real + handles_error,
    labels_real + labels_error,
    loc="upper left",
    bbox_to_anchor=(1, 1),
)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "2_b_reconstructions.png", dpi=600, bbox_inches="tight")
display(fig)

fig.set_size_inches(PAPER_WIDTH_IN, 3)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "2_b_reconstructions_large.png", dpi=600, bbox_inches="tight")
