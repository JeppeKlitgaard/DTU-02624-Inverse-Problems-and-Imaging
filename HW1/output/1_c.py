np.random.seed(42)


def x_func(s):
    return s * np.sin(s) * np.exp((s - np.pi / 2) ** 2)


x = x_func(s)
y = A @ x
unscaled_noise = np.random.randn(*y.shape)
unscaled_noise_level = np.linalg.norm(unscaled_noise)
noise_coeff = 0.02 * np.linalg.norm(y) / unscaled_noise_level
noise = noise_coeff * unscaled_noise
y_noisy = y + noise
x_naive = np.linalg.inv(A) @ y_noisy

# Prepare for Picard plot
singular_values = np.diagonal(S)
expansion_coefficients = np.abs(U.T @ y_noisy)
reconstruction_coefficients = expansion_coefficients / singular_values
plt.figure(figsize=(PAPER_WIDTH_IN, 2.3))
plt.plot(j, singular_values, marker="o", ls="-", ms=5, label="Singular, $μ_j$")
plt.plot(
    j,
    expansion_coefficients,
    marker="s",
    ls="-",
    ms=5,
    label="Expansion, $|(y^δ, y_j)|$",
)
plt.plot(
    j,
    reconstruction_coefficients,
    marker="^",
    ls="-",
    ms=5,
    label="Reconstruction, $|(y^δ, y_j)|/μ_j$",
)
plt.semilogy()
plt.legend()
plt.xlim(0, 25)
plt.ylim(1e-6, 1e4)
plt.axvline(9, color="k", ls="--")
plt.xticks(list(plt.xticks()[0]) + [9])
plt.title("Picard Plot with Noisy Data ($2\\%$ rel. noise)")
plt.xlabel("Index, $j$")
plt.ylabel("Magnitude")
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "1_c_picard_plot.png", dpi=600, bbox_inches="tight")
display(plt.gcf())

plt.gcf().set_size_inches(PAPER_WIDTH_IN, 4)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "1_c_picard_plot_large.png", dpi=600, bbox_inches="tight")

# Compute reconstruction
j_cutoff = 9  # Inclusive
S_truncated = S[:j_cutoff, :j_cutoff]
U_truncated = U[:, :j_cutoff]
V_truncated = V.T[:, :j_cutoff]
K_inv_truncated = V_truncated @ np.linalg.inv(S_truncated) @ U_truncated.T
x_reconstructed = K_inv_truncated @ y_noisy

plt.figure(figsize=(PAPER_WIDTH_IN, 2.0))
plt.plot(s, x_func(s), label="$x(s)$", ls="--", lw=1.5)
plt.plot(s, x_reconstructed, label="$x_\\text{rec}(s)$", lw=1.5)
plt.plot(s, y_noisy, label="$y_\\text{noisy}(s)$", lw=1.5)
plt.plot(s, y, label="$y(s)$", ls=":", lw=2.5)
plt.legend()

# Create new axis for naive solution
ax_main = plt.gca()
ax_noisy = plt.twinx()
ax_noisy.plot(s, x_naive, label="$x_\\text{naive}(s)$", color="C3", ls=":", lw=1)
ax_noisy.grid(False)
ax_main.set_zorder(1)
ax_main.patch.set_visible(False)
ax_noisy.set_ylabel("$x_\\text{naive}(s)$")

# Merge legends
handles1, labels1 = ax_main.get_legend_handles_labels()
handles2, labels2 = ax_noisy.get_legend_handles_labels()
ax_main.legend(handles1 + handles2, labels1 + labels2, loc="upper left")

plt.xlabel("Spatial Coordinate, $s$")
plt.title("Noisy Backward Heat Equation and Reconstruction")
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "1_c_reconstruction.png", dpi=600, bbox_inches="tight")
display(plt.gcf())

plt.gcf().set_size_inches(PAPER_WIDTH_IN, 3.5)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "1_c_reconstruction_large.png", dpi=600, bbox_inches="tight")
