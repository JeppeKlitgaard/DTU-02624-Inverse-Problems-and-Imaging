alphas = np.geomspace(1e-4, 10, 10)
alphas = [1e-4, 1e-3, 1e-2, 1e-1, 1, 10]

plt.figure(figsize=(PAPER_WIDTH_IN, 2.2))
for alpha in alphas:
    x_tikhonov = np.linalg.inv(A.T @ A + alpha * np.eye(n)) @ A.T @ y_noisy
    plt.plot(s, x_tikhonov, label=f"$x_{{{alpha}}}(s)$")

# plt.plot(s, x_naive, label='$x_{\\mathrm{naive}}(s)$', ls=":")
plt.plot(s, x_true, label="$x(s)$", ls="--", color="k")
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("Tikhonov Regularised Reconstructions")
plt.xlabel("Spatial Coordinate, $s$")
plt.ylim(-1, 1.5)
plt.tight_layout()
plt.savefig(
    OUTPUT_DIR / "2_c_tikhonov_reconstructions.png", dpi=600, bbox_inches="tight"
)
display(plt.gcf())

plt.gcf().set_size_inches(PAPER_WIDTH_IN, 4)
plt.tight_layout()
plt.savefig(
    OUTPUT_DIR / "2_c_tikhonov_reconstructions_large.png", dpi=600, bbox_inches="tight"
)
# Best TSVD
k = 20
U, S, Vh = np.linalg.svd(A)
U_truncated = U[:, :k]
S_truncated = np.diag(S)[:k, :k]
V_truncated = Vh.T[:, :k]
A_inv_truncated = V_truncated @ np.linalg.inv(S_truncated) @ U_truncated.T
x_best_tsvd = A_inv_truncated @ y_noisy

# Best Tikhonov
alpha = 4e-2
x_best_tikhonov = np.linalg.inv(A.T @ A + alpha * np.eye(n)) @ A.T @ y_noisy

plt.figure(figsize=(PAPER_WIDTH_IN, 2.1))
plt.plot(s, x_best_tikhonov, label=f"$x_{{α={alpha}}}(s)$", lw=1)
plt.plot(s, x_best_tsvd, label=f"$x_{{k={k}}}(s)$", lw=1)
plt.plot(s, x_true, label="$x(s)$", ls="--", color="k")
plt.legend()
plt.title("TSVD and Tikhonov Comparison")
plt.xlabel("Spatial Coordinate, $s$")
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "2_c_tsvd_vs_tikhonov.png", dpi=600, bbox_inches="tight")
display(plt.gcf())

plt.gcf().set_size_inches(PAPER_WIDTH_IN, 4)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "2_c_tsvd_vs_tikhonov_large.png", dpi=600, bbox_inches="tight")
