

def heaviside_step(x):
    return 0.5 * (np.sign(x) + 1)


def x_func_discontinuous(s):
    return x_func(s) + 4 * heaviside_step(s - 1 / 2)


x_true = x_func_discontinuous(s)
y_true = A @ x_true

# Add noise
np.random.seed(1234)
unscaled_noise = np.random.randn(*y_true.shape)
unscaled_noise_level = np.linalg.norm(unscaled_noise)
noise_coeff = 0.05 * np.linalg.norm(y_true) / unscaled_noise_level
noise = noise_coeff * unscaled_noise
y_noisy = y_true + noise

# Best TSVD
k_best = 50
U, S, Vh = np.linalg.svd(A)
U_truncated = U[:, :k_best]
S_truncated = np.diag(S)[:k_best, :k_best]
V_truncated = Vh.T[:, :k_best]
A_inv_truncated = V_truncated @ np.linalg.inv(S_truncated) @ U_truncated.T
x_best_tsvd = A_inv_truncated @ y_noisy

# Slow TSVD
k_slow = 20
U, S, Vh = np.linalg.svd(A)
U_truncated = U[:, :k_slow]
S_truncated = np.diag(S)[:k_slow, :k_slow]
V_truncated = Vh.T[:, :k_slow]
A_inv_truncated = V_truncated @ np.linalg.inv(S_truncated) @ U_truncated.T
x_slow_tsvd = A_inv_truncated @ y_noisy

# Best Tikhonov
alpha = 1e-2
x_best_tikhonov = np.linalg.inv(A.T @ A + alpha * np.eye(n)) @ A.T @ y_noisy

fig, (ax_main, ax_zoom) = plt.subplots(
    1, 2, figsize=(PAPER_WIDTH_IN, 2.3), gridspec_kw={"width_ratios": [3, 1]}
)
ax_main.plot(s, x_best_tikhonov, label=f"$x_{{α={alpha}}}(s)$", lw=1)
ax_zoom.plot(s, x_best_tikhonov, label=f"$x_{{α={alpha}}}(s)$", lw=1)
ax_main.plot(s, x_best_tsvd, label=f"$x_{{k={k_best}}}(s)$", lw=1)
ax_zoom.plot(s, x_best_tsvd, label=f"$x_{{k={k_best}}}(s)$", lw=1)
ax_main.plot(s, x_slow_tsvd, label=f"$x_{{k={k_slow}}}(s)$", lw=1)
ax_zoom.plot(s, x_slow_tsvd, label=f"$x_{{k={k_slow}}}(s)$", lw=1)
ax_main.plot(s, x_true, label="$x(s)$", ls="--", color="k")
ax_zoom.plot(s, x_true, label="$x(s)$", ls="--", color="k")
ax_main.legend()
fig.suptitle(
    "TSVD and Tikhonov Regularisation of Discontinuous Function", y=0.90, fontsize=10
)
ax_main.set_xlabel("Spatial Coordinate, $s$")
ax_zoom.set_xlabel("Spatial Coordinate, $s$")
ax_zoom.set_xlim(0.45, 0.55)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "2_d_discontinuous.png", dpi=600, bbox_inches="tight")
display(fig)

fig.set_size_inches(PAPER_WIDTH_IN, 3.5)
fig.suptitle(fig._suptitle.get_text(), y=0.95, fontsize=12)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "2_d_discontinuous_large.png", dpi=600, bbox_inches="tight")
