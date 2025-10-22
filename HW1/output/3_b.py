# Deblurring problem
n = 2**8
h = 12 / n
n4 = n // 4
delta = 1e-7
alpha = delta**3

r1 = np.zeros(n)

# Cosine values for the formula
c = np.cos(np.arange(-1, n4 + 1) * 4 * np.pi / n)

# Fill in the first n/4 entries
r1[:n4] = h + 9 / (h * np.pi**2) * (2 * c[1 : n4 + 1] - c[:n4] - c[2 : n4 + 2])

# Middle entry
r1[n4] = h / 2 + 9 / (h * np.pi**2) * (np.cos(4 * np.pi / n) - 1)

A = toeplitz(r1)


# Set up true image
def tophat_func(s, width, center):
    return heaviside_step(s - (center - width / 2)) - heaviside_step(
        s - (center + width / 2)
    )


def x_deblur_func(s):
    return np.exp(-((s + 2) ** 2) / 0.5) + tophat_func(s, 2, 2)


s = np.linspace(-6, 6, n)

x_true = x_deblur_func(s)
y_true = A @ x_true

# Add noise
np.random.seed(42)
unscaled_noise = np.random.randn(*y_true.shape)
unscaled_noise_level = np.linalg.norm(unscaled_noise)
noise_coeff = delta * np.linalg.norm(y_true) / unscaled_noise_level
noise = noise_coeff * unscaled_noise
y_noisy = y_true + noise

# Compute regularisation
U, S, Vh = np.linalg.svd(A)

Y_basis = U
X_basis = Vh.T
mu = S
y_yj_ip = Y_basis.T @ y_noisy
filter_factors = mu**2 / (mu**3 + alpha)
coefficients = filter_factors * y_yj_ip
x_deblurred = X_basis @ coefficients

# Naive solution
x_true_naive = np.linalg.inv(A) @ y_true
x_naive = np.linalg.inv(A) @ y_noisy

fig, ax_x = plt.subplots(figsize=(PAPER_WIDTH_IN, 2.0))
ax_y = ax_x.twinx()
# ax_x.plot(s, x_true_naive, label='$x_{\\mathrm{naive}}(s)$', lw=1)
# ax_x.plot(s, x_naive, label='$x_{\\mathrm{naive}}^δ(s)$', lw=1)
ax_x.plot(s, x_deblurred, label="$x_α(s)$", lw=1)
ax_x.plot(s, x_true, label="$x(s)$", lw=1, c="k")
ax_y.plot(s, y_true, label="$y(s)$", lw=1, c="tab:orange")

ax_y.grid(False)

# Legend
handles_x, labels_x = ax_x.get_legend_handles_labels()
handles_y, labels_y = ax_y.get_legend_handles_labels()
ax_x.legend(handles_x + handles_y, labels_x + labels_y, loc="upper left")
ax_x.set_xlabel("Spatial Coordinate, $s$")
ax_x.set_ylabel("$x(s)$")
ax_y.set_ylabel("$y(s)$")

ax_x.set_title(
    f"Regularisation of Deblurring Problem given $δ={sci_notation(delta, 0)}$"
)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "3_b_validation.png", dpi=600, bbox_inches="tight")
display(fig)

fig.set_size_inches(PAPER_WIDTH_IN, 3.5)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "3_b_validation_large.png", dpi=600, bbox_inches="tight")
# Investigate convergence
deltas = np.geomspace(1.0, 1e-7, 50)

machine_prec = np.finfo(float).eps
machine_prec_transformed = np.cbrt(machine_prec)

convergence_data = []
for delta in deltas:
    alpha = delta**3

    # Compute regularisation
    y_yj_ip = Y_basis.T @ y_noisy
    filter_factors = mu**2 / (mu**3 + alpha)
    coefficients = filter_factors * y_yj_ip
    x_deblurred = X_basis @ coefficients

    error_norm = np.linalg.norm(x_deblurred - x_true, ord=2)
    convergence_data.append({"delta": delta, "alpha": alpha, "error_norm": error_norm})


df = pd.DataFrame(convergence_data)

fig, ax = plt.subplots(figsize=(PAPER_WIDTH_IN, 2.0))
ax.xaxis.set_inverted(True)
plt.plot(df["delta"], df["error_norm"], marker="o", ls="-", label="$\\|x_α - x\\|_2$")
plt.plot(
    df["delta"],
    4 * df["delta"] ** (1 / 5),
    ls="--",
    lw=1.5,
    label="$\\mathcal{O}(δ^{1/5})$",
)
plt.axvline(
    machine_prec_transformed, color="k", ls="--", lw=1.5, label="$\\sqrt[3]{ε}$"
)
plt.loglog()
plt.grid(axis="y", which="minor", ls=":")
plt.legend()
plt.title("Convergence of Regularisation for Deblurring Problem")
plt.xlabel("Noise Level, $δ$")
plt.ylabel("Reconstruction Error\n$\\|x_α - x\\|_2$")
plt.tight_layout()
fig.savefig(OUTPUT_DIR / "3_b_convergence.png", dpi=600, bbox_inches="tight")
display(fig)

fig.set_size_inches(PAPER_WIDTH_IN, 3.5)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "3_b_convergence_large.png", dpi=600, bbox_inches="tight")
