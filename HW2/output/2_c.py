# Abel Problem
N = 1024
a = 0
b = 1
h = (b - a) / N
s = a + (np.arange(N) + 0.5) * h
t = s + h / 2  # Shift t values to prevent singularity at t=s

s_mesh, t_mesh = np.meshgrid(s, t)

A = np.zeros((N, N))
indices = s_mesh < t_mesh  # Select according to integral
A[indices] = h / np.sqrt((t_mesh[indices] - s_mesh[indices]))


def x_2_func(s):
    return s ** (3 / 4) - 1


def grad_x_2_func(s):
    return (3 / 4) * s ** (-1 / 4)


x_true = x_2_func(s)
grad_x_true = grad_x_2_func(s)
y_true = A @ x_true

# Add noise, 10%
target_noise_level = 0.10
unscaled_noise = np.random.randn(*y_true.shape)
unscaled_noise_level = np.linalg.norm(unscaled_noise)
noise_coeff = target_noise_level * np.linalg.norm(y_true) / unscaled_noise_level
noise = noise_coeff * unscaled_noise
noise_level = np.linalg.norm(noise)
print(f"Achieved noise level: {noise_level}")
y_noisy = y_true + noise

# Plot noisy measurements
plt.figure(figsize=(PAPER_WIDTH_IN, 2.2))
plt.plot(t, y_noisy, label="$y^δ(t)$", lw=1)
plt.plot(t, y_true, label="$y(t)$", lw=1.5, ls="--")
plt.title("Noisy Measurements for Abel Problem given $x(s) = s^2$")
plt.xlabel("Temporal Coordinate, $t$")
plt.legend()
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "2_a.png", dpi=600, bbox_inches="tight")
display(plt.gcf())

plt.gcf().set_size_inches(PAPER_WIDTH_IN, 3.5)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "2_a_large.png", dpi=600, bbox_inches="tight")

# --- 


def create_differentiation_matrix(n: int, h: float) -> npt.NDArray:
    D = np.zeros((n, n))

    # Use forward difference for the first point
    D[0, 0] = -1.0 / h
    D[0, 1] = 1.0 / h

    # Central difference for interior points
    for i in range(1, n - 1):
        D[i, i - 1] = -0.5 / h
        D[i, i + 1] = 0.5 / h

    # Backward difference for the last point
    D[n - 1, n - 2] = -1.0 / h
    D[n - 1, n - 1] = 1.0 / h

    return D


D = create_differentiation_matrix(N, h)

grad_x_estimate = D @ x_true
plt.plot(s, grad_x_true, label="True Gradient", lw=1.5)
plt.plot(s, grad_x_estimate, label="Estimated Gradient", lw=1)
plt.figure()
plt.plot(s, np.abs(grad_x_true - grad_x_estimate), label="Absolute Error", lw=1)

# --- 


def tikhonov_ordinary(alpha: float) -> npt.NDArray:
    n = A.shape[1]
    x_alpha = np.linalg.inv(A.T @ A + alpha * np.eye(n)) @ A.T @ y_noisy
    return x_alpha


def tikhonov_generalised_smoothness(alpha: float) -> npt.NDArray:
    # TODO: Derive in report
    n = A.shape[1]
    x_alpha = np.linalg.inv(A.T @ A + alpha * (np.eye(n) + D.T @ D)) @ A.T @ y_noisy
    return x_alpha


def F_1(alpha, x_alpha):
    return np.linalg.norm(A @ x_alpha - y_noisy)


alphas = np.linspace(0, 0.2, 200)
data = []
for alpha in tqdm(alphas):
    x_alpha_ordinary = tikhonov_ordinary(alpha)
    F_1_val_ordinary = F_1(alpha, x_alpha_ordinary)
    x_diff_ordinary = x_alpha_ordinary - x_true

    x_alpha_generalised = tikhonov_generalised_smoothness(alpha)
    F_1_val_generalised = F_1(alpha, x_alpha_generalised)
    x_diff_generalised = x_alpha_generalised - x_true

    data.append(
        {
            "alpha": alpha,
            "F_1_ordinary": F_1_val_ordinary,
            "x_alpha_ordinary": x_alpha_ordinary,
            "x_diff_ordinary": x_diff_ordinary,
            "x_diff_ordinary_L2": np.linalg.norm(x_diff_ordinary),
            "F_1_generalised": F_1_val_generalised,
            "x_alpha_generalised": x_alpha_generalised,
            "x_diff_generalised": x_diff_generalised,
            "x_diff_generalised_L2": np.linalg.norm(x_diff_generalised),
        }
    )

df = pd.DataFrame(data)

df_2c = df.copy()

# --- 
# Safety factor, referred to as nu in lecture notes
safety_factor = 1.05

# Find optimal alpha for ordinary Tikhonov
df_2c["ordinary_search_target"] = df_2c["F_1_ordinary"] - safety_factor * noise_level
idx_ordinary = df_2c["ordinary_search_target"].abs().idxmin()
df_2c["generalised_search_target"] = (
    df_2c["F_1_generalised"] - safety_factor * noise_level
)
idx_generalised = df_2c["generalised_search_target"].abs().idxmin()
idx_ordinary_true = df_2c["x_diff_ordinary_L2"].abs().idxmin()
idx_generalised_true = df_2c["x_diff_generalised_L2"].abs().idxmin()
print(
    f"Optimal alpha (ordinary): Estimate: {df_2c.iloc[idx_ordinary]['alpha']}, True min: {df_2c.iloc[idx_ordinary_true]['alpha']}"
)
print(
    f"Optimal alpha (generalised): Estimate: {df_2c.iloc[idx_generalised]['alpha']}, True min: {df_2c.iloc[idx_generalised_true]['alpha']}"
)

ordinary_estimate_alpha = df_2c.iloc[idx_ordinary]["alpha"]
ordinary_estimate_residual = df_2c.iloc[idx_ordinary]["F_1_ordinary"]
ordinary_estimate_error = df_2c.iloc[idx_ordinary]["x_diff_ordinary_L2"]
generalised_estimate_alpha = df_2c.iloc[idx_generalised]["alpha"]
generalised_estimate_residual = df_2c.iloc[idx_generalised]["F_1_generalised"]
generalised_estimate_error = df_2c.iloc[idx_generalised]["x_diff_generalised_L2"]
ordinary_true_alpha = df_2c.iloc[idx_ordinary_true]["alpha"]
ordinary_true_residual = df_2c.iloc[idx_ordinary_true]["F_1_ordinary"]
ordinary_true_error = df_2c.iloc[idx_ordinary_true]["x_diff_ordinary_L2"]
generalised_true_alpha = df_2c.iloc[idx_generalised_true]["alpha"]
generalised_true_residual = df_2c.iloc[idx_generalised_true]["F_1_generalised"]
generalised_true_error = df_2c.iloc[idx_generalised_true]["x_diff_generalised_L2"]

ordinary_best_x = df_2c.iloc[idx_ordinary]["x_alpha_ordinary"]
ordinary_best_alpha = df_2c.iloc[idx_ordinary]["alpha"]

# --- 
# Discrepancy plot
plt.figure(figsize=(PAPER_WIDTH_IN, 2))
sns.lineplot(data=df, x="alpha", y="F_1_ordinary", label="$||Ax_α- yδ||$ (ordinary)")
# sns.lineplot(data=df, x="alpha", y="F_1_generalised")

plt.axhline(noise_level, color="black", ls=":", label="Noise level")
plt.axhline(
    noise_level * safety_factor,
    color="black",
    ls="--",
    label="Noise level (safety factor)",
)

plt.plot(
    ordinary_estimate_alpha,
    ordinary_estimate_residual,
    "o",
    color="C0",
    markersize=7,
    label="$α_\\mathrm{DP}$",
)
plt.plot(
    ordinary_true_alpha,
    ordinary_true_residual,
    "*",
    color="tab:red",
    markersize=9,
    label="$α_\\mathrm{optimal}$",
)
plt.legend()

plt.title("Discrepancy Principle for Ordinary Tikhonov Regularisation")
plt.xlabel("Regualarisation Parameter, $α$")
plt.ylabel("Residual Norm, $||Ax_α - y^δ||$")

plt.figure()
sns.lineplot(data=df, x="alpha", y="x_diff_ordinary_L2")
plt.plot(
    ordinary_estimate_alpha,
    ordinary_estimate_error,
    "o",
    color="C0",
    markersize=7,
    label="$α_\\mathrm{DP}$",
)
plt.plot(
    ordinary_true_alpha,
    ordinary_true_error,
    "*",
    color="tab:red",
    markersize=9,
    label="$α_\\mathrm{optimal}$",
)
# sns.lineplot(data=df, x="alpha", y="x_diff_generalised_L2")
plt.semilogy()

# --- 
# Discrepancy plot
fig, (ax_resid, ax_error) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, 2.5))

sns.lineplot(ax=ax_resid, data=df_2c, x="alpha", y="F_1_ordinary")
# sns.lineplot(data=df_2c, x="alpha", y="F_1_generalised")

ax_resid.axhline(noise_level, color="black", ls=":", label="Noise level, $δ$")
ax_resid.axhline(
    noise_level * safety_factor, color="black", ls="--", label="Augmented level, $νδ$"
)

ax_resid.plot(
    ordinary_estimate_alpha,
    ordinary_estimate_residual,
    "o",
    color="C1",
    markersize=7,
    label="$α_\\mathrm{DP}$",
)
ax_resid.plot(
    ordinary_true_alpha,
    ordinary_true_residual,
    "*",
    color="C2",
    markersize=9,
    label="$α_\\mathrm{optimal}$",
)

ax_resid.set_xlabel("Regualarisation Parameter, $α$")
ax_resid.set_ylabel("Residual Norm, $||Ax_α - y^δ||$")

sns.lineplot(ax=ax_error, data=df_2c, x="alpha", y="x_diff_ordinary_L2", color="C0")
ax_error.plot(
    ordinary_estimate_alpha, ordinary_estimate_error, "o", color="C1", markersize=7
)
ax_error.plot(ordinary_true_alpha, ordinary_true_error, "*", color="C2", markersize=9)
# sns.lineplot(data=df_2c, x="alpha", y="x_diff_generalised_L2")
ax_error.set_yscale("log")
ax_error.grid(axis="y", which="both")
ax_error.set_xlabel("Regualarisation Parameter, $α$")
ax_error.set_ylabel("Error Norm, $||x_α - x||$")

# Get legends and handles from both axes
handles1, labels1 = ax_resid.get_legend_handles_labels()
handles2, labels2 = ax_error.get_legend_handles_labels()

handles = [handles1[0]] + handles2 + handles1[1:]
labels = [labels1[0]] + labels2 + labels1[1:]

ax_error.legend(handles, labels)
ax_resid.legend().remove()

fig.suptitle(
    "Discrepancy Principle for Ordinary Tikhonov Regularisation", fontsize=12, y=0.92
)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "2_c.png", dpi=600, bbox_inches="tight")

# --- 
x_alpha_optim_ordinary = df_2c.iloc[idx_ordinary]["x_alpha_ordinary"]
x_alpha_optim_generalised = df_2c.iloc[idx_generalised]["x_alpha_generalised"]
plt.plot(s, x_true, label="True solution", lw=2.0)
plt.plot(s, x_alpha_optim_ordinary, label="Ordinary Tikhonov", lw=1.5)
plt.plot(s, x_alpha_optim_generalised, label="Generalised Tikhonov", lw=1.5)
plt.title("Tikhonov Regularisation Solutions for Abel Problem")
plt.legend()
