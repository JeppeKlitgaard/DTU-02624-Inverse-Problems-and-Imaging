# Parameters
N = 128
N_int = N - 2
T = 0.1
a = 0
b = np.pi

gamma = 1e-1
sigma = 1e-3

# Get spatial discretisation
s_grid = np.linspace(a, b, N, endpoint=True)
s_internal = s_grid[1:-1]  # Internal points only


def analytical_1(s: npt.NDArray) -> npt.NDArray:
    return np.sin(4 * s)


x_initial = analytical_1(s_grid)
x_initial_internal = x_initial[1:-1]

# Get operator
A = forward_heat_operator(N_int, T)

# We only observe some of the terminal conditions
obs_min = np.pi / 4
obs_max = 3 / 4 * np.pi
obs_mask = (obs_min < s_internal) & (s_internal < obs_max)

y_full = A @ x_initial_internal
y_obs = y_full[obs_mask]

A_obs = A[obs_mask, :]  # Observed rows only

# --- 
M = A_obs.shape[0]  # Number of observations
I_obs = np.eye(M)
I_state = np.eye(N_int)

posterior_mean_obs = (
    gamma**2
    * A_obs.T
    @ np.linalg.inv(gamma**2 * A_obs @ A_obs.T + sigma**2 * I_obs)
    @ y_obs
)
posterior_covariance_obs = (
    gamma**2 * I_state
    + gamma**4
    * A_obs.T
    @ np.linalg.inv(gamma**2 * A_obs @ A_obs.T + sigma**2 * I_obs)
    @ A_obs
)
posterior_variance_obs = np.diag(posterior_covariance_obs)

# Sample posterior
N_simulations = 50
np.random.seed(42)
x_samples_obs = np.random.multivariate_normal(
    posterior_mean_obs, posterior_covariance_obs, size=N_simulations
)

sample_mean_obs = np.mean(x_samples_obs, axis=0)
sample_variance_obs = np.var(x_samples_obs, axis=0)

fig, ax = plt.subplots(figsize=(PAPER_WIDTH_IN, 3))
ax.plot(
    s_internal,
    x_initial_internal,
    color="black",
    linewidth=3.0,
    label="Initial condition",
)
ax.plot(s_internal, sample_mean_obs, color="C0", lw=2.0, label="Sample mean")
ax.plot(s_internal, posterior_mean_obs, color="C1", ls="--", label="Posterior mean")

# Plot 95% credible intervals
sample_upper_bound = sample_mean_obs + 1.96 * np.sqrt(sample_variance_obs)
sample_lower_bound = sample_mean_obs - 1.96 * np.sqrt(sample_variance_obs)
posterior_upper_bound = posterior_mean_obs + 1.96 * np.sqrt(posterior_variance_obs)
posterior_lower_bound = posterior_mean_obs - 1.96 * np.sqrt(posterior_variance_obs)
ax.fill_between(
    s_internal, sample_lower_bound, sample_upper_bound, color="C0", alpha=0.3
)
ax.fill_between(
    s_internal, posterior_lower_bound, posterior_upper_bound, color="C1", alpha=0.3
)

# Plot observation region
ax.axvspan(obs_min, obs_max, color="gray", alpha=0.2, label="Observation region")

legend_elements = [
    mpl.patches.Patch(facecolor="gray", alpha=0.2),
    mpl.lines.Line2D([0], [0], color="black", linewidth=3.0),
    (
        mpl.lines.Line2D([0], [0], color="C0", linewidth=2.0, ls="-"),
        mpl.patches.Patch(facecolor="C0", alpha=0.3),
    ),
    (
        mpl.lines.Line2D([0], [0], color="C1", ls="--"),
        mpl.patches.Patch(facecolor="C1", alpha=0.3),
    ),
]

legend_labels = [
    "Observed",
    "$x_\\mathrm{true}$",
    "$x_\\mathrm{sample}$\n(95% CI)",
    "$x_\\mathrm{post}$\n(95% CI)",
]

ax.legend(legend_elements, legend_labels, loc="center left", bbox_to_anchor=(1, 0.5))

ax.set_xlabel("Spatial coordinate, $s$")
ax.set_ylabel("Initial Temperature, $x(s)$")

fig.suptitle(
    "Analytical and Sampled Bayesian Estimation of Initial Temperature",
    fontsize=12,
    y=0.92,
)
fig.tight_layout()
plt.savefig(OUTPUT_DIR / "1_c.png", dpi=EXPORT_DPI, bbox_inches="tight")

# --- 
# Double plot

# Parameters
N = 128
N_int = N - 2
T = 0.1
a = 0
b = np.pi

gamma = 1e-1
sigma = 1e-3

# Get spatial discretisation
s_grid = np.linspace(a, b, N, endpoint=True)
s_internal = s_grid[1:-1]  # Internal points only


def analytical_1(s: npt.NDArray) -> npt.NDArray:
    return np.sin(4 * s)


x_initial = analytical_1(s_grid)
x_initial_internal = x_initial[1:-1]

np.random.seed(42)

fig, (ax_long, ax_short) = plt.subplots(
    1, 2, figsize=(PAPER_WIDTH_IN, 3), sharey=True, sharex=True
)
for T, ax in [(0.1, ax_long), (0.005, ax_short)]:
    # Get operator
    A = forward_heat_operator(N_int, T)

    # We only observe some of the terminal conditions
    obs_min = np.pi / 4
    obs_max = 3 / 4 * np.pi
    obs_mask = (obs_min < s_internal) & (s_internal < obs_max)

    y_full = A @ x_initial_internal
    y_obs = y_full[obs_mask]

    A_obs = A[obs_mask, :]  # Observed rows only

    M = A_obs.shape[0]  # Number of observations
    I_obs = np.eye(M)
    I_state = np.eye(N_int)

    posterior_mean_obs = (
        gamma**2
        * A_obs.T
        @ np.linalg.inv(gamma**2 * A_obs @ A_obs.T + sigma**2 * I_obs)
        @ y_obs
    )
    posterior_covariance_obs = (
        gamma**2 * I_state
        + gamma**4
        * A_obs.T
        @ np.linalg.inv(gamma**2 * A_obs @ A_obs.T + sigma**2 * I_obs)
        @ A_obs
    )
    posterior_variance_obs = np.diag(posterior_covariance_obs)

    # Sample posterior
    N_simulations = 50
    x_samples_obs = np.random.multivariate_normal(
        posterior_mean_obs, posterior_covariance_obs, size=N_simulations
    )

    sample_mean_obs = np.mean(x_samples_obs, axis=0)
    sample_variance_obs = np.var(x_samples_obs, axis=0)

    ax.plot(
        s_internal,
        x_initial_internal,
        color="black",
        linewidth=1.5,
        label="Initial condition",
    )
    ax.plot(s_internal, sample_mean_obs, color="C0", lw=2.0, label="Sample mean")
    ax.plot(s_internal, posterior_mean_obs, color="C1", ls="--", label="Posterior mean")

    # Plot 95% credible intervals
    sample_upper_bound = sample_mean_obs + 1.96 * np.sqrt(sample_variance_obs)
    sample_lower_bound = sample_mean_obs - 1.96 * np.sqrt(sample_variance_obs)
    posterior_upper_bound = posterior_mean_obs + 1.96 * np.sqrt(posterior_variance_obs)
    posterior_lower_bound = posterior_mean_obs - 1.96 * np.sqrt(posterior_variance_obs)
    ax.fill_between(
        s_internal, sample_lower_bound, sample_upper_bound, color="C0", alpha=0.3
    )
    ax.fill_between(
        s_internal, posterior_lower_bound, posterior_upper_bound, color="C1", alpha=0.3
    )

    # Plot observation region
    ax.axvspan(obs_min, obs_max, color="gray", alpha=0.2, label="Observation region")

    ax.set_title(f"$T = {T}$")

legend_elements = [
    mpl.patches.Patch(facecolor="gray", alpha=0.2),
    mpl.lines.Line2D([0], [0], color="black", linewidth=1.5),
    (
        mpl.lines.Line2D([0], [0], color="C0", linewidth=2.0, ls="-"),
        mpl.patches.Patch(facecolor="C0", alpha=0.3),
    ),
    (
        mpl.lines.Line2D([0], [0], color="C1", ls="--"),
        mpl.patches.Patch(facecolor="C1", alpha=0.3),
    ),
]

legend_labels = [
    "Observed",
    "$x_\\mathrm{true}$",
    "$x_\\mathrm{sample}$\n(95% CI)",
    "$x_\\mathrm{post}$\n(95% CI)",
]

ax.legend(legend_elements, legend_labels, loc="center left", bbox_to_anchor=(1, 0.5))

ax_short.set_xlabel("Spatial coordinate, $s$")
ax_long.set_xlabel("Spatial coordinate, $s$")
ax_long.set_ylabel("Initial Temperature, $x(s)$")

fig.suptitle(
    "Analytical and Sampled Bayesian Estimation of Initial Temperature",
    fontsize=12,
    y=0.92,
)
fig.tight_layout()
plt.savefig(OUTPUT_DIR / "1_c_double.png", dpi=EXPORT_DPI, bbox_inches="tight")
