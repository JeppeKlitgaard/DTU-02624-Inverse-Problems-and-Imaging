def forward_heat_operator(N_int: int, T: float) -> npt.NDArray:
    """
    Creates the discrete forward operator A for the 1D heat equation.

    The operator maps the initial state u(s, 0) to the final state u(s, T).
    It is constructed using the Discrete Sine Transform (DST).

    Note that the spatial domain is assumed to be (0, π).

    Arguments:
        N_int: The number of internal grid points.
        T: The final time.

    Returns:
        Matrix A of shape (N_int, N_int) that can be applied to the _internal_ state vector.
    """

    k = np.arange(N_int) + 1  # Wavenumbers
    decay_factors = np.exp(-(k**2) * T)  # Decay factors for each mode

    A = np.zeros((N_int, N_int))
    for j in range(N_int):
        # Basis vector
        e_j = np.zeros(N_int)
        e_j[j] = 1.0

        # Type I DST corresponding to sine series expansion
        coeffs = dst(e_j, type=1, norm="ortho")
        decayed_coeffs = decay_factors * coeffs
        col_j = idst(decayed_coeffs, type=1, norm="ortho")

        A[:, j] = col_j

    return A


# Parameters
N = 128
N_int = N - 2
T = 0.1
a = 0
b = np.pi

# Get spatial discretisation
s_grid = np.linspace(a, b, N, endpoint=True)
s_internal = s_grid[1:-1]  # Internal points only


def analytical_1(s: npt.NDArray) -> npt.NDArray:
    return np.sin(4 * s)


x_initial = analytical_1(s_grid)
x_initial_internal = x_initial[1:-1]

# Get operator
A = forward_heat_operator(N_int, T)
y_final_numeric_internal = A @ x_initial_internal

# --- 
np.random.seed(42)

gamma = 1e-1
sigma = 1e-3

N_simulations = 10
I = np.eye(*A.shape)

fig, (ax_coarse, ax_fine) = plt.subplots(
    1, 2, figsize=(PAPER_WIDTH_IN, 3), sharey=True, sharex=True
)
for N_simulations, ax in [(5, ax_coarse), (500, ax_fine)]:
    posterior_covariance = (
        gamma**2 * I
        + gamma**4 * A.T @ np.linalg.inv(gamma**2 * A @ A.T + sigma**2 * I) @ A
    )
    posterior_mean = (
        gamma**2
        * A.T
        @ np.linalg.inv(gamma**2 * A @ A.T + sigma**2 * I)
        @ y_final_numeric_internal
    )
    posterior_variance = np.diag(posterior_covariance)
    x_sample = np.random.multivariate_normal(
        posterior_mean, posterior_covariance, size=N_simulations
    )

    sample_mean = np.mean(x_sample, axis=0)
    sample_variance = np.var(x_sample, axis=0)

    ax.plot(
        s_internal,
        x_initial_internal,
        color="black",
        linewidth=3.0,
        label="Initial condition",
    )
    ax.plot(s_internal, sample_mean, color="C0", lw=2.0, label="Sample mean")
    ax.plot(s_internal, posterior_mean, color="C1", ls="--", label="Posterior mean")

    # Plot 95% credible intervals
    sample_upper_bound = sample_mean + 1.96 * np.sqrt(sample_variance)
    sample_lower_bound = sample_mean - 1.96 * np.sqrt(sample_variance)
    posterior_upper_bound = posterior_mean + 1.96 * np.sqrt(posterior_variance)
    posterior_lower_bound = posterior_mean - 1.96 * np.sqrt(posterior_variance)
    ax.fill_between(
        s_internal, sample_lower_bound, sample_upper_bound, color="C0", alpha=0.3
    )
    ax.fill_between(
        s_internal, posterior_lower_bound, posterior_upper_bound, color="C1", alpha=0.3
    )

    ax.set_title(f"$N_\\mathrm{{samples}} = {N_simulations}$")

legend_elements = [
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
    "$x_\\mathrm{true}$",
    "$x_\\mathrm{sample}$\n(95% CI)",
    "$x_\\mathrm{post}$\n(95% CI)",
]

ax.legend(legend_elements, legend_labels, loc="center left", bbox_to_anchor=(1, 0.5))

ax_coarse.set_xlabel("Spatial coordinate, $s$")
ax_fine.set_xlabel("Spatial coordinate, $s$")
ax_coarse.set_ylabel("Initial Temperature, $x(s)$")

fig.suptitle(
    "Analytical and Sampled Bayesian Estimation of Initial Temperature",
    fontsize=12,
    y=0.92,
)
fig.tight_layout()
plt.savefig(OUTPUT_DIR / "1_b.png", dpi=EXPORT_DPI, bbox_inches="tight")

# --- 
# Visualise temperature progression
for T_final in np.linspace(0.001, 0.25, 5):
    A_original = forward_heat_operator(N_int, T_final)
    y = A_original @ x_initial_internal

    plt.plot(s_internal, y, label=f"T={T_final:.3f}")
