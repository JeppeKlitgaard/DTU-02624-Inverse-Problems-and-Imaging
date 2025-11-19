

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


alphas = np.linspace(0, 0.015, 200)
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
df_2d = df.copy()

# --- 
# Find optimal alpha for ordinary Tikhonov
df_2d["ordinary_search_target"] = df_2d["F_1_ordinary"] - safety_factor * noise_level
idx_ordinary = df_2d["ordinary_search_target"].abs().idxmin()
df_2d["generalised_search_target"] = (
    df_2d["F_1_generalised"] - safety_factor * noise_level
)
idx_generalised = df_2d["generalised_search_target"].abs().idxmin()
idx_ordinary_true = df_2d["x_diff_ordinary_L2"].abs().idxmin()
idx_generalised_true = df_2d["x_diff_generalised_L2"].abs().idxmin()
print(
    f"Optimal alpha (ordinary): Estimate: {df_2d.iloc[idx_ordinary]['alpha']}, True min: {df_2d.iloc[idx_ordinary_true]['alpha']}"
)
print(
    f"Optimal alpha (generalised): Estimate: {df_2d.iloc[idx_generalised]['alpha']}, True min: {df_2d.iloc[idx_generalised_true]['alpha']}"
)

ordinary_estimate_alpha = df_2d.iloc[idx_ordinary]["alpha"]
ordinary_estimate_residual = df_2d.iloc[idx_ordinary]["F_1_ordinary"]
ordinary_estimate_error = df_2d.iloc[idx_ordinary]["x_diff_ordinary_L2"]
generalised_estimate_alpha = df_2d.iloc[idx_generalised]["alpha"]
generalised_estimate_residual = df_2d.iloc[idx_generalised]["F_1_generalised"]
generalised_estimate_error = df_2d.iloc[idx_generalised]["x_diff_generalised_L2"]
ordinary_true_alpha = df_2d.iloc[idx_ordinary_true]["alpha"]
ordinary_true_residual = df_2d.iloc[idx_ordinary_true]["F_1_ordinary"]
ordinary_true_error = df_2d.iloc[idx_ordinary_true]["x_diff_ordinary_L2"]
generalised_true_alpha = df_2d.iloc[idx_generalised_true]["alpha"]
generalised_true_residual = df_2d.iloc[idx_generalised_true]["F_1_generalised"]
generalised_true_error = df_2d.iloc[idx_generalised_true]["x_diff_generalised_L2"]

generalised_best_x = df_2d.iloc[idx_generalised]["x_alpha_generalised"]
generalised_best_alpha = df_2d.iloc[idx_generalised]["alpha"]

# --- 
# Discrepancy plot
fig, (ax_resid, ax_error) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, 2.5))

sns.lineplot(ax=ax_resid, data=df_2d, x="alpha", y="F_1_generalised")

ax_resid.axhline(noise_level, color="black", ls=":", label="Noise level, $δ$")
ax_resid.axhline(
    noise_level * safety_factor, color="black", ls="--", label="Augmented level, $νδ$"
)

ax_resid.plot(
    generalised_estimate_alpha,
    generalised_estimate_residual,
    "o",
    color="C1",
    markersize=7,
    label="$α_\\mathrm{DP}$",
)
ax_resid.plot(
    generalised_true_alpha,
    generalised_true_residual,
    "*",
    color="C2",
    markersize=9,
    label="$α_\\mathrm{optimal}$",
)

ax_resid.set_xlabel("Regualarisation Parameter, $α$")
ax_resid.set_ylabel("Residual Norm, $||Ax_α - y^δ||$")

sns.lineplot(ax=ax_error, data=df_2d, x="alpha", y="x_diff_generalised_L2", color="C0")
ax_error.plot(
    generalised_estimate_alpha,
    generalised_estimate_error,
    "o",
    color="C1",
    markersize=7,
)
ax_error.plot(
    generalised_true_alpha, generalised_true_error, "*", color="C2", markersize=9
)
# sns.lineplot(data=dfs, x="alpha", y="x_diff_generalised_L2")
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
    "Discrepancy Principle for Generalised Tikhonov Regularisation", fontsize=12, y=0.92
)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "2_d.png", dpi=600, bbox_inches="tight")

# --- 
plt.figure(figsize=(PAPER_WIDTH_IN, 2.5))
plt.plot(
    s,
    ordinary_best_x,
    label=f"Ordinary Tikhonov, $x_{{α={ordinary_best_alpha:.4f}}}$",
    lw=1,
)
plt.plot(
    s,
    generalised_best_x,
    label=f"Generalised Tikhonov, $x_{{α={generalised_best_alpha:.4f}}}$",
    lw=2,
)

plt.plot(s, x_true, label="True Solution, $x$", lw=2.0, ls="--", color="black")
plt.legend()

plt.xlabel("Spatial Coordinate, $s$")
plt.ylabel("Solution, $x(s)$")
plt.title("Tikhonov Regularisation Solutions for Abel Problem")

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "2_d_solutions.png", dpi=600, bbox_inches="tight")
