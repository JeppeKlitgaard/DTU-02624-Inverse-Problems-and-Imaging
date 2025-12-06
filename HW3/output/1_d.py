# Since operator is self adjoint, we can easily implement Landweber

def landweber_iteration(x_k, y_noisy, T, omega):
    A_x_k = forward_heat_operator(x_k, T=T)
    return x_k - omega * forward_heat_operator(A_x_k - y_noisy, T=T)

# Largest singular value corresponds to the lowest frequency mode
# exp(-(k^2 + l^2) * T) with k=l=1
singular_value_max = np.exp(-2 * (1**2) * T)
omega = 0.9 * (2 / (singular_value_max**2))  # Safety margin
print(omega)
x_k = y_noisy.copy()
history = [x_k.copy()]
# N_iter = 1500
N_iter = 1500
for i in trange(N_iter):
    x_k = landweber_iteration(x_k, y_noisy, T, omega)
    history.append(x_k.copy())

# --- 
fig, (ax_recon, ax_error) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2 - 0.1), sharex=True)

img_recon = ax_recon.imshow(history[-1], cmap='gray', norm=gray_norm)
img_error = ax_error.imshow(history[-1] - x, cmap='jet')
plt.colorbar(img_recon, ax=ax_recon, fraction=0.046, pad=0.04)
plt.colorbar(img_error, ax=ax_error, fraction=0.046, pad=0.04)

for ax in (ax_recon, ax_error):
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

ax_recon.set_title("Reconstructed Image")
ax_error.set_title("Reconstruction Error")

fig.tight_layout()
plt.savefig(OUTPUT_DIR / "1_d_images.png", dpi=EXPORT_DPI, bbox_inches='tight')

# --- 
data = []
for i, x_k in enumerate(tqdm(history)):
    error_L2 = np.linalg.norm(x_k - x, ord=2) / np.linalg.norm(x, ord=2)
    error_Linf = np.linalg.norm(x_k - x, ord=np.inf) / np.linalg.norm(x, ord=np.inf)
    data.append({"Iteration": i, "Relative Error L2": error_L2, "Relative Error Linf": error_Linf})

df = pd.DataFrame(data)

# --- 
plt.figure(figsize =(PAPER_WIDTH_IN, 2.5))

ax1 = plt.gca()
ax2 = ax1.twinx()

sns.lineplot(data=df, x="Iteration", y="Relative Error L2", label="$L_2$ norm", c="C0", ax=ax1)
sns.lineplot(data=df, x="Iteration", y="Relative Error Linf", label="$L_∞$ norm", c="C1", ax=ax2)

ax1.set_ylabel(r"Relative Error $\|x_k - x_\text{true}\|_2$")
ax2.set_ylabel(r"Relative Error $\|x_j - x_\text{true}\|_∞$")

# Collect labels in one legend
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')
ax2.legend_.remove()

ax2.grid(False)

plt.title("Landweber Iteration Convergence")
ax1.set_xlabel("Iteration, $j$")

plt.savefig(OUTPUT_DIR / "1_d_convergence.png", dpi=EXPORT_DPI, bbox_inches='tight')
