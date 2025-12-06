# Generate mask
mask = x == 0.0

x_inpaint = x.copy()
# x_inpaint[mask] = 0.5  # Initial guess, start at gray
plt.figure()
plt.imshow(mask)
plt.figure()
plt.imshow(x_inpaint, cmap='gray', norm=gray_norm)

# x_

# --- 
u = cp.Variable((N, N))

z = x_inpaint.copy()
M = (~mask).astype(np.float64)

fidelity = 0.5 * cp.sum_squares(cp.multiply(M, u - z))
regularisation = cp.tv(u)

lambda_ = 0.01
objective = cp.Minimize(fidelity + lambda_ * regularisation)
problem = cp.Problem(objective)
result = problem.solve()

u_sol_tv = u.value.copy()

# --- 
plt.imshow(u_sol_tv, cmap='gray', norm=gray_norm)
plt.xticks([])
plt.yticks([])
plt.title("TV Inpainting Result")
plt.savefig(OUTPUT_DIR / "2_b_tv.png", dpi=EXPORT_DPI, bbox_inches='tight')


# --- 
fig, (ax_orig, ax_inpaint) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2))

img = ax_orig.imshow(x_inpaint, cmap="gray", norm=gray_norm)
ax_orig.set_title("Original Image")
ax_orig.set_aspect('equal')
ax_orig.set_xticks([])
ax_orig.set_yticks([])

img_inpaint = ax_inpaint.imshow(u_sol_tv, cmap="gray", norm=gray_norm)
ax_inpaint.set_title("TV Inpainting")
ax_inpaint.set_aspect('equal')
ax_inpaint.set_xticks([])
ax_inpaint.set_yticks([])

fig.tight_layout()
plt.savefig(OUTPUT_DIR / "2_b_tv_noisy_comparison.png", dpi=EXPORT_DPI, bbox_inches='tight')

# --- 
u = cp.Variable((N, N))
objective = cp.Minimize(cp.tv(u))
constraint = [u[~mask] == x_inpaint[~mask]]
problem = cp.Problem(objective, constraint)
result = problem.solve()

u_sol_tv_noiseless = u.value.copy()

# --- 
fig, (ax_orig, ax_inpaint) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2))

img = ax_orig.imshow(x_inpaint, cmap="gray", norm=gray_norm)
ax_orig.set_title("Original Image")
ax_orig.set_aspect('equal')
ax_orig.set_xticks([])
ax_orig.set_yticks([])

img_inpaint = ax_inpaint.imshow(u_sol_tv_noiseless, cmap="gray", norm=gray_norm)
ax_inpaint.set_title("TV Inpainting (Noiseless)")
ax_inpaint.set_aspect('equal')
ax_inpaint.set_xticks([])
ax_inpaint.set_yticks([])

fig.tight_layout()
plt.savefig(OUTPUT_DIR / "2_b_tv_noiseless_comparison.png", dpi=EXPORT_DPI, bbox_inches='tight')

# --- 
fig, axs = plt.subplots(1, 3, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2))

img = ax_orig.imshow(x_inpaint, cmap="gray", norm=gray_norm)
axs[0].imshow(x_inpaint, cmap="gray", norm=gray_norm)
axs[0].set_title("Original Image")

axs[1].imshow(u_sol_tv, cmap="gray", norm=gray_norm)
axs[1].set_title("TV Inpainting ($λ=0.01$)")

axs[2].imshow(u_sol_tv_noiseless, cmap="gray", norm=gray_norm)
axs[2].set_title("TV Inpainting (Noiseless)")

for ax in axs:
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

fig.tight_layout()
plt.savefig(OUTPUT_DIR / "2_b_tv_comparison.png", dpi=EXPORT_DPI, bbox_inches='tight')

# --- 
# Harmonic
u = cp.Variable((N, N))

z = x_inpaint.copy()
M = (~mask).astype(np.float64)
fidelity = 0.5 * cp.sum_squares(cp.multiply(M, u - z))

grad_u_x = cp.diff(u, axis=0)
grad_u_y = cp.diff(u, axis=1)
regularisation =  0.5 * (cp.sum_squares(grad_u_x) + cp.sum_squares(grad_u_y))

lambda_ = 0.1
objective = cp.Minimize(fidelity + lambda_ * regularisation)
problem = cp.Problem(objective)
result = problem.solve()

u_sol_harmonic = u.value.copy()

# --- 
plt.figure(figsize = (60 * mm, 60 * mm))
plt.imshow(u_sol_harmonic, cmap='gray')
plt.xticks([])
plt.yticks([])
plt.title("Harmonic Inpainting Result")
plt.savefig(OUTPUT_DIR / "2_b_harmonic.png", dpi=EXPORT_DPI, bbox_inches='tight')

# --- 
fig, (ax_orig, ax_inpaint) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2))

img = ax_orig.imshow(x_inpaint, cmap="gray", norm=gray_norm)
ax_orig.set_title("Original Image")
ax_orig.set_aspect('equal')
ax_orig.set_xticks([])
ax_orig.set_yticks([])

img_inpaint = ax_inpaint.imshow(u_sol_harmonic, cmap="gray", norm=gray_norm)
ax_inpaint.set_title("Harmonic Inpainting")
ax_inpaint.set_aspect('equal')
ax_inpaint.set_xticks([])
ax_inpaint.set_yticks([])

fig.tight_layout()
plt.savefig(OUTPUT_DIR / "2_b_harmonic_noisy_comparison.png", dpi=EXPORT_DPI, bbox_inches='tight')

# --- 
u = cp.Variable((N, N))
grad_u_x = cp.diff(u, axis=0)
grad_u_y = cp.diff(u, axis=1)
objective = cp.Minimize(cp.sum_squares(grad_u_x) + cp.sum_squares(grad_u_y))
constraint = [u[~mask] == x_inpaint[~mask]]
problem = cp.Problem(objective, constraint)
result = problem.solve()

u_sol_harmonic_noiseless = u.value.copy()

# --- 
fig, (ax_orig, ax_inpaint) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2))

img = ax_orig.imshow(x_inpaint, cmap="gray", norm=gray_norm)
ax_orig.set_title("Original Image")
ax_orig.set_aspect('equal')
ax_orig.set_xticks([])
ax_orig.set_yticks([])

img_inpaint = ax_inpaint.imshow(u_sol_harmonic_noiseless, cmap="gray", norm=gray_norm)
ax_inpaint.set_title("Harmonic Inpainting (Noiseless)")
ax_inpaint.set_aspect('equal')
ax_inpaint.set_xticks([])
ax_inpaint.set_yticks([])

fig.tight_layout()
plt.savefig(OUTPUT_DIR / "2_b_harmonic_noiseless_comparison.png", dpi=EXPORT_DPI, bbox_inches='tight')

# --- 
fig, axs = plt.subplots(1, 3, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2))

img = ax_orig.imshow(x_inpaint, cmap="gray", norm=gray_norm)
axs[0].imshow(x_inpaint, cmap="gray", norm=gray_norm)
axs[0].set_title("Original Image")

axs[1].imshow(u_sol_tv, cmap="gray", norm=gray_norm)
axs[1].set_title("Harmonic Inpainting, ($λ=0.1$)")

axs[2].imshow(u_sol_tv_noiseless, cmap="gray", norm=gray_norm)
axs[2].set_title("Harmonic Inpainting (Noiseless)")

for ax in axs:
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

fig.tight_layout()
plt.savefig(OUTPUT_DIR / "2_b_harmonic_comparison.png", dpi=EXPORT_DPI, bbox_inches='tight')
