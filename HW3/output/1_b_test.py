# Construct test case from task
N = 256
k = 4
l = 3
T = 1

# Homogenous so don't store endpoints
ds = np.pi / (N + 1)
s_1 = np.arange(ds, np.pi, ds)
s_2 = s_1.copy()
s1_grid, s2_grid = np.meshgrid(s_1, s_2, indexing="ij")
x_initial = np.sin(k * s1_grid) * np.sin(l * s2_grid)
y_final = forward_heat_operator(x_initial, T=T)

# Test:
# Since only the k-lth mode is non-zero, only this mode is propagated
decay = np.exp(-(k**2 + l**2) * T)

# Optionally, if we care about being normalised:
# integral_scaling = (np.pi / 2)**2
# y_expected = x_initial * decay * integral_scaling
y_expected = x_initial * decay

plt.figure()
plt.imshow(x_initial)
plt.colorbar()
plt.title("Initial condition")

plt.figure()
plt.imshow(y_final)
plt.colorbar()
plt.title("Final condition (computed)")

plt.figure()
plt.imshow(y_expected)
plt.colorbar()
plt.title("Final condition (expected)")

plt.figure()
plt.imshow(y_final - y_expected)
plt.colorbar()
plt.title("Difference (computed - expected)")

# --- 
fig, (ax, ax_error) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2 - 0.5))

img = ax.imshow(y_expected, plt.cm.viridis)
ax.set_title("Measurements, $y_\\text{exact}(t_1, t_2; T)$")
ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])

# https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
plt.colorbar(img, ax=ax, fraction=0.046, pad=0.04)
img_error = ax_error.imshow(y_final - y_expected, plt.cm.jet)
ax_error.set_title("Error, $e(t_1, t_2; T)$")
ax_error.set_aspect('equal')
ax_error.set_xticks([])
ax_error.set_yticks([])
plt.colorbar(img_error, ax=ax_error, fraction=0.046, pad=0.04)

fig.tight_layout()
fig.suptitle("Forward Heat Operator Validation", y=1.02)

plt.savefig(OUTPUT_DIR / "1_b.png", dpi=EXPORT_DPI, bbox_inches='tight')
