img = image.open(HW3_DIR / "KarinHarbo512.png").convert("L")
x = np.asarray(img) / 255.0

T = 0.001
y = forward_heat_operator(x, T=T)

plt.figure()
plt.imshow(x, cmap='gray', norm=gray_norm)
plt.grid(False)
plt.figure()
plt.grid(False)
plt.imshow(y, cmap='gray', norm=gray_norm)

np.random.seed(42)
noise = np.random.normal(0, 0.01, size=y.shape)
y_noisy = y + noise
plt.figure()
plt.grid(False)
plt.imshow(y_noisy, cmap='gray', norm=gray_norm)

# plt.figure()
# plt.grid(False)
# plt.imshow(y - y_noisy, cmap='gray')

# --- 
fig, (ax_x, ax_y_noisy) = plt.subplots(1, 2, figsize=(PAPER_WIDTH_IN, PAPER_WIDTH_IN / 2 - 0.1), sharex=True)

img_x = ax_x.imshow(x, cmap='gray', norm=gray_norm)
img_y_noisy = ax_y_noisy.imshow(y_noisy, cmap='gray', norm=gray_norm)

for ax in (ax_x, ax_y_noisy):
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

# Add colorbars to both axes so they have equal size
# https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
divider_x = make_axes_locatable(ax_x)
cax_x = divider_x.append_axes("right", size="5%", pad=0.05)
cax_x.axis('off')  # Hide the left colorbar axis

divider_y = make_axes_locatable(ax_y_noisy)
cax_y = divider_y.append_axes("right", size="5%", pad=0.05)
plt.colorbar(img_y_noisy, cax=cax_y)

ax_x.set_title("Initial condition, $x(s_1, s_2)$")
ax_y_noisy.set_title("Noisy measurements, $y_\\text{noisy}(t_1, t_2; T)$")

fig.tight_layout()
plt.savefig(OUTPUT_DIR / "1_c.png", dpi=EXPORT_DPI, bbox_inches='tight')
