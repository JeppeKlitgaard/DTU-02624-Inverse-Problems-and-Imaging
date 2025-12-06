N = 100
img = image.open(HW3_DIR / "KarinHarboStripe.png").convert("L")
x = np.asarray(img) / 255.0

assert x.shape == (N, N)

# --- 
plt.figure()
plt.imshow(x, cmap='gray', norm=gray_norm)
plt.xticks([])
plt.yticks([])
plt.title("KarinHarboStripe.png")

plt.savefig(OUTPUT_DIR / "2_a.png", dpi=EXPORT_DPI, bbox_inches='tight')
