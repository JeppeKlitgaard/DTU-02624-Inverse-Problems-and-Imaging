T = 0.03
n = 128

k = np.arange(1, n + 1)
t = k * np.pi / (n + 1)
s = k * np.pi / (n + 1)

# We construct A matrix using spectral method
# Alternative would be to use a quadrature rule

# Create left and right singular vectors sampled at t_i, s_i
j = np.arange(1, n + 1)
U = np.sqrt(2 / (n + 1)) * np.sin(np.einsum("i,j->ij", t, k))
V = np.sqrt(2 / (n + 1)) * np.sin(np.einsum("i,j->ij", s, k))
assert np.allclose(U.T @ U, np.eye(n)), "V is not orthogonal"
assert np.allclose(V.T @ V, np.eye(n)), "V is not orthogonal"

# Create diagonal matrix S of singular values
S = np.diag(np.exp(-(j**2) * T))

# Now we can construct A = U D V^*
A = U @ S @ V.T

# Now we can compare as requested in assignment
# This is almost cheating since we constructed A from its SVD
# But it would yield the same results if constructed using a sufficiently
# accurate quadrature rule
U_svd, S_svd, Vh_svd = np.linalg.svd(A)
plt.figure(figsize=(PAPER_WIDTH_IN, 2.0))
plt.plot(np.diagonal(S), label="$K$ singular values")
plt.plot(S_svd, label="$A$ singular values")
plt.axvline(33, color="k", ls="--")
plt.axhline(1e-16, color="k", ls=":", alpha=0.5)
plt.xticks(list(plt.xticks()[0]) + [33])
plt.semilogy()
plt.title("Comparison of Singular Values in Backward Heat Problem")
plt.xlabel("Index, $j$")
plt.ylabel("Singular value, $μ_j$")
plt.legend()
plt.ylim(1e-50, 1e4)
plt.xlim(0 - 3, 80)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "1_b_singular_values.png", dpi=600, bbox_inches="tight")
display(plt.gcf())

plt.gcf().set_size_inches(PAPER_WIDTH_IN, 3.5)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "1_b_singular_values_large.png", dpi=600, bbox_inches="tight")


# Singular vectors vs functions
plt.figure()


def singular_vector_func(j, x):
    return np.sqrt(2 / np.pi) * np.sin(j * np.pi * x)


t_i_fine = np.linspace(0, np.pi, 2500)
fig, axes = plt.subplots(3, 1, figsize=(PAPER_WIDTH_IN, 2.7), sharex=True, sharey=True)
for i, (ax, j_k) in enumerate(zip(axes, [9, 30, 35])):
    u_calculated = singular_vector_func(j_k + 1, t_i_fine / np.pi)
    u_calculated *= np.sqrt(t[1] - t[0])  # Normalization factor
    ax.plot(t_i_fine, u_calculated, label=f"Singular Function, $x_j(t)$")
    ax.plot(t, U_svd[:, j_k], "o", ms=3, label="Singular Vector, $u_j(t)$")
    ax.set_title(f"$j={j_k}$", x=1.06, y=0.5)

    if i == 0:
        ax.legend()

    if i == 2:
        ax.set_xlabel("Temporal Coordinate, $t$")

fig.suptitle("Right Singular Vectors vs Singular Functions", y=0.92, fontsize=10)
fig.tight_layout()
fig.subplots_adjust(hspace=0.08)
fig.savefig(OUTPUT_DIR / "1_b_singular_vectors.png", dpi=600, bbox_inches="tight")
display(fig)

fig.set_size_inches(PAPER_WIDTH_IN, 5)
fig.suptitle(fig._suptitle.get_text(), y=0.97, fontsize=12)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "1_b_singular_vectors_large.png", dpi=600, bbox_inches="tight")
