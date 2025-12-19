#import "../presentation_template.typ": *
#show: doc => setup(1, doc)

= Introduction - Homework 1

- Backwards Heat Equation
- Abel Problem
- Deblurring Problem
- Quick overview of math along the way
- Likely will not be able to cover all
  - Will go through most part quite quickly
  - Feel free to ask more detailed questions during or after
- Feedback received on report
  - Use _FFT_ for better efficiency (backwards heat problem)
  - Use _augmented least-squares_ for improved stability (Tikhonov)

= Backwards Heat Equation

Forward model:
$
  K x(t) = ∫_0^π x(s) k(t,s) dif s wide wide k(t,s) = 2/π sum_(n=1)^∞ e^(-T n^2) sin(n t) sin(n s)
$

Notice higher-order modes $n≫1$ are particularly attenuated\
⇒ inverse problem ill-posed (sensitivity).

== Backwards Heat Problem — Singular System

$
  K x(t) = ∫_0^π x(s) k(t,s) dif s wide wide k(t,s) = 2/π sum_(n=1)^∞ e^(-T n^2) sin(n t) sin(n s)
$

Is given system the singular system?
$
  x_j (t) = sqrt(2\/π) sin(j t), wide y_j (s) = sqrt(2\/π) sin(j s), wide μ_j = e^(-T j^2)
$
We find this satisfies requirements for singular system:
- $x_j, y_j$ span $L^2 (0, π)$ (fourier basis)
- $K x_j (t)
  &= ∫_0^π sqrt(2/π) sin(j s) (2/π ∑_(n=1)^∞ #purple[$e^(-T n^2)$] sin(n t) sin(n s)) dif s\
  &= (2/π)^(3/2) (∑_(n=1)^∞ #purple[$e^(-T n^2)$] sin(n t) ∫_0^π #blue[$sin(j s) sin(n s)$] dif s)
  = (2/π)^(3/2) (∑_(n=1)^∞ #purple[$e^(-T n^2)$] sin(n t) π/2 #blue[$δ_(j n)$])\
  &= #purple[$e^(-T n^2)$] sqrt(2/π) sin(j t)
  = μ_j y_j (t)$
  - Tricks: #blue[orthogonality $sin(j s) sin(n s) = δ_(j n)$]\; #purple[Gaussian decay ⇒ uniform convergence]
- Similarly $K^* y_j = μ_j x_j$ upon finding $K = K^*$ due to symmetric, real kernel: $k^*(s, t) = k(t, s)$
- $μ_j > 0$ and $μ_j → 0 "as" j → ∞$

⇒  _singular system_ ✔️

== Backwards Heat Problem — Discretisation

We discretise operator $K$ into $A ∈ ℝ^(n×n)$ with $n=128$ using either *FFT* $∼ cal(O)(n log n)$, or *sampling of kernel* at grid points $∼ cal(O)(n^2)$.

For this assignment we did the latter, constructing $A = U Σ V^TT$
$
  U_(i j) (t) = sqrt(2/(n+1)) sin (j ⋅ t_j)
  wide
  V_(i j) (s) = sqrt(2/(n+1)) sin (j ⋅ s_j)
  wide
  Σ_(j j ) = e^(-T j^2)
$

*Notes*:
-  Possible as we know singular system analytically.
  - Generally requires quadrature.
- Normalisation in pre-factor
- For other assignments, the FFT approach was used instead
  - Matrix-free implementation used in HW3.

== Backwards Heat Problem — Accuracy of Discretisation

#side-by-side[
  #image(
    "output/1_b_singular_values.png",
  )

  Note cut-off at machine precision, $∼10^(-16)$, corresponding to $j≈33$
][
  #image(
    "output/1_b_singular_vectors.png",
  )
  Last mode is past $j=33$ cut-off, so noise dominates decomposition
]

== Backwards Heat Problem — Truncated SVD
#progress-marker-lower

Given function
$x(s) = s sin(s) exp((s - π / 2)^2)$
invert noisy measurements $y^δ$ where noise level $δ = 2%$

Truncated SVD (TSVD) to regularise inversion, only first $j^⋆$ modes are used.\
$j^⋆=9$ found using *Picard Plot* and *Discrete Picard Condition*: _expansion coefficients must decay faster than singular values._

#align(center, image(
  "output/1_c_picard_plot.png", width: 90%
))

== Backwards Heat Problem — Reconstruction Result

#align(center, image(
  "output/1_c_reconstruction.png",
))

- *Note*: $x_"naive"$ is read on right axis. It represents naive inversion using $A^(-1)$

= Abel Problem
#progress-marker-upper

We consider the Abel Problem defined in $L^2 (0, 1)$:
$
  K x(t) = ∫_0^t x(s) / sqrt(t - s) dif s
$
Proof of compactness skipped here, but outlined as:
- Singularity handled by offset kernel $K_ε$
- Proof then becomes $K_ε → K "as" ε → 0$. Use difference operator, Young's convolution inequality

Noisy measurements are produced analytically using $y^δ = y+n$ where $n$ corresponds to 5% Gaussian noise.

We attempt to reconstruct using both TSVD and Tikhonov regularisation:
$
  x_α = (K^* K + α I)^(-1) K^* y^δ = (A^T A + α I)^(-1) A^T y^δ
$
#speaker-note[Tikhonov has 3 formulations: Variational, Normal Equations, Bayesian. Here we use Normal Equations form.]

Note: improved numerical stability when using *augmented least-squares* formulation, which may be solved using `np.linalg.lstsq`:
$
  tilde(A) = [A; sqrt(α) I], wide tilde(y) = [y^δ; 0]
$

This arises from the Tikhonov variational formulation, $min_x norm(A x - y^δ)_2^2 + α norm(x)_2^2$
#speaker-note[Simply a trick where the two norms are ]

== Abel Problem — TSVD Reconstruction

#side-by-side(columns: (80%, 1fr))[
#image("output/2_b_picard.png")
#image("output/2_b_reconstructions.png")
][
  Notes:
  - Picard ⇒ $j^⋆ ≈ 20$
  - Strong bias as $s→1$
]

== Abel Problem — Tikhonov Reconstruction
#side-by-side(columns: (65%, 1fr))[
#image("output/2_c_tikhonov_reconstructions.png")
#image("output/2_c_tsvd_vs_tikhonov.png")
][
  - Both biased, Tikhonov less so
  - Very smooth true solution, so truncation works well
  - As always, *variance-bias trade-off*
]

== Abel Problem — Proving Bias at $s=1$
#progress-marker-lower

It should be possible to probe bias at $s=1$ analytically.
- First attempt in discrete case failed
- Continuous case easier, outline:
  - Consider $(K^* K x_j)(s=1) = μ_j^2 x_j (s=1)$
  - Find adjoint operator $K^* y(s) = ∫_s^1 y(t)/(sqrt(t-s)) dif t$ using inner-products and Fubini's theorem
  - Thus we find: $(K^* K x_j)(s=1) = μ_j^2 x_j (s=1) = underbrace(#red[$∫_(s=1)^1$], =0) ∫_0^t 1/sqrt(t-τ) y(τ) dif τ = 0$
  - Since $μ_j ≠ 0$, we must have $x_j (s=1) = 0$ for all $j$
- Argue that discrete case must tend to this as $n→∞$

== Abel Problem — Discontinuous Function
We also investigate a discontinuous function:

#image("output/2_d_discontinuous.png")

- High frequency components required to capture discontinuity
- This means less effective noise reduction

= Deblurring Problem
#progress-marker-upper

We are given a regularisation operator $tilde(R)_α$ such that
$
  x_α = tilde(R)_α y = sum_j μ_j^2 / (μ_j^3 + α) (y, y_j) x_j
$

For this to be a regulariser, we require $tilde(R)_α K x = x "as" α → 0$. Outline of proof:
- Decompose $x$ in singular system: $x = sum_j (x, x_j) x_j, wide y = K x = sum_j μ_j (x, x_j) y_j$ (complete basis, injectivity)
- Find $(y, y_j) = μ_j (x, x_j)$ by linearity of inner product (sum-integral swap) and orthonormality of basis
- Consider norm of difference, $norm(tilde(R)_α K x - x)^2 = norm(sum_j [ (μ_j^3 / (μ_j^3 + α)) - 1] (x, x_j) x_j)^2$
- Apply Plancherel theorem: $sum_(k=1)^∞ abs((x, x_k)_X)^2 = norm(x)_X^2$
- Take limit $α →0^+$
- Use _Lebesque's Dominated Convergence Theorem_ to interchange limit and sum
- Thus we find that the term factor $μ_j^3 / (μ_j^3 + α) - 1 → 0$ as $α → 0, ∀j$
- Hence $norm(tilde(R)_α K x - x) → 0$ as $α → 0 wide qed$

== Deblurring Problem — Convergence

We test convergence of regulariser using blurring problem
$
  K x(t) = ∫_(-6)^6 k (t-s) x(s) dif s wide "where" wide k(z) = cases(1 + cos(π z / 3) &quad "if" quad abs(z) < 3, 0 &quad "if" quad abs(z) ≥ 3)
$

In report, we picked $α(δ) = δ^3$, though it should be $α(δ) = δ^2$, which we investigate the convergence of here for $δ=10^(-7)$, which determines the noise level.

#align(center, image("output/3_b_validation.png", width: 80%))

== Deblurring Problem — Convergence Results
#side-by-side(columns: (65%, 1fr))[
#image("output/3_b_validation.png")
#image("output/3_b_convergence.png")
][
  - Good convergence observed
  - Changing from $α=δ^3$ to $α=δ^2$ does not change rate of convergence or results
  - Diverges after $δ≈6×10^(-6)$ due to numerical stability
]
