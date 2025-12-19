#import "../presentation_template.typ": *
#import "@preview/physica:0.9.7": TT, pdv, dv
#show: doc => setup(2, doc)

= Introduction - Homework 2

- Backwards Heat Equation
  - Bayesian Framework
  - Reduced Observational Region
- Abel Problem and Generalised Tikhonov
  - Uniqueness and Existence
  - Morozov Discrepancy Principle
- Feedback received on report
  - Use _augmented least-squares_ for improved stability (Tikhonov)
  - Discuss null-space of central difference matrix

= Backwards Heat Equation

We seek to model the 1D backwards heat equation
$
  pdv(u, t) - pdv(u, x, 2) = 0, wide 0 < x < π, quad t ∈ (0, T]\
  u(0, t) = u(π, t) = 0, wide u(x, 0) = f(x).
$

Which we model using a forward operator $K$ such that $y(t) = K x(t)$, where $K$ may be discretised on $n=128$ grid points as $A ∈ ℝ^(n×n)$.

We consider the problem using a _Bayesian framework_, $Y = A x + E$ where $x$ is a realisation of $X∼𝒩(0, γ^2 I)$ with $γ=10^(-1)$ and the noise enters as $E∼𝒩(0, σ^2 I), quad σ=10^(-3)$. Thus, an observation is affected by noise as $y_"obs" = y + e$

== Backwards Heat Equation — Probabilistic Formulation

We obtain:
$
  π_"pr" (x) &∝ exp(-x^2/(2 γ^2)) \
  π(y_"obs"|x) &∝ exp(-(y_"obs" - A x)^2/(2 σ^2)) \
  π(y_"obs") &= ∫ π(y_"obs"|x) π_"pr" (x) dif x ∝ 1 \
$

Which using _Bayes' theorem_ gives the posterior distribution as
$
  π_"post"
  &= π(x|y_"obs") = (π(y_"obs"|x) π_"pr" (x))/(π(y_"obs")) ∝ #red[$exp(-1/2 (norm(y_"obs" - A x)^2_2/σ^2 + norm(x)^2_2/γ^2))$] \
  &∝ #blue[$exp(-1/2 (x - dash(x))^TT Γ^(-1)_"post" (x - dash(x)))$] \
$ <eq:posterior_gaussian_proposal>

Thus the posterior is also Gaussian with mean $dash(x)$ and covariance $Γ_"post"$.

== Backwards Heat Equation — Posterior Mean and Covariance

Algebraically arduous to find closed-form expressions for $dash(x)$ and $Γ_"post"$, outlined as:
- Expand #red[red] expression on last slide and compare against #blue[blue] expression
- Quadratic part becomes inverse covariance: $Γ^(-1)_"post" = (A^TT A)/σ^2 + I/γ^2$
- Rest gives mean: $dash(x) = Γ_"post" A^TT y_"obs"/σ^2$
- Obtain $Γ_"post"$ by inverting $Γ^(-1)_"post"$ using _Woodbury matrix identity_:
  - $Γ_"post" = γ^2 I - #purple[$γ^4 A^TT (σ^2 + γ^2 A A^TT)^(-1) A$]$
  - $dash(x) = γ^2 A^TT (σ^2 I + γ^2 A A^TT)^(-1) y_"obs"$
  - Note: Can recover normal form Tikhonov by letting $α = σ^2\/γ^2$ and algebra
  - $dash(x) = (A^TT A + α I)^(-1) A^TT y_"obs"$

Can also prove that incorporation of observations reduces uncertainty, $Γ_"post" < γ^2 I$
- Note #purple[purple] term is _symmetric_ and show _positive definiteness_ by considering quadratic form and transformations $z = A y$, $B = (σ^2 + γ^2 A A^TT)$
- We find that $B$ is positive definite ⇒ $B^(-1)$ is positive definite ⇒ $Γ_"post" < γ^2 I wide qed$

== Backwards Heat Equation — Monte Carlo Sampling

#image("output/1_b.png")

- For $N_"samples" ≫ 1$ we get very good approximation of expected posterior distribution
- Mean estimates surprisingly good even for low $N_"samples"$

== Backwards Heat Equation — Reduced Observational Region
#progress-marker-lower

We reduce our observational region to $s ∈ [π/4, (3π)/4]$ instead of full $[0, π]$.
- Unit diffusivity, terminal time $T$ determines ill-posedness
- Information propagates between observed and unobserved regions with time
- Smaller $T$ ⇒ less propagation ⇒ harder to infer unobserved region (falls back to prior)

#v(1fr)
#align(center, image("output/1_c_double.png"))

= Abel Problem
#progress-marker-upper

We consider Abel problem with integral operator $K$ such that $y(t) = K x(t)$:
$
  K x(t) = ∫_0^1 k(t, s) x(s) dif s wide k(t, s) = cases(1/sqrt(t - s) &quad s < t, 0 &quad s ≥ t) wide x(s) = s^(3\/4) - 1
$
which is discretised into $A∈ℝ^(N × N)$ where $N=1024$.

We introduce measurement noise of 10% relative to norm of $y_"obs"$ to obtain $y^δ = A x + n$.

We find $x(s) = s^(3\/4) - 1 ∈ H^1(0, 1)$
- $x(s) ∈ L^2 (0, 1) quad ∧ quad  x'(s) ∈ L^2 (0, 1) quad ⇒ quad x(s) ∈ H^1(0, 1)$
- Weak derivative is just derivative since $x$ is smooth
- Square integrals fairly trivial to compute

== Abel Problem — Generalised Tikhonov Regularisation

We seek to prove existence of uniqueness of solution to generalised Tikhonov problem
$
  x_α = limits("arg min")_(x∈H^1 (0, 1)) med norm(A x - y^δ)_2^2 + α (norm(x)_2^2 + norm(D x)_2^2)
$ <eq:generalised_tikhonov_problem>

Existence:
- _Weierstrass Extreme Value Theorem_: continuous function on non-empty compact set attains global extrema on set
- Continous by inspection
- Sublevel sets are bounded by coercivity of objective, i.e. $f(x) → +∞ "as" norm(x)_2 → +∞$
  - Derivative and residual terms are non-negative by construction
  - $α norm(x)_2^2$ is trivially coercive for $α > 0$

#speaker-note[
  From what I recall the trick here is to show that the minimiser cannot be on the boundary due to coercivity, so must be on a sublevel set. We pick that one, thus minimiser is on there by Weierstrass EVT.
]

Uniqueness:
- Convexity analysis (02611 Optimisation for DS, alternatively 2nd order condition)
- Sum of convex functions is convex; norms are convex
- Non-negative scaling, raising to power $p ≥ 1$ preserves convexity
- Strict convexity from Hessian: $∇^2_x norm(x)_2^2 = 2I ≻ 0$, thus strictly convex by second order condition

== Abel Problem — Morozov Discrepancy Principle
#progress-marker-lower

An optimal regularisation parameter $α^⋆$ which minimises the reconstruction error $norm(x_a - x)_2$. This is not available in practice as $x$ is unknown.

Instead use _Morozov Discrepancy Principle_, which requires only estimate of noise level $δ$:
- Find $α_"DP"$ such that $norm(K x_(α_"DP") - y^δ) ≈ ν δ$
- $ν = 1.05 ≥ 1$ is safety factor to avoid under-regularisation
- Here found using grid search, can also be found using root-finding methods

#side-by-side(columns: (1fr, 70%))[
#set align(horizon)
- Slightly undershoots optimal value
- Note flatness of curve near minimum
  - This is problem dependent
  - Ill-posed = very flat
][
#align(center, image("output/2_c.png"))
]

== Abel Problem — MDP on Generalised Tikhonov

#side-by-side(columns: (70%, 1fr))[
In previous slide we had $D = 0$, i.e. ordinary Tikhonov.

Now consider generalised Tikhonov with $D$ as central difference matrix.

Note: Null-space of $D$ are all constant vectors ⇒ Derivative term does not penalise/bias constant offsets

#v(1fr)
#align(center, image("output/2_d.png"))
][
$
  D &= 1/h mat(
    -1, 1;
    -1/2, 0, 1/2;
    , ⋱, ⋱, ⋱;
    ,, -1/2, 0, 1/2;
    ,,, -1, 1;
  ),\
quad h &= 1/N
$

#v(1fr)
- Overshoots optimal value
- Lower errors possible than with ordinary Tikhonov
  - Better able to exploit smoothness of true solution
]

== Abel Problem — Reconstruction Results
#image("output/2_d_solutions.png")

- Over-regularisation apparent in reconstruction
- Generalised Tikhonov able to attenuate HF noise using derivative term
