#import "@preview/physica:0.9.6": dv
#import "@preview/physica:0.9.4": curl, grad, tensor, pdv, dv, eval, TT
#import "preamble.typ": mref, mathformatter, update-counter
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.8": *
#import "@preview/cetz:0.4.2"


//! Styling

#set math.equation(numbering: "(1)")

#show ref: it => {
  let eq = math.equation
  let head = heading

  let el = it.element
  if el != none and el.func() == eq {
    // Override equation references.
    numbering(
      el.numbering,
      ..counter(eq).at(el.location()),
    )
  } else {
    // Other references as usual.
    it
  }
}

//! Setup plugins
#show: codly-init.with()
#codly(languages: codly-languages)

//! Declarations
#let wider = h(3em)
#let vv = mathformatter(underbars: 0, bold: true, upright: true)
#let mm = mathformatter(underbars: 0, bold: true, upright: true)
#set math.mat(delim:"[")
#set math.vec(delim:"[")

#let SVD = math.op("SVD")
#let Span = math.op("Span")

//! Cover Page
#counter(page).update(0)

#set document(title: [02624 Inverse Problems and Imaging — Homework 2])
#align(center, title())

Jeppe Klitgaard `<s250250@dtu.dk>`\
Date: 2025-10-24

#outline()

= Introduction <sec:intro>

This report presents solutions to Homework 2 for the course 02624 Inverse Problems and Imaging.

GitHub Copilot was used for tab completions in the development of the code associated with this report.
Google Gemini 2.5 Pro was used for proof-reading and improvement suggestions on the final draft of the report. The prompt input and output is included in @app:gemini_review. Only a few of the suggestions were adopted, primarily related to grammar or repeated words.

The code used to generate experiments and figures can be found in the appendix.

#pagebreak()
#counter(page).update(1)
#set page(numbering: "1/1")

#counter(heading).update(0)
#set heading(numbering: "1.a.i)")

//! Content
//* Content Setup
#set page(
  margin: (
    top: 1.5cm,
    bottom: 1.0cm,
    left: 1.5cm,
    right: 2.5cm,

  )
)

#set par(spacing: 0.8em, leading: 0.55em)
#show math.equation.where(block: true): set block(spacing: 0.6em)

#set figure(gap: 0pt)
#show figure.caption: c => {
  set text(size: 10pt)
  c
}

#show heading: it => {
  set text(size: 13pt)
  set block(above: 0.8em)
  it
}
#set text(
  size: 11pt
)

= 1D Heat Problem <sec:1>

We consider the 1D backwards heat problem on the interval $[0, π]$ over time period $t ∈ (0, T≔0.1]$
modelled by the operator $K$ defined such that $y(t) = K x(t)$, corresponding to the problem:
$
  pdv(u, t) - pdv(u, x, 2) = 0, wider 0 < x < π, quad t ∈ (0, T]
$
$
  u(0, t) = u(π, t) = 0, wider u(x, 0) = f(x) \
$

The operator $K$ may be discretised on equidistant grids of $N=128$ points in both space and time to obtain a matrix $A ∈ ℝ^(N × N)$ such that $y = A x$ approximates the continuous case.

The problem is considered in a Bayesian framework with Gaussian noise as given by:
$
  Y = A x + E
$

Where $x$ is a realisation of the prior distribution $X ∼ cal(N)(0, γ^2 I)$ with $γ ≔ 10^(-1)$ and $E$
is the Gaussian noise distributed as $E ∼ cal(N)(0, σ^2 I)$ where $σ ≔ 10^(-3)$. As such, we find that observations are determined by the realisation of the noise, $e$, as:
$
  y_"obs" = y + e
$

Thus, we obtain probability densities:
$
  π_"pr" (x) &∝ exp(-x^2/(2 γ^2)) \
  π(y_"obs"|x) &∝ exp(-(y_"obs" - A x)^2/(2 σ^2)) \
  π(y_"obs") &= ∫ π(y_"obs"|x) π_"pr" (x) dif x ∝ 1 \
$

From which we may obtain the posterior distribution through Bayes' theorem:
$
  π_"post"
  &= π(x|y_"obs") = (π(y_"obs"|x) π_"pr" (x))/(π(y_"obs"))\
  &∝ exp(-1/2 (norm(y_"obs" - A x)^2_2/σ^2 + norm(x)^2_2/γ^2)) \
  &∝ exp(-1/2 (x - dash(x))^TT Γ^(-1)_"post" (x - dash(x))) \
$ <eq:posterior_gaussian_proposal>

Where the last relation proposes that the posterior is also Gaussian with mean $dash(x)$ and covariance $Γ_"post"$.

== Deriving the Ordinary Tikhonov Solution <sec:1_a>
We may verify that @eq:posterior_gaussian_proposal holds and in doing so obtain expressions for $dash(x)$ and $Γ_"post"$ by the following treatment:
$
  exp(-1/2 (x - dash(x))^TT Γ^(-1)_"post" (x - dash(x)))
  &∝ exp(-1/2 (norm(y_"obs" - A x)^2_2/σ^2 + norm(x)^2_2/γ^2))\
  &= exp(-1/2 [((y_"obs" - A x)^TT (y_"obs" - A x))/σ^2 + (x^TT x )/γ^2])\
  &= exp(-1/2 [1/σ^2 (y_"obs"^TT y_ "obs" - y_"obs"^TT A x - x^TT A^TT y_"obs" + x^TT A^TT A x) + (x^TT x )/γ^2])\
  &= exp( -1/2 [underbrace(x^TT ((A^TT A) / σ^2 + I/γ^2) x, "Quadratic") - underbrace(2 x^TT A^TT / σ^2  y_"obs", "Linear") + underbrace((y_"obs"^TT y_"obs")/σ^2, "Constant")]) \
$

Which by comparison against the standard form of the Gaussian distribution given in @eq:posterior_gaussian_proposal reveals:
$
  Γ^(-1)_"post" &= (A^TT A) / σ^2 + I/γ^2 \
  dash(x) &= Γ_"post" A^TT y_"obs"/σ^2 \
$

We may now invert $Γ^(-1)_"post"$ to obtain $Γ_"post"$ using the Woodbury matrix identity~@matrix_cookbook[Eq. 159] to find:
$
  Γ_"post"
  &= γ^2 I - γ^2 A^TT/σ (I + A/σ γ^2 A^TT/σ)^(-1) A/σ γ^2 \
  &= γ^2 I - γ^4 A^TT (σ^2 + γ^2 A A^TT)^(-1) A\
$ <eq:posterior_covariance_derived>
Using which we obtain the final expression for the posterior mean as:
$
  dash(x)
  &= γ^2 A^TT/σ^2 y_"obs" - γ^4 A^TT (σ^2 + γ^2 A A^TT)^(-1) A A^TT/σ^2 y_"obs" \
  &= γ^2/σ^2 A^TT (I - (σ^2 + γ^2 A A^TT)^(-1) γ^2 A A^TT) y_"obs" \
  &= γ^2/σ^2 A^TT (I - (σ^2 + γ^2 A A^TT)^(-1) ((σ^2 + γ^2 A A^TT ) - σ^2 I)) y_"obs" \
  &= γ^2/σ^2 A^TT (I - I + σ^2 (σ^2 I + γ^2 A A^TT)^(-1)) y_"obs" \
  &= γ^2 A^TT (σ^2 I + γ^2 A A^TT)^(-1) y_"obs" \
$ <eq:sample_mean_derived>

We seek to put this into the standard Tikhonov form~@kirsch[Eq. 2.18] for a solution $dash(x) = R_α y_"obs"$ for the discrete operator $A$:
$
  R_α ≔ (α I + A^TT A)^(-1) A^TT \
$

We do this by manipulating @eq:sample_mean_derived and introducing $α ≔ σ^2\/γ^2$:
$
  dash(x) &= γ^2 A^TT (σ^2 I + γ^2 A A^TT)^(-1) y_"obs" \
  dash(x) &= A^TT (α I + A A^TT)^(-1) y_"obs" \
  (α I + A^TT A) dash(x) &= (α I + A^TT A) A^TT (α I + A A^TT)^(-1) y_"obs" \
  (α I + A^TT A) dash(x) &= A^TT (α I + A A^TT) (α I + A A^TT)^(-1) y_"obs" \
  (α I + A^TT A) dash(x) &= A^TT y_"obs" \
  dash(x) &= (α I + A^TT A)^(-1) A^TT y_"obs" \
$

As such, we have recovered the standard Tikhonov solution form as desired.

Lastly, we seek to show that the posterior covariance, $Γ_"post"$, satisfies $Γ_"post" < γ^2 I$, which by inspection of @eq:posterior_covariance_derived is equivalent to showing that the second term is positive definite, which we also note to be symmetric by inspection:
$
  cancel(γ^4) A^TT (σ^2 + γ^2 A A^TT)^(-1) A ≻ 0 \
$

Where $γ^4 > 0$ by construction and may be set to unity without loss of generality. We leverage the quadratic form definition of positive definiteness~@matrix_cookbook[Eq. 439] to show that for any non-zero vector $y ∈ ℝ^N$:
$
  y^TT A^TT (σ^2 + γ^2 A A^TT)^(-1) A y= z^TT B^(-1) z
$

Where we have defined $z ≔ A y$ and $B ≔ (σ^2 + γ^2 A A^TT)$ and seek to leverage the property that the inverse of a matrix $M$ is positive definite if  $M$ itself is positive definite~@matrix_cookbook[§9.6.4]:
$
  z^TT B z = σ^2 z^TT z + γ^2 z^TT A A^TT z = σ^2 norm(z)^2_2 + γ^2 norm(A^TT z)^2_2 > 0
$

As such, we find that $B$ is positive definite, thus implying $B^(-1)$ to also be positive definite, which in turn requires that $Γ_"post" < γ^2 I$ by the construction above.

This result is consistent with our intuition that incorporation of observations into the posterior distribution should only reduce uncertainty in the estimate when compared against the prior distribution.

== Analytical and Sampled Bayesian Estimation of Initial Temperature <sec:1_b>

We proceed to probe the initial temperature estimation of the problem described above on an equidistant grid of $N=128$ points over the domain $s ∈ [0, π]$ by drawing $N_"samples"$ samples from the posterior distribution $π_"post" (x)$ given in @eq:posterior_gaussian_proposal from which the sample mean and variances are computed across the domain.

Additionally, the posterior mean and variances may be computed directly using the analytical expressions found in @eq:posterior_covariance_derived and @eq:sample_mean_derived.

The results of both approaches are shown in @fig:1_b, in which we observe close agreement between the analytical and sampled approaches, with both approximating the true initial temperature profile. The effect of finite sampling is particularly visible on the credibility interval of the sampled mean as observed in the left panel of the figure, which may be improved by increasing $N_"samples"$ as seen in the right panel.

#figure(
  image("output/1_b.png"),
  caption: [
    Analytical and sampled Bayesian estimation of the initial temperature $x(s)$ from noisy observations $y_"obs"$ at time $T=0.1$ with noise level $σ = 10^(-3)$ and prior standard deviation $γ = 10^(-1)$. The black line indicates the true initial temperature.
  ]
) <fig:1_b>

== Reduced Observational Region <sec:1_c>

We now consider the case in which only the interval $(1/4 π, 3/4 π)$ of the spatial domain is observed, which limits the amount of information available for the inversion. This is done by only including the rows of $A$ and entries of $y_"obs"$ corresponding to this region when computing the posterior mean and covariance.

Given our heat problem is defined with unit diffusivity, the terminal time $T$ alone determines the degree of ill-posedness of the problem, with shorter terminal times ordinarily resulting in a less ill-posed problem.
However, by limiting the observation region, a longer terminal time allows more information to propagate into the observed region from the unobserved regions, which may improve the inversion results in the unobserved regions.

For this reason, we investigate the inversion for a longer terminal time of $T=0.1$ as compared to the shorter time of $T=0.005$ as seen in @fig:1_c.

We may compare the left panel of @fig:1_c against @fig:1_b to observe the degradation in inversion quality due to the limited observation region, where notably the uncertainty in the unobserved regions increases and the estimate is shifted towards the prior mean. This effect is more pronounced in the right panel corresponding to a shorter terminal time of $T=0.005$, in which the uncertainty away from the observations approaches the level of the prior distribution and the estimate reverts to that of the prior.

#figure(
  image("output/1_c_double.png"),
  caption: [
    Analytical and sampled Bayesian estimation of the initial temperature $x(s)$ from noisy observations $y_"obs"$ at times $T=0.1$ (left) and $T=0.005$ (right) with noise level $σ = 10^(-3)$ and prior standard deviation $γ = 10^(-1)$. Only the region $(1/4 π, 3/4 π)$ is observed. The black line indicates the true initial temperature.
  ]
) <fig:1_c>

= Abel Problem and $H^1$-Regularisation <sec:2>

We consider the Abel problem corresponding to the integral operator $K$ defined such that $y(t) = K x(t)$, where:
$
  K x(t) = ∫_0^1 k(t, s) x(s) dif s,
  wider
  k(t, s) = cases(
    1/sqrt(t-s)\, &quad s < t,
    0\, &quad s > t
    ),
  wider
  x(s) = s^(3\/4) - 1
$

We again discretise the operator $K$ by a matrix $A ∈ ℝ^(N × N)$ where $N=1024$. We introduce a Gaussian measurement noise of 10% relative to the norm of the noiseless measurements to produce noisy observations $y^δ$ from $y = A x$. Lastly, we let $D$ denote a finite difference matrix such that $D x$ approximates the first derivative of $x$.

== Show that $x(s) ∈ H^1(0, 1)$ <sec:2_a>

In order to show that $x(s) = s^{3/4} - 1$ lies in the Sobolev space $H^1(0, 1)$, we need to demonstrate that both $x(s)$ and its weak derivative $x'(s)$ are square-integrable over the interval. We first consider $x(s)$:
$
  ∫_0^1 |x(s)|^2 dif s
  = ∫_0^1 (s^(3/4) - 1)^2 dif s
  = ∫_0^1 (s^(3/2) - 2 s^(3/4) + 1) dif s\
  = [ (2/5) s^(5/2) - 2 ⋅ 4/7 s^(7/4) + s]_0^1
  = 14/35 - 40/35 + 35/35
  = 9/35 < ∞ \
$

Thus proving that $x(s)$ is square-integrable. Next, we compute the weak derivative, here equivalent to the classical derivative due to the smoothness of $x(s)$:
$
  x'(s) = 3/4 s^(-1/4) \
$
Which we show to be square-integrable as well:
$
  ∫_0^1 |x'(s)|^2 dif s
  = ∫_0^1 (3/4 s^(-1/4))^2 dif s
  = 9/16 ∫_0^1 s^(-1/2) dif s
  = 9/16 [2 s^(1/2)]_0^1
  = 9/8 < ∞ \
$

Consequently, we conclude that $x(s) ∈ H^1 (0, 1)$ as implied by $x(s) ∈ L^2 (0, 1)$ and $x'(s) ∈ L^2 (0, 1)$.

== Uniqueness and Existence of Generalised Tikhonov Solution

By approximating the $H^1$ norm by $norm(x)_2^2 + norm(D x)_2^2$, we obtain the following optimisation problem for the generalised Tikhonov solution:
$
  x_α = limits("arg min")_(x∈H^1 (0, 1)) med norm(A x - y^δ)_2^2 + α (norm(x)_2^2 + norm(D x)_2^2)
$ <eq:generalised_tikhonov_problem>

We let the objective function be denoted $f(x)$ and seek to prove that a unique solution $x_α$ exists for any $α > 0$.

To prove existence, we seek to use Weierstrass Extreme Value Theorem~@dtu_02611_L1[S. 16], which states that a continuous function defined on a compact set attains a global minimum on that set. Given that the objective function is continuous by inspection, we only need to prove that the sublevel sets $S_c = { x ∈ ℝ^n : f(x) ≤ c }$ are bounded, which is implied by the coercivity of the objective function:
$
  norm(A x - y^δ)_2^2 + α (norm(x)_2^2 + norm(D x)_2^2) → +∞ "as" norm(x)_2 → +∞ \
$

Where the derivative and residual terms are non-negative by construction, thus reducing the problem to showing that $α norm(x)_2^2$ is coercive, which is trivially true for any $α > 0$. This allows us to restrict our search for a minimiser to a sublevel set which is compact, thus allowing us to apply the Extreme Value Theorem to conclude that a global minimiser exists on an arbitrarily chosen sublevel set $S_c$.

Uniqueness may be determined by considering the convexity of the objective function. Any norm is convex by construction, and raising to a power $p >= 1$ preserves convexity, thus implying that the squared norms are convex.

Non-negative scaling likewise preserves convexity, thus preserving the convexity of the second term $α (norm(x)_2^2 + norm(D x)_2^2)$ for $α>0$. Lastly, the sum of convex functions is itself convex, thus implying that the entire objective function is convex. To guarantee uniqueness, we require _strict_ convexity, which may be guaranteed by considering the Hessian of $norm(x)_2^2$:
$
  ∇^2 norm(x)_2^2 = 2 I ≻ 0\
$
Thus implying that $norm(x)_2^2$ is strictly convex by the second order condition.

This approach relies heavily on the known convexity properties of various functions and operations as taught by the 02611 Optimisation for Data Science course. As an alternative, we may rely on the second order condition to prove strict convexity of the entire objective function by considering its Hessian:
$
  ∇^2 [norm(A x - y^δ)_2^2 + α (norm(x)_2^2 + norm(D x)_2^2)]
  = 2 A^TT A + 2 α I + 2 α D^TT D \
$
Where $A^TT A$ and $D^TT D$ are both positive semi-definite by the same treatment used in @sec:1_a:
$
  y^TT B^TT B y = norm(B y)^2_2 >= 0 wide ⇔ wide B^TT B ≽ 0
$
Where $y∈ℝ^N$ and $B$ is either $A$ or $D$, thus implying:
$
  ∇^2 [norm(A x - y^δ)_2^2 + α (norm(x)_2^2 + norm(D x)_2^2)]
  = 2 A^TT A + 2 α I + 2 α D^TT D ≻ 0
$
As such, we conclude that the objective function is strictly convex and thus admits a unique minimiser for any $α > 0$.

#context [
  #let previous-counter = counter(heading).get()
  #update-counter(counter(heading), (none, 0))
  #set heading(numbering: "1.*)")
  == Closed Form Solution for Generalised Tikhonov <sec:closed_form_tikhonov>
  #counter(heading).update(previous-counter)
  #set heading(numbering: "1.a)")
]

In order to compute the Tikhonov solutions for the subsequent sections, we will first derive
a closed form solution to the generalised Tikhonov problem from @eq:generalised_tikhonov_problem, to which end we compute the gradient and set it to zero to find the minimiser:
$
  0
  &= ∇ (norm(A x - y^δ)_2^2 + α (norm(x)_2^2 + norm(D x)_2^2))\
  &= ∇ ( (A x - y^δ)^TT (A x - y^δ) + α (x^TT x + x^TT D^TT D x) )\
  &= 2 A^TT A x - 2A^TT y^δ + 2 α x + 2α D^TT D x\
  &= (A^TT A +α I + α D^TT D) x - A^TT y^δ \
$
Which we may solve for $x$ to obtain the closed form solution for the generalised Tikhonov problem:
$
  x_α = (A^TT A + α(I + D^TT D))^(-1) A^TT y^δ \
$ <eq:generalised_tikhonov_solution>

Where the ordinary Tikhonov solution may be obtained by setting $D=0$.
$
  x_α = (A^TT A + α I)^(-1) A^TT y^δ \
$ <eq:ordinary_tikhonov_solution>

These closed-form solutions are used in the computational experiments of the subsequent sections.

== Ordinary Tikhonov and Discrepancy Principle <sec:2_c>

We seek to find an optimal regularisation parameter $α > 0$ by applying the Morozov discrepancy principle to the ordinary Tikhonov solution, @eq:ordinary_tikhonov_solution. Additionally, we assume that an estimate of the noise level $δ$ is available, in which case the discrepancy principle requires that we find an $α$ such that~@hansen[§5.2]:
$
norm(K x_α - y^δ) ≈ ν δ
$

Where $ν >= 1$, here chosen as $ν=1.05$, is a safety factor chosen to produce a conservative estimate of the regularisation parameter.
We find such an $α_"DP"$ by performing a search across $200$ linearly spaced values of the regularisation parameter $α∈[0, 0.2]$ and selecting the value that most closely satisfies the discrepancy principle to obtain @fig:2_c, in which we also display the error norm $norm(x_α - x)_2$ to gauge the performance of the discrepancy principle.

We find that the discrepancy principle provides an estimated regularisation parameter very close the optimal value. Additionally, we observe in the right panel that the error does not depend strongly on the choice of $α$ within the neighbourhood of the optimal value, thus making the discrepancy principle a robust method for choosing the regularisation parameter in this case.

#figure(
  image(
    "output/2_c.png",
  ),
  caption: [
    Application of the Morozov discrepancy principle to determine an optimal regularisation parameter $α$ for the ordinary Tikhonov solution of the Abel problem. The left panel contains the residual norm $norm(A x_α - y^δ)_2$ while the right panel contains the error norm $norm(x_α - x)_2$. The horizontal dashed line indicates the target discrepancy level $ν δ$.
  ]
) <fig:2_c>

== Generalised Tikhonov and Discrepancy Principle <sec:2_d>

We repeat the procedure of @sec:2_c for the generalised Tikhonov solution, where we let $D$ denote a central difference first order differentiation matrix such that $D x$ approximates the first derivative of $x$. At the boundaries, the first order forward and backward difference schemes are used instead:
$
  D = 1/h mat(
    -1, 1;
    -1/2, 0, 1/2;
    , ⋱, ⋱, ⋱;
    ,, -1/2, 0, 1/2;
    ,,, -1, 1;
  ),
wider h = 1/N
$

We again perform the search undertaken in @sec:2_d, now using the generalised Tikhonov solution from @eq:generalised_tikhonov_solution to obtain @fig:2_d, where notably the search interval has been changed to $α∈[0, 0.015]$ to better capture the behaviour of the generalised solution. When compared against the ordinary Tikhonov solution of @sec:2_c, we find that a lower degree of regularisation is required on the norm of $x$ when the norm of its derivative is also penalised, which is consistent with our intuition given the high degree of smoothness of the true solution, which is shown alongside the solutions below in @fig:2_d_solutions.

Notably, the generalised Tikhonov solution is able to provide a better reconstruction as evidenced by the lower error norms observed in the right panel of @fig:2_d when compared against @fig:2_c.
Interestingly, the introduction of the safety factor $ν=1.05$ appears to negatively impact the performance of the discrepancy principle in this case, as can be seen in the left panel of @fig:2_d, where the unaugmented noise level matches the optimal solution almost exactly.

Were we to remove the safety factor for both the ordinary and generalised Tikhonov solutions however, we observe from @fig:2_c that the ordinary Tikhonov solution would suffer from under-regularisation as a result. As such, the choice of safety factor appears to require a heuristic or intuition-based approach in and of itself, though notably the impact of the choice of $ν$ will not be known without access to the true solution.

#figure(
  image(
    "output/2_d.png",
  ),
  caption: [
    Morozov discrepancy principle applied to the generalised Tikhonov solution of the Abel problem.
  ]
) <fig:2_d>

Lastly, we compare the reconstructions obtained by both the ordinary and generalised Tikhonov solutions using the regularisation parameters obtained by the discrepancy principle in @fig:2_d_solutions, which reveals that the generalised solution is able to capture the smoothness of the underlying problem to significantly reduce the noise in the reconstruction when compared against the ordinary solution. Interestingly, both solutions appear biased at the left side of the domain, whereas only the generalised solution appears biased at the right side.

It would be interesting to further investigate how the performance of a generalised Tikhonov solution could be improved by introducing separate regularisation parameters for the norm of $x$ and its derivative.

#figure(
  image(
    "output/2_d_solutions.png",
  ),
  caption: [
    Comparison of the ordinary and generalised Tikhonov solutions of the Abel problem using regularisation parameters obtained by the Morozov discrepancy principle.
  ]
) <fig:2_d_solutions>

// Appendix
#set heading(numbering: "A.1)", supplement: [appendix])
#set math.equation(numbering: (c) => [(A.#(c))])
#counter(heading).update(0)

#pagebreak()

= Appendix

== Code preamble <app:code_preamble>
#raw(read("output/preamble.py"), block: true, lang: "python")

== Code for @sec:1_b <app:code_1_b>
#raw(read("output/1_b.py"), block: true, lang: "python")

== Code for @sec:1_c <app:code_1_c>
#raw(read("output/1_c.py"), block: true, lang: "python")

== Code for @sec:2_c <app:code_2_c>
#raw(read("output/2_c.py"), block: true, lang: "python")

== Code for @sec:2_d <app:code_2_d>
#raw(read("output/2_d.py"), block: true, lang: "python")

== Google Gemini 2.5 Pro Review Comments <app:gemini_review>

Note: This uses a generic _reviewer_ Gem.

#{
  for i in range(2, 9) {
    image("gemini_feedback.pdf", page: i)
  }
}

#pagebreak()
#bibliography("refs.bib")
