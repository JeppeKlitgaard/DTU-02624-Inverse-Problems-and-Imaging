#import "@preview/physica:0.9.6": dv
#import "@preview/physica:0.9.4": curl, grad, tensor, pdv, dv, eval, TT
#import "preamble.typ": mref, mathformatter
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

#align(center, text(17pt)[
  *02624 Inverse Problems and Imaging — Homework 1*
])

Jeppe Klitgaard `<s250250@dtu.dk>`\
Date: 2025-10-24

#outline()

= Introduction <sec:intro>

This report presents solutions to Homework 1 for the course 02624 Inverse Problems and Imaging.

As discussed in class with Professor Kim Knudsen, the report spans slightly more than the initially
required 7 pages in order to improve readability.

Code and larger versions of the figures can be found in the appendix alongside optional tasks and
lengthier intermediate mathematical derivations.

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

= Backwards Heat Equation <sec:1>
We are given the _forward model_ of the _backwards heat equation_:
$
K x(t) &= ∫_0^π x(s) k(t, s) dif s quad "where" quad k(t, s) &= 2/π ∑_(n=1)^∞ e^(-T n^2) sin(n t) sin(n s)
$ <eq:1>

== Singular System <sec:1.a>
Using that $K$ is compact in $L^2(0, π)$, we can demonstrate that it has the singular system $(x_j, y_j, μ_j)$ where the _singular functions_, $x_j, y_j$ are given by:
$
  x_j (t) &= sqrt(2/π) sin(j t) wider y_j (s) &= sqrt(2/π) sin(j s)
$ <eq:singular_functions>
And the singular values, $μ_j$ are given by:
$
  μ_j &= e^(-T j^2)
$ <eq:singular_values>

To verify that this is indeed the singular system for $K$, we recall the definition of a singular system~@kirsch[Theorem A.57], for which it holds for every $j ∈ J ⊆ ℕ_1$:
$
K x_j &= μ_j y_j wider
K^* y_j &= μ_j x_j
$ <eq:singular_system_def>

We first verify that $K x_j = μ_j y_j$:
$
K x_j (t)
&= ∫_0^π x_j (s) k(t, s) dif s
= ∫_0^π sqrt(2/π) sin(j s) (2/π ∑_(n=1)^∞ e^(-T n^2) sin(n t) sin(n s)) dif s\
&= (2/π)^(3/2) (∑_(n=1)^∞ e^(-T n^2) sin(n t) ∫_0^π sin(j s) sin(n s) dif s)
= (2/π)^(3/2) (∑_(n=1)^∞ e^(-T n^2) sin(n t) π/2 δ_(j n))\
&= e^(-T j^2) sqrt(2/π) sin(j t)
= μ_j y_j (t)
$

Where the summation and integral interchange was enabled by uniform convergence of the series through inspection of @eq:1, in which the Gaussian decay factor $e^(-T n^2)$ ensures rapid convergence. The orthogonality of the sine function in $L^2(0, π)$ is shown in Appendix @app:sine_ortho, where we find that the
integral over the sines of different integer frequencies $n, j$ is given by $π/2 δ_(j n)$ as shown in @eq:sine_ortho_result.

To consider the other equation in @eq:singular_system_def, we must first determine the adjoint operator, $K^*$, which is obtained through the conjugate-transpose of the kernel.
$
  k^*(s, t) = k(overline(t), s)
$
By noting that $k(s, t)$ is real-valued and symmetric by inspection of @eq:1, we have:
$
  k^*(s, t) = k(t, s)
$

Such that:
$
K^* y(s)
&= ∫_0^π y(t) k^*(s, t) dif t = ∫_0^π y(t) k(t, s) dif t =K y(s) &
$

From which we may conclude that $K$ is _self-adjoint_, $K = K^*$.

As such, we may repeat treatment above for the equation $K^* y_j = μ_j x_j$, which becomes true by inspection.

In addition to this, the singular functions must be orthonormal in $L^2(0, π)$, which they are by the treatment in @app:sine_ortho and the singular values must be non-negative, which is indeed the case for the given system.

Lastly, we note that the singular values decay to zero as $j → ∞$ since $T > 0$, as we would expect for a compact operator.

We conclude that $K$ is indeed a compact operator with singular system $(x_j, y_j, μ_j)$ as given in @eq:singular_functions and @eq:singular_values.

== Discretisation of $K$ <sec:1.b>

We seek to discretise $K$ by a matrix $A ∈ ℝ^(n × n)$ where $n=128$ and $T=0.03$, to which end the spatial and temporal domains $[0, π]$ are discretised into $n$ equidistant points:
$
  t_i &= (k π)/(n+1) wider s_i &= (k π)/(n+1) wider "where" wider k=1, …, n
$
Where we note that the ends of the interval are not included as the sine functions are zero at these points, which would not result in an orthogonal basis.

Next, we construct a matrices $U, V$ whose columns are the left and right singular vectors, $x_j, y_j$, sampled at the points $s_i$:
$
  U_(i j) &= sqrt(2/(n+1)) sin(j ⋅ t_i)
  wider
  V_(i j) &= sqrt(2/(n+1)) sin(j ⋅ s_i)
$
Where the scaling factor ensures that the columns are normalised.

From this, we may construct the diagonal matrix $Σ$ whose entries are the singular values $μ_j$:
$
  Σ_(j j) = μ_j = e^(-T j^2)
$

Which enables us to construct the discretisation of $K$ directly:
$
  A = U Σ V^TT
$

Here, we have exploited that we have analytical expressions for the components of the singular system. In a more general setting, we may have to resort to approximation of the integrals in @eq:1 using a quadrature rule.

We can compare the singular values of $A$ against those for $K$ obtained by analytical means, as seen in @fig:1_b_singular_values. We observe that they decay rapidly at the same rate up until the 33rd singular value, after which the singular values of $A$ remain constant at a value
determined by the numerical precision of the computation, which is approximately $10^(-16)$ in this case, as marked by the dotted black line.

#figure(
  image("output/1_b_singular_values.png"),
  caption: [
    Singular values of the discretised operator $A$ for $n=128$ and $T=0.03$ compared against the analytical expression for the singular values of $K$.
  ]
) <fig:1_b_singular_values>

Similarly, we can compare the singular vectors of $A$ against the singular functions of $K$, as seen in @fig:1_b_singular_vectors. We observe an excellent match between the discretised vectors and the continuous singular functions at lower _modes_. At higher frequency modes, namely those past the $j=33$ threshold, we observe that the limited machine precision affects the accuracy of the representation and noise dominates the decomposition.

#figure(
  image("output/1_b_singular_vectors.png"),
  caption: [
    3 select singular vectors and their associated singular functions for $n=128$ and $T=0.03$.
  ]
) <fig:1_b_singular_vectors>

== Inversion using Truncated SVD <sec:1.c>

Given the function:
$
  x(s) = s sin(s) exp((s - π\/2)^2)
$ <eq:1_c>

We seek to recover a solution to the inverse problem $K x = y$ given noisy measurements $y^δ$. We construct our measurements by computing the forward model $y = K x$ and adding Gaussian noise with relative size $norm(n)_(L^2(0, π)) \/ norm(y)_(L^2(0, π)) = 2%$.

A naïve inversion of the problem using the inverse of $A$ is expected to perform poorly due to the hegemony of noise at higher frequency modes. We may instead use a _truncated SVD_ approach in which only the first $j^⋆$ singular values are included in the inversion.

Finding a suitable choice of the threshold index $j^⋆$ is undertaken through the construction of a _Picard plot_ containing the singular values $μ_j$, the expansion coefficients $abs((y^δ, y_j))$ and the solution coefficients $abs((y^δ, y_j))\/μ_j$ as seen in @fig:1_c_picard_plot.

#figure(
  image("output/1_c_picard_plot.png"),
  caption: [
    Picard plot of the noisy data $y^δ$ with $2%$ relative noise.
  ]
) <fig:1_c_picard_plot>

We recall the _discrete Picard condition_~@hansen[3.16], which states that the expansion coefficients must decay faster than the singular values for a stable solution to exist. We use this in conjunction with @fig:1_c_picard_plot to choose a threshold index $j^⋆ = 9$ at which the solution coefficients begin to diverge, indicating that noise is dominating the solution. We note that this choice is somewhat subjective and given we only have a single realisation of the noise available in the plot, we conservatively place it at $j^⋆ = 9$ as opposed to $j^⋆ = 10$.

With this, we are able to compute a reconstruction of $x$ using the truncated SVD:
$
  U Σ V^† = A
$
From which we may construct the truncated matrices by reducing the number of columns of $U, V$ and the size of $Σ$ to $j^⋆$.

From this, we may reconstruct $x$ using @eq:singular_system_def:
$
  x_"rec" = ∑_(j=1)^(j^⋆) ((y^δ, y_j))/μ_j x_j = underbrace(tilde(V) tilde(Σ)^(-1) tilde(U)^†, tilde(A)^(-1)) y^δ
$ <eq:1_c_tsvd>

Where the tildes denote that the matrices are truncated to the first $j^⋆$ columns/rows.

The resulting reconstruction is shown in @fig:1_c_reconstruction, where we observe a good match between the true model $x$ and the reconstruction $x_"rec"$. The noisy measurement $y^δ$ is also shown alongside the true signal $y$.

We observe that the regularised inversion is able to recover the main features of $x$ despite the presence of significant noise in the measurements. Attempting to invert the problem without regularisation results in a completely unusable solution dominated by noise, as seen by $x_"naive"$ in @fig:1_c_reconstruction, whose values are to be read off the right vertical axis.

#figure(
  image("output/1_c_reconstruction.png"),
  caption: [
    Reconstruction of $x$ using truncated SVD with $j^⋆ = 9$ compared against the true $x$.
  ]
) <fig:1_c_reconstruction>

= Abel Problem <sec:2>

We are given the _forward model_ of the _Abel problem_ on the domain $L^2(0, 1)$:
$
K x(t) &= ∫_0^t 1/sqrt(t-s) x(s) dif s
$ <eq:2_abel_forward>

#heading([Proving Compactness of $K$], level: 2, numbering: (..nums) => [#nums.at(0).\*)])
#counter(heading).update((2, 0))

Due to space constraints, this section has been relegated to the appendix in @app:2.optional

== Computing the Forward Model and Noisy Measurements <sec:2.a>
We seek to compute the forward model $y = K x$ associated with the Abel problem given in @eq:2_abel_forward for the function:
$
  x(s) = s^2
$

We construct noisy measurements by addition of Gaussian noise with relative size $norm(n)_(L^2(0, 1)) \/ norm(y)_(L^2(0, 1)) = 5%$ such that:
$
  y^δ = y + n
$

Which produces the realisation seen in @fig:2_a.
#figure(
  image("output/2_a.png"),
  caption: [
     Noisy realisation of Abel problem with $x(s)=s^2$ perturbed by $5%$ relative white noise.
  ]
) <fig:2_a>

The code for this section can be found in @app:code_2_a.

== Truncated SVD Reconstruction <sec:2.b>
In similar fashion to @sec:1.c, we seek to reconstruct $x$ from the noisy measurements $y^δ$ using a truncated SVD approach already described in @sec:1.c.

As in @sec:1.c, we create a Picard plot as seen in @fig:2_b_picard in order to gauge the information carried by the components of the SVD.

#figure(
  image("output/2_b_picard.png"),
  caption: [
    Picard plot of the noisy data $y^δ$ with $5%$ relative noise.
  ]
) <fig:2_b_picard>

Qualitatively, we observe that components past $k≈20$ are dominated by the noise,
which helps guide our selection of truncation levels $k$ to consider for our reconstructions seen in @fig:2_b_reconstructions.

#figure(
  image("output/2_b_reconstructions.png"),
  caption: [
    Reconstructions of $x(t)$ using truncated SVD with varying truncation levels $k$ compared against the true $x(t)$. Second panel shows the relative error given by $(x - x_k)\/norm(x)_2$.
  ]
) <fig:2_b_reconstructions>

Figures @fig:2_b_picard[] and @fig:2_b_reconstructions[] both indicate that the problem under consideration is much less ill-posed than that of the backwards heat equation in @sec:1, as we are able to perform a stable, albeit noisy, inversion even at high truncation levels such as $k=100$.

Nevertheless, the analysis of the Picard plot is in agreement with the results of @fig:2_b_reconstructions, which both indicate that truncation levels around $k=20$ provide a good balance between variance and bias in the reconstruction.

We notice by inspection of the relative errors in the second panel of @fig:2_b_reconstructions that
a significant bias is introduced at the boundary $s=1$.

=== Reconstruction Bias at Boundary $s=1$ <sec:2.b.i>

Having observed that the TVSD reconstructions drop to zero at $s=1$ for all truncation levels, we investigate further by considering the Truncated SVD reconstruction as given in @eq:1_c_tsvd:
$
x_"rec" (s) = ∑_(j=1)^(k) ((y^δ, y_j))/μ_j x_j (s)
$ <eq:2_b_tsvd_reconstruction>

Since neither $(y^δ, y_j)$ nor $μ_j$ depend on $s$ nor are generally zero, we focus on the singular functions $x_j (s)$, which are related to the forward operator $K$ through the adjoint operator $K^*$ as per @eq:singular_system_def:
$
  (K^* K x_j)(s) = μ_j^2 x_j (s)
$

Which requires us to find an expression for the adjoint operator $K^*$, defined through the inner product:
$
  (K x, y) &= (x, K^* y)\
  &= ∫_0^1 (∫_0^t 1/sqrt(t-s) x(s) dif s) overline(y(t)) dif t\
  &= ∫_0^1 x(s) (∫_s^1 1/sqrt(t-s) overline(y(t)) dif t) dif s = (x, K^* y)
$
Where the interchange of the integrals is enabled by Fubini's theorem and we identify the adjoint operator as:
$
  K^* y(s) = ∫_s^1 1/sqrt(t-s) y(t) dif t
$

The transformation of the limits may be appreciated by considering the sketch below:

#figure(
  cetz.canvas({
    import cetz.draw: *
    let x-range = 1.5
    let y-range = 1.5
    let arrow-style = (mark: (end: "stealth", fill: black))
    let integral-style = (stroke: (dash: "dashed"))

    let x-offset = 0
    let x-end = (x-range +0.5 + x-offset, 0)
    let y-end = (0 + x-offset, y-range +0.5)

    line((0 + x-offset, 0), x-end, ..arrow-style)
    line((0 + x-offset, 0), y-end, ..arrow-style)
    line((0 + x-offset, 0), (x-range + x-offset, y-range))
    line((0 + x-offset, y-range), (x-range + x-offset, x-range))

    // Axis Labels
    content(x-end, $s$, anchor: "south-west")
    content(y-end, $t$, anchor: "south-west")

    for idx in range(1, 5) {
      let t = x-range * idx / 5
      line((0 + x-offset, t), (t + x-offset, t), ..integral-style)
    }

    let x-offset = 4
    let x-end = (x-range +0.5 + x-offset, 0)
    let y-end = (0 + x-offset, y-range +0.5)
    line((0 + x-offset, 0), x-end, ..arrow-style)
    line((0 + x-offset, 0), y-end, ..arrow-style)
    line((0 + x-offset, 0), (x-range + x-offset, y-range))
    line((0 + x-offset, y-range), (x-range + x-offset, x-range))

    // Axis Labels
    content(x-end, $s$, anchor: "south-west")
    content(y-end, $t$, anchor: "south-west")

    for idx in range(1, 5) {
      let s = x-range * idx / 5
      line((s + x-offset, s), (s + x-offset, y-range), ..integral-style)
    }
  })
)

This allows us to reconsider the left singular functions $x_j (s)$ at $s=1$:
$
  (K^* K x_j)(s=1) = ∫_(s=1)^1 ∫_0^t 1/sqrt(t-τ) y(τ) dif τ = 0 = μ_j^2 x_j (s=1)
$
Where the outer integrals upper and lower limits are equal, resulting in a zero value for any reconstruction evaluated at $s=1$.
Given $μ_j ≠ 0$, we find that $x_j (s=1) = 0$ for all $j$.

Returning to the reconstruction expression of @eq:2_b_tsvd_reconstruction, we find that:
$
  x_"rec" (s=1) = ∑_(j=1)^(k) ((y^δ, y_j))/μ_j x_j (s=1) = 0
$

We note that the proof does not immediately extend to the discrete case where the outer integral is
approximated by a quadrature rule, though given the discretised operator must approximate its
continuous equivalent with increasingly fine discretisations, it must hold $x_"rec" (s=1)$ is biased
towards zero as the discretisation is refined.

//  by adopting the notation that the reconstruction using truncation level $k$ is denoted $x^((k))$
// and considering this reconstruction using summation notation:
// $
//   x^((k))_i = ∑_(j=1)^k sum_(l=1)^k sum_(m=1)^n V_(i j) (Σ^(-1))_(j l) (U^TT)_(l m) y^δ_m
// $
// Noting that the diagonality of $Σ^(-1)$ implies $Σ^(-1)_(j l) = 1/μ_j δ_(j l)$ and $(U^TT)_(l m) = U_(m l)$ we simplify while restricting ourselves to the last discretisation point $i=n$:
// $
//   x^((k))_n = ∑_(j=1)^k V_(n j) 1/μ_j ∑_(m=1)^n U_(m j) y^δ_m
// $ <eq:abel_reconstruction_last_point>

// We seek to determine the value of $V_(n j)$, which we do through the discretiation of @eq:singular_system_def:
// $
//   A_(i j) V_(j k) = U_(i j) Σ_(j k) ⇔ A_(l i) U_(l j) = V_(i l) Σ_(j l) = V_(i j) μ_j
// $
// Where $A_(j i)$ are the elements of the discretised operator of @eq:2_abel_forward:
// $
//   A_(j i) = cases(
//     sqrt(h/(t_j - s_i)) quad & "if" quad s_i < t_j\
//     0 quad & "otherwise" quad
//   )
// $

// Given that $s_n > t_j$ for all $j$ we find that $A_(j n) = 0$ for all $j$, which implies:
// $
//   0 = V_(n j) μ_j
// $ <eq:V_nj_zero_constraint>

// Given that we have truncated the sum over index $j$ at a sensible level $k$, we have guaranteed that $μ_j ≠ 0$ for all $j ≤ k$, which in turn requires $V_(n j) = 0$ for all $j ≤ k$ by @eq:V_nj_zero_constraint.

// If we now use this result in @eq:abel_reconstruction_last_point we find:
// $
//   x^((k))_n = ∑_(j=1)^k 0 ⋅ 1/μ_j ∑_(m=1)^n U_(m j) y^δ_m = 0
// $

// Which proves that the reconstruction at the last discretisation point is always zero.

== Tikhonov Regularisation <sec:2.c>

We additionally consider Tikhonov regularisation parameterised by $α$ to produce reconstructions $x_α$ given noisy measurements $y^δ$. There are different equivalent formulations of Tikhonov regularisation; we consider the following closed-form formulation:
$
  x_α = tilde(R)_α y^δ
  &= (K^* K + α I) ^(-1) K^* y^δ\
  &= (A^TT A + α I) ^(-1) A^TT y^δ
$

Which is easily computed for a suitable range of parameters $α ∈ [10^(-4), 10]$
in order to obtain the reconstructions as seen in @fig:2_c_tikhonov_reconstructions.

#figure(
  image("output/2_c_tikhonov_reconstructions.png"),
  caption: [
    Reconstructions of $x(t)$ using Tikhonov regularisation with varying regularisation parameters $α$ compared against the true $x(t)$.
  ]
) <fig:2_c_tikhonov_reconstructions>

We observe that $α ≈ 10^(-2)$ appears to be a good balance between noise suppression and the bias introduced by regularisation which drives the reconstruction towards zero at $s=1$ as also observed for the TSVD regularisation in @sec:2.b.

#figure(
  image("output/2_c_tsvd_vs_tikhonov.png"),
  caption: [
    Comparison of TSVD reconstruction with $k=20$ and Tikhonov reconstruction with $α=4×10^(-2)$.
  ]
) <fig:2_c_tsvd_vs_tikhonov>

By comparing the two types of reconstruction in @fig:2_c_tsvd_vs_tikhonov, we observe that they have rather different characteristics. The Tikhonov reconstruction is not subject to the systematic bias introduced by the crude truncation of the higher modes in the TSVD approach, but also does not
suppress the high-frequency noise as effectively as observed in the TSVD reconstruction. As such, the choice between the two methods will depend on the
application at hand and poses the familiar bias-variance trade-off associated with regularisation.

== Regularisation of a Discontinuous Function <sec:2.d>

We attempt the above regularisation methods on a discontinuous function constructed from the original quadratic function and a scaled Heaviside step function:
$
  x(s) = s^2 + 4 ⋅ H(s - 1/2)
$
Where $H(x)$ is the Heaviside step function defined as:
$
  H(x) = cases(
    0 quad & "if" quad x < 0\
    1/2 quad & "if" quad x = 0\
    1 quad & "if" quad x > 0
  )
$

Which allows us to produce another comparison of the two regularisation methods as seen in @fig:2_d_discontinuous.

#figure(
  image("output/2_d_discontinuous.png"),
  caption: [
    Comparison of TSVD reconstructions with $k∈{20, 50}$ and Tikhonov reconstruction with $α=10^(-2)$ for a discontinuous function.
  ]
) <fig:2_d_discontinuous>

We observe that the regularisation methods struggle to accurately reconstruct the discontinuity at $s=1/2$, which is to be expected given that both methods implicitly assume some degree of smoothness in the underlying function being reconstructed. Increasing $k$ or decreasing $α$ allows for better reconstruction of the high-frequency components associated with the discontinuity at the cost of decreased noise suppression.

= Blurring Problem <sec:3>

We are given an operator $tilde(R)_α$ defined as:
$
  x_a = tilde(R)_α y = sum_j μ_j^2 / (μ_j^3 + α) (y, y_j) x_j
$

== Showing that $tilde(R)_α$ is a Regularisation Operator <sec:3.a>
We seek to show that:
$
  lim_(α→0) tilde(R)_α K x = x
$

By injectivity the singular functions $x_j, y_j$ span their respective spaces such that any $x ∈ X$ may be decomposed as:
$
  x = sum_j (x, x_j) x_j
$
Similarly, using $K x_j = μ_j y_j$ we decompose $y$ as:
$
  y = K x = sum_j μ_j (x, x_j) y_j
$

Which allows us to consider the inner product $(y, y_j)$:
$
  (y, y_j) = (K x, y_j) &= (sum_k μ_k (x, x_k) y_k, y_j) = sum_k [μ_k (x, x_k) (y_k, y_j)] = sum_k [μ_k (x, x_k) δ_(k j)]\
  &= μ_j (x, x_j)
$
Where the summation and inner product interchange by linearity of the inner product and $(y_k, y_j) = δ_(k j)$ by orthonormality of the basis functions.

Now we consider the norm of the difference between $tilde(R)_α K x$ and $x$:
$
  norm(tilde(R)_α K x - x)^2
  = norm(sum_j [μ_j^2 / (μ_j^3 + α) (y, y_j) x_j] - sum_j (x, x_j) x_j)^2
  = norm(sum_j [(μ_j^3 / (μ_j^3 + α) -1) (x, x_j) x_j ])^2\
$

Onto which we leverage that ${x_j}$ is complete by injectivity and orthonormal in order to apply Plancherel theorem / Parseval's Identity~@kirsch[Theorem A.15]:
$
  sum_(k=1)^∞ |(x, x_k)_X|^2 = norm(x)_X^2
$
In order to recover:
$
  norm(tilde(R)_α K x - x)^2
  &= sum_j [abs(μ_j^3 / (μ_j^3 + α) -1)^2 abs((x, x_j))^2]
$
Of which we take the limit $α → 0^+$:
$
  lim_(α→0^+) norm(tilde(R)_α K x - x)^2
  &= lim_(α→0^+) sum_j [abs(μ_j^3 / (μ_j^3 + α) -1)^2 abs((x, x_j))^2]
  &= sum_j [0 ⋅ abs((x, x_j))^2] = 0
$
Where the interchange of the sum and limit are justified by satisfying the conditions of _ Lebesque's Dominated Convergence Theorem_, which requires _pointwise convergence_ and _domination by an integrable function_:
$
"Pointwise Convergence:" &wide lim_(α→0^+) abs(μ_j^3 / (μ_j^3 + α) -1)^2 abs((x, x_j))^2 = 0 &&quad ∀ quad j\
"Domination:" &wide abs(μ_j^3 / (μ_j^3 + α) -1)^2 abs((x, x_j))^2 ≤ abs((x, x_j))^2 &&quad ∀ quad j, α > 0
$

As such, given that the norm of their difference goes to zero as $α → 0^+$, we conclude that $lim_(α→0^+) norm(tilde(R)_α K x - x) = 0$.

== Convergence in Deblurring Problem <sec:3.b>

We consider the blurring problem:
$
K x(t) = ∫_(-6)^6 k(t-s) x(s) dif s
wide "where" wide
k(z) = cases(
  1 + cos(π z \/ 3) quad & "if" quad abs(z) < 3\
  0 quad & "if" quad abs(z) ≥ 3
)
$

We select a test function $x(s)$ given by the sum of a Gaussian and a _top hat_ function:
$
  x(s) = exp(-2(s+2)^2) + H(s - 1) - H(s - 3)
$

And introduce measurement noise according to parameter $δ$ such that $y^δ = y+d$ where $norm(d) < δ$ and construct regularisation parameter $α=δ^3$ using the regularisation technique described in @sec:3.

In order to validate our implementation we construct a realisation for $δ=10^(-7)$ as seen in @fig:3_b_validation.

#figure(
  image("output/3_b_validation.png"),
  caption: [
    Validation of deblurring implementation with $δ=10^(-7)$ and $α=δ^3$.
  ]
) <fig:3_b_validation>

Satisfied that our implementation is correct, we proceed to investigate the convergence of the the regularisation regime by considering $norm(x_α - x)_2$ where $x_a = tilde(R)_(α(δ)) y^δ$ over a range of noise levels $δ ∈ [1, 10^(-7)]$ as seen in @fig:3_b_convergence.

#figure(
  image("output/3_b_convergence.png"),
  caption: [
    Convergence of regularised solution $x_α$ to true solution $x$ as noise level $δ → 0$ with $α=δ^3$.
  ]
) <fig:3_b_convergence>

We notice that the convergence rate appears to be of order $cal(O)(δ^(1\/5))$ and that the convergence in the two norm breaks at $δ≈5×10^(-6)$, which corresponds to $α=δ^3$ reaching the numerical precision of the computation, $ε≈2.22×10^(-16)$.

As such, we were able to numerically demonstrate convergence of the regularised solution to the true solution with decreasing noise levels.


// Appendix
#set heading(numbering: "A.1)")
#set math.equation(numbering: (c) => [(A.#(c))])
#counter(heading).update(0)

#pagebreak()
= Appendix
== Orthogonality of Sine Functions <app:sine_ortho>
We wish to show that the sine functions are orthogonal in $L^2(0, π)$:
$
∫_0^π sin(j s) sin(n s) dif s
&= ∫_0^π 1/2 [cos((n-j)s) - cos((n+j)s)] \
&= 1/2 [sin((n-j)s)/(n-j) - sin((n+j)s)/(n+j)]_0^π \
&= 1/2 [sin((n-j)π)/(n-j) - sin((n+j)π)/(n+j)]
$ <eq:sine_orthog>
Where $n, j ∈ ℕ_1$

We now consider the case $j ≠ n$ such that $n-j=p ∈ ℤ$ and $n+j=q ∈ ℤ$:
$
∫_0^π sin(j s) sin(n s) dif s
&= 1/2 [sin(p π)/p - sin(q π)/q] = 1/2[0 - 0] = 0
$

The other case, $j = n$ is more interesting and gives:
$
∫_0^π sin(j s) sin(n s) dif s
&= 1/2 [sin((n-j)π)/(n-j) - sin((2n)π)/(2n)]\
&= 1/2 sin((n-j)π)/(n-j)\
&= 1/2 ((n-j)π - (n-j)^3 π^3 + (n-j)^5 π^5 - …)\
&= π/2 + 0
$

Where we have used the power series definition of the sine function.

As such, we find that the integral in @eq:sine_orthog evaluates to:
$
∫_0^π sin(j s) sin(n s) dif s = π/2 δ_(j n)
$ <eq:sine_ortho_result>

Where $δ_(j n)$ is the Kronecker delta function.

== Optional Material for @sec:2 <app:2.optional>

We seek to prove compactness of $K$ by the integral criterion, which states that an integral operator is _compact_ if its domain is finite and the _kernel_ is square-integrable:
$
∫_0^1 ∫_0^1 |k(t, s)|^2 dif s dif t = ∫_0^1 ∫_0^t (1/sqrt(t-s))^2 dif s dif t = ∫_0^1 ∫_0^t 1/(t-s) dif s dif t < ∞
$

Where by inspection we observe a singularity at the line $s=t$, which we handle by introducing $K_ε$ with a small parameter $ε > 0$:
$
K_ε x(t) &= ∫_0^t 1/sqrt(t-s+ε) x(s) dif s
$

Such that the integral criterion becomes:
$
  ∫_0^1 ∫_0^1 |k_ε(t, s)|^2 dif s dif t = ∫_0^1 ∫_0^t 1/(t-s+ε) dif s dif t < ∞
$

Using substitution $u=t-s+ε$ we solve the inner integral:
$
∫_0^t 1/(t-s+ε) dif s
&= -∫_(t+ε)^ε 1/u dif u = [ln abs(u)]_ε^(t+ε) = ln(t+ε) - ln(ε)
$

Which is integrable and finite for $t ∈ [0, 1]$ and $ε > 0$ by inspection.

As such, we now need to show that $K_ε → K$ as $ε → 0$. We do this by considering the difference between the two operators:
$
D_ε = K - K_ε = ∫_0^t d_ε(t, s) x(s) dif s
$
Where $d_ε(t, s) = 1/sqrt(t-s) - 1/sqrt(t-s+ε)$ is the difference operator kernel.

We consider the following definition of the operator norm of the difference operator:
$
  norm(D_ε)_"op" = sup_(norm(x)_(L^2) = 1) norm(D_ε x)_(L^2)
$
We thus seek to investigate the 2-norm of difference operator subject to $norm(x)_(L^2) = 1$, to which end we employ _Young's Convolution Inequality_:
$
norm(f * g)_r ≤ norm(f)_p ⋅ norm(g)_q
$
Where $1\/p + 1\/q = 1\/r + 1$ and $1≤p,q,r≤∞$.

As such, with choices $r = p = 2, q=1$ we find:
$
  norm(D_ε)_(L^2) ≤ norm(d_ε (t, s))_(L^1) ⋅ norm(x)_(L^2) = norm(d_ε (t, s))_(L^1)
$

Which finally leads us to consider the $L^1$-norm of the difference operator kernel:
$
norm(d_ε (t, s))_(L^1) = ∫_0^t abs(1/sqrt(t-s) - 1/sqrt(t-s+ε)) dif s
$
Which we again solve using substitution $u = t-s$:
$
norm(d_ε (t, s))_(L^1) = ∫_0^1 (1/sqrt(u) - 1/sqrt(u+ε)) dif u = [2 sqrt(u) - 2 sqrt(u+ε)]_0^1 = 2 - 2sqrt(1+ε)) + 2 sqrt(ε)
$
Combining the results above, we are able to bound the operator norm:
$
norm(D_ε)_"op" ≤ sup (2 - 2sqrt(1+ε) + 2 sqrt(ε))
$

Taking the limit $ε → 0^+$ we find:
$
  norm(K - K_ε)_"op" → 0 quad "as" quad ε → 0^+
$

Since $K_ε → K$ in the operator norm as $ε → 0^+$ and $K_ε$ is compact, we conclude that $K$ is compact.

== Large Figures for @sec:1 <app:large_figures_1>

#figure(
  image("output/1_b_singular_values_large.png"),
  caption: [
    Singular values of the discretised operator $A$ for $n=128$ and $T=0.03$ compared against the analytical expression for the singular values of $K$.
  ]
)

#figure(
  image("output/1_b_singular_vectors_large.png"),
  caption: [
    3 select singular vectors and their associated singular functions for $n=128$ and $T=0.03$.
  ]
)

#figure(
  image("output/1_c_picard_plot_large.png"),
  caption: [
    Picard plot of the noisy data $y^δ$ with $2%$ relative noise.
  ]
)

#figure(
  image("output/1_c_reconstruction_large.png"),
  caption: [
    Reconstruction of $x$ using truncated SVD with $j^⋆ = 9$ compared against the true $x$.
  ]
)

== Large Figures for @sec:2 <app:large_figures_2>

#figure(
  image("output/2_a_large.png"),
  caption: [
     Noisy realisation of Abel problem with $x(s)=s^2$ perturbed by $5%$ relative white noise.
  ]
)

#figure(
  image("output/2_b_picard_large.png"),
  caption: [
    Picard plot of the noisy data $y^δ$ with $5%$ relative noise.
  ]
)

#figure(
  image("output/2_b_reconstructions_large.png"),
  caption: [
    Reconstructions of $x(t)$ using truncated SVD with varying truncation levels $k$ compared against the true $x(t)$.
  ]
)

#figure(
  image("output/2_c_tikhonov_reconstructions_large.png"),
  caption: [
    Reconstructions of $x(t)$ using Tikhonov regularisation with varying regularisation parameters $α$ compared against the true $x(t)$.
  ]
)

#figure(
  image("output/2_c_tsvd_vs_tikhonov_large.png"),
  caption: [
    Comparison of TSVD reconstruction with $k=20$ and Tikhonov reconstruction with $α=4×10^(-2)$.
  ]
)

#figure(
  image("output/2_d_discontinuous_large.png"),
  caption: [
    Comparison of TSVD reconstructions with $k∈{20, 50}$ and Tikhonov reconstruction with $α=10^(-2)$ for a discontinuous function.
  ]
)

== Large Figures for @sec:3 <app:large_figures_3>

#figure(
  image("output/3_b_validation_large.png"),
  caption: [
    Validation of deblurring implementation with $δ=10^(-7)$ and $α=δ^3$.
  ]
)

#figure(
  image("output/3_b_convergence_large.png"),
  caption: [
    Convergence of regularised solution $x_α$ to true solution $x$ as noise level $δ → 0$ with $α=δ^3$.
  ]
)

== Code for @sec:1.b <app:code_1_b>
#raw(read("output/1_b.py"), block: true, lang: "python")

== Code for @sec:1.c <app:code_1_c>
#raw(read("output/1_c.py"), block: true, lang: "python")

== Code for @sec:2.a <app:code_2_a>
#raw(read("output/2_a.py"), block: true, lang: "python")

== Code for @sec:2.b <app:code_2_b>
#raw(read("output/2_b.py"), block: true, lang: "python")

== Code for @sec:2.c <app:code_2_c>
#raw(read("output/2_c.py"), block: true, lang: "python")

== Code for @sec:2.d <app:code_2_d>
#raw(read("output/2_d.py"), block: true, lang: "python")

== Code for @sec:3.b <app:code_3_b>
#raw(read("output/3_b.py"), block: true, lang: "python")

#bibliography("refs.bib")

