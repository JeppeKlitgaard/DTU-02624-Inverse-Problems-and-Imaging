#import "@preview/physica:0.9.6": dv
#import "@preview/physica:0.9.4": curl, grad, tensor, pdv, dv, eval, TT
#import "preamble.typ": mref, mathformatter
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.8": *

//! Styling
#set page(
  numbering: "1/1",
  margin: (
    top: 1.5cm,
    bottom: 1.5cm,
    left: 1.5cm,
    right: 2.5cm,

  )
)

#set heading(numbering: "1.a.1)")
#set math.equation(numbering: "(1)")

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

#align(center, text(17pt)[
  *02624 Inverse Problems and Imaging — Homework 1*
])

Jeppe Klitgaard `<s250250@dtu.dk>`\
Date: 2025-10-09

//! Content
= Backwards Heat Equation <sec:1>
We are given the _forward model_ of the _backwards heat equation_:
$
K x(t) &= ∫_0^π x(s) k(t, s) dif s\
k(t, s) &= 2/π ∑_(n=1)^∞ e^(-T n^2) sin(n t) sin(n s)
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
K x_j &= μ_j y_j \
K^* y_j &= μ_j x_j
$ <eq:singular_system_def>

We first verify that $K x_j = μ_j y_j$:
$
K x_j (t)
&= ∫_0^π x_j (s) k(t, s) dif s\
&= ∫_0^π sqrt(2/π) sin(j s) (2/π ∑_(n=1)^∞ e^(-T n^2) sin(n t) sin(n s)) dif s\
&= (2/π)^(3/2) (∑_(n=1)^∞ e^(-T n^2) sin(n t) ∫_0^π sin(j s) sin(n s) dif s)\
&= (2/π)^(3/2) (∑_(n=1)^∞ e^(-T n^2) sin(n t) ∫_0^π sin(j s) sin(n s) dif s)\
&= (2/π)^(3/2) (∑_(n=1)^∞ e^(-T n^2) sin(n t) π/2 δ_(j n))\
&= e^(-T j^2) sqrt(2/π) sin(j t)\
&= μ_j y_j (t)
$

Where the summation and integral interchange was enabled by uniform convergence of the series through inspectation of @eq:1, in which the Gaussian decay factor $e^(-T n^2)$ ensures rapid convergence. The orthogonality of the sine function in $L^2(0, π)$ is shown in Appendix @app:sine_ortho, where we find that the
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

From which we may conclude that $K$ is _self-adjoint_:
$
  K = K^*
$

As such, we may repeat treatment above for the equation $K^* y_j = μ_j x_j$, which becomes true by inspection.

In addition to this, the singular functions must be orthonormal in $L^2(0, π)$, which they are by the treatment in @app:sine_ortho and the singular values must be non-negative, which is indeed the case for the given system.

Lastly, we note that the singular values decay to zero as $j → ∞$ since $T > 0$, as we would expect for a compact operator.

We conclude that $K$ is indeed a compact operator with singular system $(x_j, y_j, μ_j)$ as given in @eq:singular_functions and @eq:singular_values.

== Discretisation of $K$ <sec:1.b>

We seek to discretise $K$ by a matrix $A ∈ ℝ^(n × n)$ where $n=128$ and $T=0.03$, to which end the spatial and temporal domains $[0, π]$ are discretised into $n$ equidistant points:
$
  t_i &= (k π)/(n+1) wider k=1, …, n\
  s_i &= (k π)/(n+1) wider k=1, …, n\
$
Where we note that the ends of the interval are not included as the sine functions are zero at these points, which would not result in an orthogonal basis.

Next, we construct a matrices $U, V$ whose columns are the left and right singular vectors, $x_j, y_j$, sampled at the points $s_i$:
$
  U_(i j) &= sqrt(2/(n+1)) sin(j ⋅ t_i)\
  V_(i j) &= sqrt(2/(n+1)) sin(j ⋅ s_i)\
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
determined by the numerical precision of the computation, which is approximately $10^(-16)$ in this case, as demarked by the dotted black line.

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

Finding a suitable choice of the threshold index $j^⋆$ is undertaking through the construction of a _Picard plot_ containing the singular values $μ_j$, the expansion coefficients $abs((y^δ, y_j))$ and the solution coefficients $abs((y^δ, y_j))\/μ_j$ as seen in @fig:1_c_picard_plot.

#figure(
  image("output/1_c_picard_plot.png"),
  caption: [
    Picard plot of the noisy data $y^δ$ with $2%$ relative noise.
  ]
) <fig:1_c_picard_plot>

We recall the discrete Picard condition~@hansen[3.16], which states that the expansion coefficients must decay faster than the singular values for a stable solution to exist. We use this in conjunction with @fig:1_c_picard_plot to choose a threshold index $j^⋆ = 9$ at which the solution coefficients begin to diverge, indicating that noise is dominating the solution. We note that this choice is somewhat subjective and given we only have a single realisation of the noise available in the plot, we conservatively place it at $j^⋆ = 9$ as opposed to $j^⋆ = 10$.

With this, we are able to compute a reconstruction of $x$ using the truncated SVD:
$
  U Σ V^† = A
$
From which we may construct the truncated matrices by reducing the number of columns of $U, V$ and the size of $Σ$ to $j^⋆$.

From this, we may reconstruct $x$ using @eq:singular_system_def:
$
  x_"rec" = ∑_(j=1)^(j^⋆) ((y^δ, y_j))/μ_j x_j = underbrace(tilde(V) tilde(Σ)^(-1) tilde(U)^†, tilde(A)^(-1)) y^δ
$

Where the tildes denote that the matrices are truncated to the first $j^⋆$ columns/rows.

The resulting reconstruction is shown in @fig:1_c_reconstruction, where we observe a good match between the true model $x$ and the reconstruction $x_"rec"$. The noisy measurement $y^δ$ is also shown alongside the true signal $y$.

We observe that the regularised inversion is able to recover the main features of $x$ despite the presence of significant noise in the measurements. Attempting to invert the problem without regularisation results in a completely unusable solution dominated by noise, as seen by $x_"naive"$ in @fig:1_c_reconstruction, whose values are to be read off the right vertical axis.

#figure(
  image("output/1_c_reconstruction.png"),
  caption: [
    Reconstruction of $x$ using truncated SVD with $j^⋆ = 9$ compared against the true $x$.
  ]
) <fig:1_c_reconstruction>

= Abel Problem

We are given the _forward model_ of the _Abel problem_:
$
K x(t) &= ∫_0^t 1/sqrt(t-s) x(s) dif s
$ <eq:2_abel_forward>

Defined on the domain $L^2(0, 1)$.

#counter(heading).update((2, ))
==

// Appendix
#set heading(numbering: "A.1")
#set math.equation(numbering: (c) => [(A.#(c))])
#counter(heading).update(0)

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

== Code for @sec:1.b <app:code_1_b>

#raw(read("output/1_b.py"), block: true, lang: "python")

== Code for @sec:1.c <app:code_1_c>

#raw(read("output/1_c.py"), block: true, lang: "python")

#bibliography("refs.bib")
