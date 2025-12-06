#import "@preview/physica:0.9.6": dv
#import "@preview/physica:0.9.6": iprod as ip
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

#set document(title: [02624 Inverse Problems and Imaging — Homework 3])
#align(center, title())

Jeppe Klitgaard `<s250250@dtu.dk>`\
Date: 2025-12-05

#outline()

= Introduction <sec:intro>

This report presents solutions to Homework 3 for the course 02624 Inverse Problems and Imaging.

GitHub Copilot was used for tab completions in the development of the code associated with this report.

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

= 2D Heat Problem <sec:1>

We consider the 2D backwards heat problem on the domain $(0, π)^2$ with homogenous Dirichlet boundary conditions as given by the forward operator

$
  K x(t_1, t_2) = ∫_0^π ∫_0^π k(t_1, t_2, s_1, s_2) x(s_1, s_2) dif s_1 dif s_2,
$ <eq:1_forward_operator>
with kernel
$
  k(t_1, t_2, s_1, s_2) = sum_(n,m) e^(-(n^2 + m^2) T) sin(n s_1) sin(m s_2) sin(n t_1) sin(m t_2)
$ <eq:1_kernel>
where $x(s_1, s_2)$ is the image constituting the initial condition. The output of the model, $y(t_1, t_2) = K x(t_1, t_2)$ represents the observed heat distribution at a time $T$.

== Ill-posedness of the 2D Backwards Heat Problem <sec:1a>

While the backwards heat problem may in some circumstances abide by both the existence and uniqueness criterion of the Hadamard definition of well-posedness, it does not lend itself to stable inversion, which can be appreciated by considering the exponential decay of the singular values with respect to time in the forward operator $K$. Upon inspection of @eq:1_kernel, it is clear that the higher frequency components are rapidly attenuated as time $T$ increases, which leads to instability with respect to even minute perturbations.

== Implementation of the Forward Operator <sec:1b>

Upon noticing that @eq:1_forward_operator constitutes a _convolution_, we realise that significant computational savings are available by leveraging the favourable scaling of the Fast Fourier Transform (FFT), by which we may transform to and from the frequency domain with $n log(n)$ complexity. This allows us to carry out the convolution as pointwise multiplications in the frequency domain, which circumvents $n^2$ complexity entirely.

Given that we are working on a domain $(0, π)^2$ with homogenous Dirichlet boundary conditions, we choose to model only the interior of the domain, which also enables us to use the sine polynomial basis.

The implementation of the forward operator is given in @app:1b_code_forward_operator.

We are able to validate our implementation by constructing a test case in which we only activate a single frequency component:
$
  x_"test" (s_1, s_2) = sin(k s_1) sin(l s_2)
$
where we let $k=4, l=3$.

Inspection of the forward operator and @eq:1_kernel confirms that the inactive modes will remain zero, while the active mode will be attenuated by a factor
$
  α = e^(-(k^2 + l^2) T).
$

This allows us to compute an analytical solution $y_"exact"$ and $y_"numerical"$, from which we can compute the error, $e$:
$
  y_"exact" (t_1, t_2) &= α sin(k t_1) sin(l t_2),\
  y_"numerical" (t_1, t_2) &= K x_"test" (t_1, t_2),\
  e &= y_"exact" - y_"numerical".
$

The test case is visualised in @fig:1b, which reveals that the error is of order $e~10^(-20)$, which is negligible and may be attributed to numerical round-off errors.

As such, we conclude that our implementation of the forward operator is correct.

#figure(
  image("output/1_b.png"),
  caption: [
    Test case for FFT-based implementation of the forward operator for the 2D heat equation.
  ]
) <fig:1b>

== Application on Image <sec:1c>

We apply the forward operator on an image, `KarinHarbo512.png`, which is converted to a grayscale representation and normalised such that the pixel values lie in the interval $[0, 1]$.

This image is then passed through the forward operator with terminal time $T=0.001$ and a Gaussian noise of standard deviation $δ=0.01$ is added to produce $y^δ$, which is depicted in @fig:1c.

#figure(
  image("output/1_c.png"),
  caption: [
    Application of the forward operator on image `KarinHarbo512.png` with added Gaussian noise.
  ]
) <fig:1c>

== Self-adjointness of the Forward Operator and Inversion using Landweber Iterations <sec:1d>

=== Self-adjointness of the Forward Operator <sec:1d_adjointness>
Using a similar treatment as in Homework 1, we may prove self-adjointness of the forward operator $K$ defined in @eq:1_forward_operator by considering its action on $x ∈ L^2((0, π)^2)$:
$
  K x (t_1, t_2)
  &= ∫_0^π ∫_0^π k(t_1, t_2, s_1, s_2) x(s_1, s_2) dif s_1 dif s_2\
  &= ∫_0^π ∫_0^π k^*(dash(s)_1, dash(s)_2, t_1, t_2) x(s_1, s_2) dif s_1 dif s_2
  = K^* x (t_1, t_2),
$ <eq:1_self_adjointness_long>

Where $k^*$ is the adjoint of the kernel $k$ and $dash(s)_i$ denotes the complex-conjugate of $s_i$. The second equality follows from symmetry of the kernel $k$:
$
  k(t_1, t_2, s_1, s_2)
  &= sum_(n, m) β_(n m) sin(n s_1) sin(m s_2) sin(n t_1) sin(m t_2)\
  &= k(s_1, s_2, t_1, t_2)
  = k^*(dash(s)_1, dash(s)_2, t_1, t_2)
$
where $β_(n m) = e^(-(n^2 + m^2) T)$ and $dash(s)_i = s_i$ since $s_i ∈ ℝ$.

As such, from @eq:1_self_adjointness_long, we conclude that $K = K^*$ — that is, we conclude that the forward operator is self-adjoint.

An alternative approach would be to prove self-adjointness by equating the inner products $ip(K x, y)$ and $ip(x, K y)$ for arbitrary functions $x, y ∈ L^2((0, π)^2)$.

=== Inversion using Landweber Iterations <sec:1d_landweber>
This allows us to easily compute Landweber iterations, also referred to as _gradient descent_ within optimisation theory, which we recall to be a semi-convergent iterative algorithm for solving ill-posed linear inverse problems. It may be expressed in the form
$
  x_(j+1) &= x_j + ω K^* (K x_j - y^δ)
$

Which given self-adjointness of $K$ and approximation of $K$ by a matrix $A$ may be expressed as
$
  x_(j+1) &= x_j + ω A (A x_j - y^δ)
$
where $j$ denotes the iteration index and $ω$ is a relaxation parameter satisfying $0 < ω < 2/σ_"max"^2$ with $σ_"max"$ being the largest singular value of the operator $K≈A$.

Given the sine polynomials as our basis, we find that the singular values are the eigenvalues of the operator, which may be read off directly from @eq:1_kernel as
$
  σ_(n m) = e^(-(n^2 + m^2) T).
$
where $n, m ∈ ℕ_(1)$ implies that the largest singular value is given by $n = m = 1$, yielding
$
  σ_"max" = e^(-2 T).
$

We choose a relaxation parameter $ω = 9/10 ⋅ 2/σ_"max"^2$ where $9/10$ acts as a safety margin.

This lets us carry out $N_"iter" = 1500$ Landweber iterations, which produces @fig:1d_images, a somewhat grainy reconstruction of the original image shown the left panel of @fig:1c, though notably sharper than the blurred version shown in the right panel of the same figure.

We note that the reconstruction error is concentrated around the sharp edges of the original image, which is consistent with the analysis of @sec:1a, wherein we noted that higher frequency components are more severely attenuated by the forward operator, leading to more increased instability in their reconstruction.

Additionally, we observe a wave-like artefact around the borders of the image which we attribute to the fact that the original image does not tend to zero at the boundaries as required by the homogenous Dirichlet boundary conditions of the model. This continuity mismatch leads to the observed oscillations, which are reminiscent of the Gibbs phenomenon.

#figure(
  image("output/1_d_images.png"),
  caption: [
    Landweber reconstruction of the initial condition from noisy observations at time $T=0.001$ after $N_"iter" = 1500$ iterations.
  ]
) <fig:1d_images>

Given the semi-convergence of Landweber iterations we choose to inspect the error history over the
iterations calculated in the $L^2$ and $L^∞$-norms, which is depicted in @fig:1d_convergence.

We can appreciate the form of @fig:1d_convergence by again considering the Landweber iterations as a gradient descent scheme, which makes it clear that it is semi-convergent in that it does not necessarily converge to the true solution, but rather to a local minimum of the residual norm, or in this case error norm. This is contrasted by the behaviour of other iterative methods such as the Conjugate Gradient (CG) method, in which the semi-convergence leads to 'overshooting' the minimum and divergence of the error norms after a certain number of iterations is observed. Lastly, we hypothesise that this could also happen in the case of Landweber iterations if instability is promoted through numerical errors or significant high-frequency noise.

This could be further investigated by increasing the number of allotted of iterations, $N_"iter"$, though at the cost of increased computation time and was thus not pursued further in this report.

Further iterations and alternative implementations using Conjugate Gradient (CG) may be investigated as part of the preparatory work for the examination.

#figure(
  image("output/1_d_convergence.png"),
  caption: [
    Error history of Landweber iterations in the $L^2$ and $L^∞$-norms.
  ]
) <fig:1d_convergence>

= Total Variation Inpainting <sec:2>

We recall the following variational formulation for total variation (TV) inpainting:
$
  min_u 1/2 underbrace(norm(𝟙_(Ω \\ cal(M)) (z - u))^2, "data fidelity") + λ underbrace(∫_Ω abs(∇ u (x))_2 dif x, "regularisation") wider wider "TV Inpainting"
$ <eq:2_tv_inpainting>
where $λ$ is a regularisation parameter and $z$ is the the remaining parts of the corrupted image as given by
$
  z = cases(
    dash(u) + η wider& "in" Ω \\ cal(M),
    "unknown" wider& "in" cal(M)
  )
$
with $dash(u)$ being the original uncorrupted image, $η$ being noise and $cal(M)$ being the mask of corrupted pixels inside the full domain $Ω$.

We note that in the noiseless case, $η=0$, the data fidelity term vanishes and the regularisation parameter, $λ$, may be set to unity without loss of generality.

We may also consider _harmonic_ inpainting, which may be expressed on variational form as:
$
  min_(u ∈H_0^1(Ω)) 1/2 norm(𝟙_(Ω \\ cal(M)) (z - u))^2 + λ/2 norm(∇ u (x))^2 dif x wider wider "Harmonic Inpainting"
$ <eq:2_harmonic_inpainting>
Where $H_0^1(Ω)$ denotes the Sobolev space of functions with square-integrable first derivatives that vanish on the boundary of $Ω$. and can be viewed as a generalised Tikhonov regularisation with $H^1$-norm penalty.

== Displaying `KarinHarboStripe.png` <sec:2a>

We display the partially corrupted image `KarinHarboStripe.png` by the same procedure used in @sec:1c, wherein the image is converted to grayscale and normalised to the interval $[0, 1]$ in order to produce @fig:2a.

#figure(
  image("output/2_a.png", height: 70mm),
  caption: [
    Display of the corrupted image `KarinHarboStripe.png`. The black region is the corrupted part, $cal(M)$, to be inpainted.
  ]
) <fig:2a>

== TV and Harmonic Inpainting <sec:2b>

We may solve the TV inpainting problem defined in @eq:2_tv_inpainting by formulating the problem in `cvxpy`, and qualitatively selecting $λ=0.01$ as an optimal regularisation parameter which preserves most of the features of the uncorrupted image while inpainting the masked region.

With this, we are able to produce the reconstruction shown in the centre panel of @fig:2b_tv_comparison.

Alternatively, under the no-noise assumption, $η=0$, we may set $λ=1$ and solve the problem with the hard constraint $u = z ∉ cal(M)$ to produce the right panel of @fig:2b_tv_comparison.

#figure(
  image("output/2_b_tv_comparison.png", height: 50mm),
  caption: [
    Total Variation inpainting given regularisation parameter $λ=0.01$. The result under the noiseless assumption is shown in the right panel.
  ]
) <fig:2b_tv_comparison>

Similarly, we may solve the harmonic inpainting problem defined in @eq:2_harmonic_inpainting by again formulating the problem in `cvxpy`, and produce the reconstruction shown in the centre panel of @fig:2b_harmonic_comparison with $λ=0.1$, while the noiseless case is shown in the right panel of the same figure.

#figure(
  image("output/2_b_harmonic_comparison.png", height: 50mm),
  caption: [
    Harmonic inpainting given regularisation parameter $λ=0.1$.
  ]
) <fig:2b_harmonic_comparison>

#pagebreak()
== Inpainting Conclusions <sec:2c>

We may collate the results from @sec:2b into a single figure for easier comparison, which is shown in @fig:2c_comparison.

It is evident that the noiseless cases produce better results while having the additional benefit that correctly selecting the regularisation parameter is no longer required. However, this is only under the assumption that the image is truly noiseless to begin with, which is not likely to be the case.

In this case, the corrupted region contains a single, relatively straight line of high constrast, which lends itself to the TV inpainting method, which is able to reconstruct the line with high fidelity. However, it has the notable disadvantage that within both the light and dark regions of the image, an unnatural boxy artefact is introduced, which is smoothed out in the harmonic inpainting method. This can be understood by appreciating that the TV method seeks to reduce the total variation of the image, which promotes piecewise-constant solutions, while the harmonic method seeks to reduce the $L^2$-norm of the gradient, which promotes smooth solutions.

The harmonic method is not able to reconstruct the high-contrast line nearly as well as the TV method and instead lets the dark and light regions bleed into each other.

#figure(
  image("output/2_c.png", height: 150mm),
  caption: [
    Comparison of inpainting methods.
  ]
) <fig:2c_comparison>


// Appendix
#set heading(numbering: "A.1)", supplement: [appendix])
#set math.equation(numbering: (c) => [(A.#(c))])
#counter(heading).update(0)

#pagebreak()

= Appendix

== Code preamble <app:code_preamble>
#raw(read("output/preamble.py"), block: true, lang: "python")

== Code for @sec:1b forward operator <app:1b_code_forward_operator>
#raw(read("output/1_b_forward_operator.py"), block: true, lang: "python")

== Code for @sec:1b validation <app:1b_code_validation>
#raw(read("output/1_b_test.py"), block: true, lang: "python")

== Code for @sec:1c validation <app:1c>
#raw(read("output/1_c.py"), block: true, lang: "python")

== Code for @sec:1d validation <app:1d>
#raw(read("output/1_d.py"), block: true, lang: "python")

== Code for @sec:2a validation <app:2_a>
#raw(read("output/2_a.py"), block: true, lang: "python")

== Code for @sec:2b validation <app:2_b>
#raw(read("output/2_b.py"), block: true, lang: "python")

== Code for @sec:2c validation <app:2_c>
#raw(read("output/2_c.py"), block: true, lang: "python")

// #pagebreak()
// #bibliography("refs.bib")
