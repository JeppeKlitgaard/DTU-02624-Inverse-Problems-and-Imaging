#import "../presentation_template.typ": *
#import "@preview/physica:0.9.7": TT, pdv, dv
#show: doc => setup(3, doc)

= Introduction - Homework 3

- 2D Backwards Heat Equation
  - Application on Image
  - Self-adjointness
  - Landweber Iterations
- Total Variation Inpainting
  - TV Inpainting
  - Harmonic Inpainting
  - Noise vs Noiseless
- Feedback received on report
  - Mis-stated heat equation as convolution
    - Not strictly a convolution since kernel is not shift-invariant
    - Speed-up is similar to convolution, we turn sum-of-sums into pointwise multiplication
= 2D Backwards Heat Equation
#progress-marker-lower
We consider the 2D backwards heat equation on domain $Ω = (0, π) × (0, π)$ with homogeneous Dirichlet boundary conditions as given by forward operator:
$
  K x (t_1, t_2) = ∫_0^π ∫_0^π x(s_1, s_2) k(t_1, t_2, s_1, s_2) dif s_1 dif s_2 \
  k(t_1, t_2, s_1, s_2) = cancel(4/π^2) ∑_(n=1)^∞ ∑_(m=1)^∞ e^(-(n^2 + m^2) T) sin(n t_1) sin(m t_2) sin(n s_1) sin(m s_2)
$

- $T$ is terminal time
- Consider Gaussian decay in kernel:
  - Higher $T$ ⇒ more ill-posed
  - Higher modes $n, m$ ⇒ more ill-posed
- Problem is satisfies uniqueness and existence in Hadamard sense, but *not* stability

== 2D Backwards Heat Equation — Implementation
- Homogenous Dirichlet BCs work well with sine basis
- We model only the interior of domain since boundary is fixed
- Can be done as matrix-free operator with 2D FFT ($n log n$ complexity)
  - This requires odd-extension to produce sine basis: "Unfold" input by mirroring in both dimensions
  - Apply 2D FFT, multiply in frequency domain, apply inverse 2D FFT
  - Alternative is to use Discrete Sine Transform (DST), could prevent odd-extension and lower memory allocations
  - Math simple: $y = ℱ^(-1) (e^(- (n^2 + m^2) T) ⋅ ℱ(x))$
  - Implementation a bit tricker due to wave-numbers, orderings, odd-extension

== 2D Backwards Heat Equation — Validation of Implementation

#side-by-side(columns: (70%, 1fr))[
#align(center + horizon, image("output/1_b.png"))
][
- Test using single mode
- $x_"test" (s_1, s_2) = sin(k s_1) sin(l s_2)$
- $k=4, l=3$
- Forward operator just attenuates by $e^(- (k^2 + l^2) T)$
- Observe excellent agreement
- Error $e ∼ 10^(-20)$
]

== 2D Backwards Heat Equation — Apply to Image

#side-by-side(columns: (70%, 1fr))[
#align(center + horizon, image("output/1_c.png"))
][
- Apply to an image
- Forward model just blurs image according to diffusion
- Terminal time $T=0.001$
- Noise level $σ=10^(-2)$
]

== 2D Backwards Heat Equation — Self-Adjointness

Self-adjointness of operator $K$ is easily proven from symmetry of kernel
$
  K x (t_1, t_2)
  &= ∫_0^π ∫_0^π k(t_1, t_2, s_1, s_2) x(s_1, s_2) dif s_1 dif s_2\
  &= ∫_0^π ∫_0^π k^*(dash(s)_1, dash(s)_2, t_1, t_2) x(s_1, s_2) dif s_1 dif s_2
  = K^* x (t_1, t_2),
$ <eq:1_self_adjointness_long>
where $k^*$ is adjoint of the kernel
$
  k(t_1, t_2, s_1, s_2)
  &= sum_(n, m) β_(n m) sin(n s_1) sin(m s_2) sin(n t_1) sin(m t_2)\
  &= k(s_1, s_2, t_1, t_2)
  = k^*(dash(s)_1, dash(s)_2, t_1, t_2)
$
with $dash(s)_i = s_i$ since $s_i ∈ ℝ$ and $β_(n m) = e^(-(n^2 + m^2) T)$

Thus, $K = K^*$

== 2D Backwards Heat Equation — Landweber Iterations

Landweber iterations are also known as _gradient descent_ and are given by the iteration
$
  x_(j+1) = x_j + ω K^* (K x_j - y^δ)
$
- Discretisation and self-adjointness of $K$ yields $x_(j+1) = x_j + ω A (A x_j - y^δ)$
- Relaxation parameter $ω$ must satisfy $0 < ω < 2 \/ σ_"max"^2$
  - We may find this analytically as $σ_"max" = e^(-2 T)$
  - We let $ω = 18/10 e^(4 T)$ to have safety margin of $9/10$
- Method is _semi-convergent_
  - Finds local minimum
  - Early stopping regularises solution and prevents overfitting to noise

#v(1em)
With this, we are able to reconstruct image from blurred and noisy observation.

== 2D Backwards Heat Equation — Landweber Reconstruction
#progress-marker-lower

#side-by-side(columns: (1fr, 65%))[
- Note wave-like artefacts near boundary due to image not satisfying homogenous Dirichlet BCs
- These are reminiscent of the familiar Gibbs phenomenon
- Error concentrated on high-frequency edges as expected
- Does not converge to exact solution
- Tried CGLS, quickly diverges (≈3 iterations)
][
#image("output/1_d_images.png")
#image("output/1_d_convergence.png")
]

= Inpainting
#progress-marker-upper

#side-by-side(columns: (70%, 1fr))[
We recall the two inpainting models considered:
$
  min_u& 1/2 norm(𝟙_(Ω \\ cal(M)) (z - u))^2 + λ ∫_Ω abs(∇ u (x))_2 dif x &&wide "TV Inpainting"\
  min_(u ∈H_0^1(Ω))& 1/2 norm(𝟙_(Ω \\ cal(M)) (z - u))^2 + λ/2 norm(∇ u (x))^2 dif x &&wide "Harmonic Inpainting"
$
Where $λ$ is a regularisation parameter and $z$ is the the remaining parts of the corrupted image as given by
$
  z = cases(
    dash(u) + η wide& "in" Ω \\ cal(M),
    "unknown" wide& "in" cal(M)
  )
$
Note that in the noiseless case $η=0$ the data-fidelity term vanishes if imposing known pixels as hard constraints, which allows the regularisation parameter $λ$ to be absorbed.
][
#figure(
  image("KarinHarboStripe.png", width: 100%),
  caption: [
    Image to be inpainted
  ]
)
]

== TV Inpainting
#progress-marker-lower
The TV inpainting is carried out using `cvxpy` in both the regular case and the case of assumed noiselessness.

#align(center, image("output/2_b_tv_comparison.png", width: 80%))

- TV tries to minimise total variation, i.e., the length of edges
- TV liable to staircasing. $∂ℳ$ is one _"stair"_
- Good candidate for TV due to sharp, straight edge
- Noiseless case seems preferable here

== Harmonic Inpainting
#progress-marker-lower
Harmonic inpainting is carried out in a similar fashion using `cvxpy`, again for both cases.

#align(center, image("output/2_b_harmonic_comparison.png", width: 80%))

- Harmonic inpainting minimises gradients ⇒ smooth solutions
- This is better at smoothing transition along $∂ℳ$
- Significantly worse recovery of sharp, straight edge
- Noiseless case again seems preferable
