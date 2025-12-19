Excellent solutions, very impressive! Great theory, computations and plots. Approved.

Problem 1

a) Good. We need that x_j and y_j spans L2. That's known from Fourier series, I would say. Observe also that mu_j>0 (and tend to zero as expected).

b) Nice discretization. When evaluating Ax given x, I strongly recommend an FFT based approach (nlog(n) instead of n^2 complexity). Important in higher dimensions.I am surprised to see that the amplitude of x_j and u_j are the same, since they 

c) Excellent. Is the naive solution really not doing worse - no blow up?

Problem 2

a) We can analytically compute y(t)

b) The problem is very mildly ill-posed. Very good analysis. Thanks for the justification of x_j(1) = 0; I have learned something here!

c) The normal equation is great for analysis, but when computing the Tikhonov solution, always use the augmented least-squares problem.

d) Fine

Problem 3

a) Excellent

b) We actually need to make sure that delta^3/alpha(delta) goes to zero as delta goes to zero. My instruction was not clear at this point. Try for instance alpha(delta) = delta^2 in Figure 13.
