% Example of Picard plot for Numerical differentiation problem
% Kim Knudsen, DTU Compute, 2023

%Call NumericalDifferentiation to get matrix A, yvec and xvec
NumericalDifferentiation;
y = yvec;

% Additive noise.
delta = 0.1;
d=randn(size(y));
d=delta*norm(y)/norm(d)*d;

% Add noise.
yd = y+d;

% Compute the SVD.
[U,S,V] = svd(A);
s = diag(S);

% Expansion coefficients of yd in U basis.
ycoeff = U'*yd;

%Picard plot
figure(3) , clf

semilogy(s)
hold on
semilogy(abs(ycoeff),'*','markersize',4);
semilogy(abs(ycoeff./s),'+','markersize',4)

%xlim([1,400])

legend('\mu_j','|u_i^T y_d|','|u_i^T y_d/\mu_j|')
title('Picard plot')
