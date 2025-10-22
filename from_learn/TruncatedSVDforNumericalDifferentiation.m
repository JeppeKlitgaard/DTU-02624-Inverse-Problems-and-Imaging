% Example of Truncated SVD reconstruction for Numerical differentiation problem
% For 02624
% Kim Knudsen, DTU Compute, 2025

% Call NumericalDifferentiation to get matrix A, yvec, yvecnoise and xvec
NumericalDifferentiation;

close all;


% SVD analysis
[U,S,V] = svd(A);

% Inspect singular vectors/functions
figure(1)
plot(s,V(:,1:4))
title('First singular functions for numerical differentiation')

%True singular values
figure(2)
semilogy(diag(S));
strue=1./(((1:N)-1/2)*pi);
hold on;
semilogy(strue);

legend('Computed \mu_j','True \mu_j')
xlabel('j')
title('Singular values')

% Truncated SVD
% k determines number of terms in SVD sum; should be chosen in accordance
% with noise level
k = 10;
xk = V(:,1:k)*S(1:k,1:k)^(-1)*U(:,1:k)'*yvecnoise;

%Inspect result
figure(3)
plot(s,xvec)
hold on;
plot(s,xk);
legend('x(s)','xk(s)')
title('TSVD solution')
