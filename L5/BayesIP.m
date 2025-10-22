%% Illustration of Randomize then Optimize method for approximation of solution to Beyaes Inverse Problem
%% Kim Knudsen, 02624, Sept 2025.

clear 
close all

%Abel integral operator
% Course 02624
% Problem: Approximation of the Abel integral operator Kx(t) =  int_0^t x(s)/sqrt(t-s)ds as a matrix A using numerical quadrature on 
% off set s and t grids to avoid the singularity.
% Kim Knudsen 2025

  
%clear;
%close all;
N= 256;
%N = 512;
h= 1/N;

svec = h/2:h:(1-h/2);
tvec = (svec(1:end)+h/2)';

sMatrix = repmat(svec,N,1);
tMatrix = repmat(tvec,1,N);

%Evaluate matrix A
A = zeros(N,N);
index = sMatrix<tMatrix;
A(index) = h * 1./sqrt(abs(tMatrix(index)-sMatrix(index)));


%NumericalDifferentiation;
AbelIntegralForward;
%Provides matrix A, xvec and yvec.

%Prior variance
gamma = 0.02;

%Noise variance
sigma = 0.001;

%Regularization parameter
alpha= sigma^2/gamma^2;

%Sweep through different noise realizations and solve Tikhonov problem
figure(1)
hold on

%No samples
K=100;
xalpha = zeros(N,K);

for k = 1:K
    noise = gamma*randn(N,1);
    yvecnoise = yvec + noise;
    norm(noise)/norm(yvec)
    yext = [yvecnoise;zeros(N,1)];
    Aalpha  = [A;sqrt(alpha)*eye(N)];
    xalpha(:,k) = Aalpha \ yext; 
    plot(svec,xalpha(:,k))
    %pause
end;

plot(svec,xvec,'*')

pause
xCM = sum(xalpha,2)/K;

figure(2);
plot(svec,xCM)
hold on
plot(svec,xvec)

