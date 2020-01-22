function [sigma, A, K] = levinson(correlation)
n=20;
K = zeros(1, n-1);
sigma = zeros(1, n);
A = 1;
sigma(1) = correlation(1);

for i=2:n
    delta = correlation(i:-1:2)*A;
    K(i-1) = -delta/sigma(i-1);
    A_new = [A, 0] + K(i-1)*[0, A((i-1):-1:1)];
    sigma(i) = sigma(i-1)*(1-K(i-1)^2);
    A = A_new;
end