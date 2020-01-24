function [A, sigma, K] =myLevinson(correlation)
n=20;
K = zeros(1, n+1); % remember to cut away the first element
A = [1];

sigma = zeros(1, n);
sigma(1) = correlation(1);

for i=1:n
    delta = correlation(i+1:-1:2)*A';
    K(i+1) = -delta/sigma(i);
    A = [A, 0] + K(i+1)*[0, A(i:-1:1)];
    sigma(i+1) = sigma(i)*(1-K(i+1)^2);
end
K = K(2:n+1)';
