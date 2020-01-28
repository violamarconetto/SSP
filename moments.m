function [r]=moments(A, f, noise_variance)
r = zeros(1, 3);
r(1) = (A^2)/2 + noise_variance;
r(2) = (A^2)/2 * cos(2*pi*f);
r(3) = (A^2)/2 * cos(4*pi*f);