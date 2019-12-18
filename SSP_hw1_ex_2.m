clear all
close all 
clc
%% ex.2 point f
y_lim = [15:50];
n = size(y_lim, 2);

%plot of the function 
fplot(@(alpha) sum(log(y_lim-alpha)) + n * log(2*n/power(sum(y_lim-alpha), 2))- n, [0 min(y_lim)])
grid on 

%plot of the derivative
figure
fplot(@(alpha) -sum(1./(y_lim-alpha)) + 2*n * (sum(y_lim-alpha))/(sum(power(y_lim-alpha, 2))), [0 min(y_lim)])
grid on

%% computing the optimal value of alpha based on the function

[alpha_hat, sigma_hat] = findOptimal(y_lim)

%% Iterative algorithm to compute the zero of the derivative
