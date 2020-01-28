%% EXERCISE F -  Sinusoid in White Noise: Covariance Matching
clear all
close all
clc

%% parameters
sigma_sqr = 1;
A_1 = sqrt(2);
phi_1 = 0;
f_1 = 1/8;
n = 32;
realizations = 100;

%% y
noise_variance = 1;
noise_mean = 0;
noise = sqrt(noise_variance)*randn(realizations, n) + noise_mean;

k_vect = [0:1:n-1];
y_k = A_1*cos(2*pi*f_1*k_vect + phi_1) + noise;

%% moments
r_p = moments(A_1, f_1, noise_variance); %theoretical correlations r_p 
k = 3;
r = zeros(realizations, n); %sample correlations r
for i=1:realizations
    r(i, :) = correlation(n, y_k(i, :));
end

%% estimates of the parameters
X = zeros(1, realizations);
A1_hat_CM = zeros(1, realizations);
f1_hat_CM = zeros(1, realizations);
variance_hat_CM = zeros(1, realizations);

for i=1:realizations
    X(i) =((r(i, 3)+sqrt(r(i, 3)^2+8*r(i, 2)^2))/(4*r(i, 2)))*(r(i, 2)~=0)+zeros(size(r(i, 2)))*(r(i, 2)==0);
    A1_hat_CM(i) =sqrt(2*r(1, 2)/X(i))*(r(i, 2)~=0)+sqrt(-2*r(i, 3))*(r(i, 2)==0);
    f1_hat_CM(i) =1/(2*pi)*acos(X(i))*(r(i, 2)~=0)+1/4*(r(i, 2)==0);
    variance_hat_CM(i) =r(i, 1)-A1_hat_CM(i)^2/2;
end

%% mean and variance of the estimates
m_variance_hat_CM=mean(variance_hat_CM);
v_variance_hat_CM=var(variance_hat_CM);
m_A1_hat_CM=mean(A1_hat_CM);

%% table
table=[sigma_sqr, m_variance_hat_CM, v_variance_hat_CM, CRB_variance;...
    A_1, m_A1_hat_CM, v_A1_hat_CM, CRB_A_1;...
    f_1, m_f1_hat_CM, v_f1_hat_CM, CRB_f_1];
fprintf('\n(f) Covariance Matching\n')
fprintf('\t\t\ttheta\tmean\tvariance\tCRB\n')
fprintf('variance\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(1,:))
fprintf('A1\t\t\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(2,:))
fprintf('f1\t\t\t%.3f\t%.3f\t%.3e\t%.3e\n\n',table(3,:))
v_A1_hat_CM=var(A1_hat_CM);
m_f1_hat_CM=mean(f1_hat_CM);
v_f1_hat_CM=var(f1_hat_CM);