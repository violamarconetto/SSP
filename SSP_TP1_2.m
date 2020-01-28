%% EXERCISE E - Sinusoid in White Noise: ML Estimates and CRBs
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

%% CRB
CRB_f_1 = 6*sigma_sqr/((pi^2)*(n^3)*(A_1^2));
CRB_sigma_sqr = 2*(sigma_sqr^2)/n;
CRB_A_1 = 2*sigma_sqr/n;
CRB_phi_1= 8*sigma_sqr/(n*(A_1)^2);

%% y
noise_variance = 1;
noise_mean = 0;
noise = sqrt(noise_variance)*randn(realizations, n) + noise_mean;

k_vect = [0:1:n-1];
y_k = A_1*cos(2*pi*f_1*k_vect + phi_1) + noise;

%% m best value (m>n)
m = ceil(1/sqrt(CRB_f_1))+ mod(ceil(1/sqrt(CRB_f_1)),2);

%% ML estimates

% Initialize of the variables
sigma_sqr_hat = zeros(1, realizations);
A_1_hat = zeros(1, realizations);
phi_1_hat = zeros(1, realizations);
f_1_hat = zeros(1, realizations);
Y_f_1_hat = zeros(1, realizations);

% Compute the Fourier transform of the signal
Y_f = fft([y_k zeros(realizations, m-n)]'); %padding or the 100 realizations
Y_f = Y_f(1:m/2+1, :);
Y_f = Y_f';

for i=1:realizations
    [argvalue, argmax] = max(Y_f(i, :));
    f_1_hat(i) = argmax/m;
    Y_f_1_hat(i) = argvalue;
    A_1_hat(i) = 2/n*abs(Y_f_1_hat(i));
    phi_1_hat(i) = angle(Y_f_1_hat(i));
    y_hat = A_1_hat(i)*cos(2*pi*f_1_hat(i).*k_vect + phi_1_hat(i));
    sigma_sqr_hat(i) = sum((y_k(i, :)-y_hat).^2)/n;
end

%% Mean and variances

A_1_mean = mean(A_1_hat);
A_1_variance = var(A_1_hat);

f_1_mean = mean(f_1_hat);
f_1_variance = var(f_1_hat);

phi_1_mean = mean(phi_1_hat);
phi_1_variance = var(phi_1_hat);

sigma_sqr_mean = mean(sigma_sqr_hat);
sigma_hat_variance = var(sigma_sqr_hat);

%% table
% table
table=[sigma_sqr, sigma_sqr_mean, sigma_hat_variance, CRB_sigma_sqr;...
    A_1, A_1_mean, A_1_variance, CRB_A_1;...
    phi_1, phi_1_mean, phi_1_variance, CRB_phi_1;...
    f_1, f_1_mean, f_1_variance, CRB_f_1];
fprintf('\n(e) ML estimates and CRBs\n')
fprintf('\t\t\ttheta\tmean\tvariance\tCRB\n')
fprintf('variance\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(1,:))
fprintf('A1\t\t\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(2,:))
fprintf('phi1\t\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(3,:))
fprintf('f1\t\t\t%.3f\t%.3f\t%.3e\t%.3e\n\n',table(4,:))


%% EXERCISE F -  Sinusoid in White Noise: Covariance Matching

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
sigma_sqr_mean_CM=mean(variance_hat_CM);
sigma_hat_variance_CM=var(variance_hat_CM);
A_1_mean_CM=mean(A1_hat_CM);
A_1_variance_CM=var(A1_hat_CM);
f_1_mean_CM=mean(f1_hat_CM);
f_1_variance_CM=var(f1_hat_CM);

%% table
table=[sigma_sqr, sigma_sqr_mean_CM, sigma_hat_variance_CM, CRB_sigma_sqr;...
    A_1, A_1_mean_CM, A_1_variance_CM, CRB_A_1;...
    f_1, f_1_mean_CM, f_1_variance_CM, CRB_f_1];

fprintf('\n(f) Covariance Matching\n')
fprintf('\t\t\ttheta\tmean\tvariance\tCRB\n')
fprintf('variance\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(1,:))
fprintf('A1\t\t\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(2,:))
fprintf('f1\t\t\t%.3f\t%.3f\t%.3e\t%.3e\n\n',table(3,:))

