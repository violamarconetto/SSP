clear all
close all
n = 100;
y = pingstats('www.google.com', n, 'v');

y_lim = [min(y)-2 : max(y)+2];

%% Gaussian distribution
mu_ML_g = 1/n * sum(y);
sSquared_ML_g = 1/n * sum((y-mu_ML_g).^2);

f_g = (1/sqrt(2*pi*sSquared_ML_g)) * exp(-(y_lim-mu_ML_g).^2/(2*sSquared_ML_g));
%% Rayleigh distribution
sSquared_ML_r = 1/(2*n) * sum(power(y, 2));

for i = 1:size(y_lim, 2)
    if y_lim(i) < 0
        f_r(i) = 0;
    else
        f_r(i) = y_lim(i)/sSquared_ML_r * exp(-power(y_lim(i), 2)/(2*sSquared_ML_r));
    end
end

%% Erlang distribution
m = [0, 1, 2];
lambda_e = n*(m+1)/sum(y);

for i = 1:size(y_lim, 2)
    if y_lim(i) < 0
        f_e_0(i) = 0;
    else
        f_e_0(i) = lambda_e(1)^(m(1)+1)/factorial(m(1)) * y_lim(i).^m(1) .* exp(-lambda_e(1)*y_lim(i));
    end
end

for i = 1:size(y_lim, 2)
    if y_lim(i) < 0
        f_e_1(i) = 0;
    else
        f_e_1(i) = lambda_e(2)^(m(2)+1)/factorial(m(2)) * y_lim(i).^m(2) .* exp(-lambda_e(2)*y_lim(i));
    end
end

for i = 1:size(y_lim, 2)
    if y_lim(i) < 0
        f_e_2(i) = 0;
    else
        f_e_2(i) = lambda_e(3)^(m(3)+1)/factorial(m(3)) * y_lim(i).^m(3) .* exp(-lambda_e(3)*y_lim(i));
    end
end


%% Shifted exponential
alpha_exp = min(y);
lambda_exp = n/(sum(y)-n*alpha_exp);

size_y = size(y_lim, 2);
f_exp = zeros(1, size_y);

for i = 1:size_y
    if y_lim(i) < alpha_exp
        f_exp(i) = 0;
    else
        f_exp(i) = lambda_exp * exp(-lambda_exp * (y_lim(i)-alpha_exp));
    end
end

%% Shifted Rayleigh
[alpha_hat, sigma_hat] = findOptimal(y);

f_SR = zeros(1, size_y);

for i = 1:size_y
    if y_lim(i) < alpha_hat
        f_SR(i) = 0;
    else
        f_SR(i) = (y_lim(i)-alpha_hat)/sigma_hat * exp(-power(y_lim(i)-alpha_hat, 2)/(2*sigma_hat));
    end
end

%% Plots
nbins = 100;

histogram(y, nbins, 'Normalization', 'probability')
hold on
plot(y_lim, f_g);
plot(y_lim, f_r);
plot(y_lim, f_e_0);
plot(y_lim, f_e_1);
plot(y_lim, f_e_2);
plot (y_lim, f_exp);
plot (y_lim, f_SR);
legend('observations', 'gaussian', 'rayleigh', 'erlang0', 'erlang1', 'erlang2', 'shiftedExp', 'shiftedRay')
xlabel('observed values')
ylabel('f')
title('Comparison of different distributions')

%% Maximization of the log likelihood
log_L_g = -n/2*log(2*pi)-n/2*log(sSquared_ML_g)-1/(2*sSquared_ML_g)*sum(power((y-mu_ML_g), 2));
log_L_r = sum(log(y))-n*log(sSquared_ML_r)-1/(2*sSquared_ML_r)*sum(y);
log_L_e_0 = sum(-lambda_e(1)*y);
log_L_e_1 = sum(log(lambda_e(2)^(m(2)+1)/factorial(m(2)))+log(y.^(m(2)))-lambda_e(2)*y);
log_L_e_2 = sum(log(lambda_e(3)^(m(3)+1)/factorial(m(3)))+log(y.^(m(3)))-lambda_e(3)*y);
log_L_exp = n * log(lambda_exp) + n*lambda_exp*alpha_exp - sum(lambda_exp*y);
log_L_SR = sum(log(y-alpha_hat)) + n * log(1/sigma_hat) - sum((y-alpha_hat).^2/(2*sigma_hat));
% distributions = [gaussian; rayleigh; erlang_0; erlang_; erlang_2; exponential];
log_likelihoods = [log_L_g, log_L_r, log_L_e_0, log_L_e_1, log_L_e_2, log_L_exp, log_L_SR];

[argvalue, argmax] = max(log_likelihoods)
% distributions(argmax)
