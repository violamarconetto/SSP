function [alpha_hat, sigma_hat] = findOptimal(y)

n = size(y, 2);
alpha_range = linspace(0, min(y), 1000);
log_likelihood_values = zeros(1, 1000);

for i=1:1000
     log_likelihood_values(i) = sum(log(y-alpha_range(i))) + n * log(2*n/(sum(y-alpha_range(i)).^2))- n;
end

[max_value , arg_max] = max (log_likelihood_values);
alpha_hat = alpha_range(arg_max);
sigma_hat = 1/(2*n)*sum(power(y-alpha_hat, 2));

format;