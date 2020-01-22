clear all
close all
clc
%% question b
N = 256;
signal = sig(N);
N_prime = [64, 128, 256, 512, 1024];
figure;
subplot(2, 3, 1)
plot(signal)
grid on
for i=1:length(N_prime) 
    N_actual = N_prime(i);
    if N_actual < N % Why should this be done explicitely if it is already done by the function??
        subplot(2, 3, i+1);
        periodo(signal(1:N_actual), N_actual)
    else
        subplot(2, 3, i+1);
        periodo(signal, N_actual)
    end
end

%% question c
reversed_signal = signal(N:-1:1);
convolution = conv(signal, reversed_signal);
correlation = zeros(1, N);
for i=1:N
    correlation(i) = convolution(256+i-1);
end
%% question d

[error, filter, parcors] = levinson (correlation)