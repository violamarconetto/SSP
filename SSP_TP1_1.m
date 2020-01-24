clear all
close all
clc
%% question b
N = 256;
s=sig(N);
w = window('boxcar', N);
signal = s.*w;
N_prime = [64, 128, 256, 512, 1024];
figure;
subplot(2, 3, 1)
plot(s)
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
reversed_signal = s(N:-1:1);
convolution = conv(s, reversed_signal)/numel(s);
correlation = zeros(1, N);
for i=1:N
    correlation(i) = convolution(256+i-1);
end
%% question d

[filter, error, parcors] = myLevinson (correlation);
%[a,b,c]=levinson(correlation,20);

N = 1024
filter_f = fft(filter, N);
filter_f = filter_f(1:(N/2)+1);
f = 0:1/N:0.5;
AR_spectrum = error(end)./(abs(filter_f).^2);
%AR_spectrum = b./(abs(fft(a,N)).^2);

figure;

subplot(4, 1, 1);
plot(f, 10*log10(AR_spectrum))
hold on;
periodo(signal, 1024)

subplot(4, 1, 2);
plot(0:20, filter)
axis([0,20, -inf, inf])
title('Filters coefficients')
xlabel('')
ylabel('A(0:20)')
grid on

subplot(4, 1, 3);
plot(0:20, error)
axis([0,20, -inf, inf])
title('Evolution of the error')
ylabel('sigma^2(0:20)')
grid on

subplot(4, 1, 4);
plot(1:20, parcors)
axis([1,20, -inf, inf])
title('Parcors')
ylabel('K(0:20)')
grid on

sgtitle('Riccardo I love you')






