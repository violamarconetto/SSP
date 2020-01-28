function [r]=correlation(N, signal)
reversed_signal = signal(N:-1:1);
convolution = conv(signal, reversed_signal)/numel(signal);
r = zeros(1, N);
for i=1:N
    r(i) = convolution(N+i-1);
end