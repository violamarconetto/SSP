function [r]=sigcor(N,sigpar)
%function [r]=sigcor(N,sigpar)
%
%	Computes N lags of correlation sequence of 2 sinusoids
%	in white noise of unit variance.
%	two sinusoids: amplitudes A1, A2, frequencies f1, f2,
%	sigpar = vector of 4 signal parameters (f1,f2,A1,A2)

f1=sigpar(1);
f2=sigpar(2);
A1=sigpar(3);
A2=sigpar(4);
r=(A1^2/2)*cos(2*pi*((f1)*[0:N-1]'))+(A2^2/2)*cos(2*pi*((f2)*[0:N-1]'));
r(1)=r(1)+1;
r=r(:);	% force a column vector

