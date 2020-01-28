%% TP1 SSP - Riccardo SCHIAVONE and Andrea SENACHERIBBE
clear variables; close all; %clc
fprintf('TP1 SSP\n\n')
% 
% %% PART I : SPECTRUM ESTIMATION
% fprintf('Part I\n')
% %% Parameters
% N=256;
% A1=20;
% A2=20;
% f1=0.057;
% f2=0.082;
% N_prime=2.^[6:10];
% 
% %% (a)
% % (a) -> maybe at least to have a resolution to see a frequency of 0.057
% % and 0.082? so the space 0.082-0.057=0.025 so since we have a freq domain
% % of length 1, maybe it is N=1/0.025=250. The first power of 2 near 250 is
% % 256 -> N'=256
% 
% %% (b) periodogram
% y=sig(N); % signal
% w=window('boxcar',N); % rectangular signal
% signal_windowed=y.*w;
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% for i_Nprime=1:numel(N_prime)
%     subplot(numel(N_prime),1,i_Nprime)
%     periodo(signal_windowed(1:min(N_prime(i_Nprime),N)),N_prime(i_Nprime))
%     title(strcat('N=',int2str(N_prime(i_Nprime))))
% end
% sgtitle('(b) Periodogram')
% 
% % from plot it is clear that we need an N'>=N
% 
% %% (c) correlation
% correlation=@(y)conv(y,flipud(y))/numel(y);
% r=correlation(y);
% r=r(end-numel(y)+1:end);
% 
% %% (d) Levinson
% order=20;
% [variance_f,A,K]=levinson(r,order);
% A_f=fft(A,1024); % A(f)
% A_f=A_f(1:numel(A_f)/2+1);
% S_AR=variance_f(end)./abs(A_f).^2;
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% % plot periodogram + S_{AR}
% subplot(4,1,1), periodo(signal_windowed,1024),hold on
% plot(linspace(0,1/2,numel(S_AR)),10*log10(S_AR),'r');
% title('Periodogram vs S_{AR}(f)')
% legend('Periodogram','S_{AR}(f)')
% % plot A
% subplot(4,1,2), plot(0:order,A,'-o'), grid on
% xlabel('n'),ylabel('A_n'),axis([-inf,inf,1.2*min(A),1.2*max(A)])
% title('A_n prediction error filter coefficients')
% % plot variance_f
% subplot(4,1,3), plot(0:order,variance_f,'-o'), grid on
% xlabel('n'),ylabel('\sigma^2_{f,n}')
% title('\sigma^2_{f,n} prediction error variance')
% % plot K
% subplot(4,1,4), plot(K,'-o'), grid on
% xlabel('n'),ylabel('K_n'),axis([-inf,inf,1.2*min(K),1.2*max(K)])
% title('K_n PARCOR sequence')
% 
% sgtitle('(d) Levinson')

%% PART II : PARAMETER ESTIMATION
fprintf('Part II\n')

%% Parameters
realizations=1000;
variance=1;
A1=sqrt(2);
phi1=0;
f1=1/8;
n=32;

%% (e) ML estimates and CRBs
y=A1*cos(2*pi*f1*[0:n-1]+phi1)+ sqrt(variance)*randn(realizations,n);
y=y.';

% CRB
CRB_f1=6*variance/((pi^2)*(n^3)*(A1^2));
CRB_variance=2*(variance^2)/n;
CRB_A1=2*variance/n;
CRB_phi1=8*variance/(n*(A1)^2);

% using zero-padding (m>n)
m=ceil(1/sqrt(CRB_f1))+mod(ceil(1/sqrt(CRB_f1)),2);

Y=fft([y;zeros(ceil(m)-n,realizations)]);
Y=Y(1:m/2+1,:);
xaxis=[0:m/2]/m;
[~,pos]=max(abs(Y));
f1_hat=xaxis(pos.');
A1_hat=2/n*abs(Y(pos));
phi1_hat=angle(Y(pos));
y_hat=((A1_hat.').*cos(2*pi*(f1_hat.').*[0:n-1]+phi1_hat.')).';
variance_hat=1/n*sum((y-y_hat).^2);

% mean and variance of the estimates
m_variance_hat=mean(variance_hat);
v_variance_hat=var(variance_hat);
m_A1_hat=mean(A1_hat);
v_A1_hat=var(A1_hat);
m_phi1_hat=mean(phi1_hat);
v_phi1_hat=var(phi1_hat);
m_f1_hat=mean(f1_hat);
v_f1_hat=var(f1_hat);

% table
table=[variance, m_variance_hat, v_variance_hat, CRB_variance;...
    A1, m_A1_hat, v_A1_hat, CRB_A1;...
    phi1, m_phi1_hat, v_phi1_hat, CRB_phi1;...
    f1, m_f1_hat, v_f1_hat, CRB_f1];
fprintf('\n(e) ML estimates and CRBs\n')
fprintf('\t\t\ttheta\tmean\tvariance\tCRB\n')
fprintf('variance\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(1,:))
fprintf('A1\t\t\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(2,:))
fprintf('phi1\t\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(3,:))
fprintf('f1\t\t\t%.3f\t%.3f\t%.3e\t%.3e\n\n',table(4,:))

%% (f) Covariance Matching
r=zeros(2*n-1,realizations);
for i_realization=1:realizations
    r(:,i_realization)=correlation(y(:,i_realization));
end
r=r(end-n+1:end,:);
X=((r(3,:)+sqrt(r(3,:).^2+8*r(2,:).^2))./(4*r(2,:))).*(r(2,:)~=0)+zeros(size(r(2,:))).*(r(2,:)==0);
A1_hat_CM=sqrt(2*r(2,:)./X).*(r(2,:)~=0)+sqrt(-2*r(3,:)).*(r(2,:)==0);
f1_hat_CM=1/(2*pi)*acos(X).*(r(2,:)~=0)+1/4*(r(2,:)==0);
variance_hat_CM=r(1,:)-A1_hat_CM.^2/2;

% mean and variance of the estimates
m_variance_hat_CM=mean(variance_hat_CM);
v_variance_hat_CM=var(variance_hat_CM);
m_A1_hat_CM=mean(A1_hat_CM);
v_A1_hat_CM=var(A1_hat_CM);
m_f1_hat_CM=mean(f1_hat_CM);
v_f1_hat_CM=var(f1_hat_CM);

%table
table=[variance, m_variance_hat_CM, v_variance_hat_CM, CRB_variance;...
    A1, m_A1_hat_CM, v_A1_hat_CM, CRB_A1;...
    f1, m_f1_hat_CM, v_f1_hat_CM, CRB_f1];
fprintf('\n(f) Covariance Matching\n')
fprintf('\t\t\ttheta\tmean\tvariance\tCRB\n')
fprintf('variance\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(1,:))
fprintf('A1\t\t\t%.3f\t%.3f\t%.3f\t\t%.3f\n',table(2,:))
fprintf('f1\t\t\t%.3f\t%.3f\t%.3e\t%.3e\n\n',table(3,:))

