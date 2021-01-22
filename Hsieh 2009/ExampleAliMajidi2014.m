%x(k+1)=A*x(k)+B*u(k)+G*d(k)+w(k)
%y(k)=C*x(k)+H*d(k)+v(k)
clear all;
Nsample=1000; Nstate=2; Nmeas=4; NUI=2; NKI=1; NUP=6;
A=[0 1; -.7 1.5]; %A(2*2)
% Ah=[.1 .8;-.5 1.3]; Bh=[1.6; 2.1];
% Ah=[1 -1; .7 -1.5]; Bh=[-1.8;-2.3];
Ah=[0 0; 0 0]; Bh=[0;0];
% B=[1 0;0 1]; %B(2*1)
B=[1.8;2.3];
G=[1 0;0 1];% F(2*2)
C=[2 1;1 0;1 2;2 1]; %H(3*3)

% matrix Dx = d
D_c = [1 -tan(pi/3); 1 0.1];
d_c = [0;0];

H=zeros(Nmeas,NUI); %=[0 0;0 1]
u=1.*randn(NKI,Nsample);
% u=ones(NKI,Nsample);
Qu=cov(u.');
% Q=[0.0036 0.;0 0.0249];
% R=[0.51 0.91;0.7 0.09];
Q=[0.0036 0.00342;0.00342 0.03249];
R=[0.51 0 0 0;0 0.26 0 0;0 0 .15 0; 0 0 0 .1];
Qnp=0.00001.*eye(NUP);
% Q=Q+Qu;
Qd=[0.025 0;0 0.016];
w=sqrt(Q)*randn(Nstate,Nsample);
v=sqrt(R)*randn(Nmeas,Nsample);
x=[1;1];
xn=[1;1];
de=zeros(NUI,Nsample);
y=zeros(Nmeas,Nsample);
de(1,200:500)=.5;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
de(2,500:700)=-.4; %de(2,700:900)=-.4;
for k=1:Nsample-1
%    de(:,k)=(A-Ah)*xn(:,k);
   x(:,k+1)=A*xn(:,k)+B*u(:,k)+G*de(:,k);
   xn(:,k+1)=x(:,k+1)+w(:,k);
   y(:,k)=C*xn(:,k)+v(:,k);
%         if abs(xn(1,k+1) - tan(pi/3) * xn(2,k+1)) > 2
%             xn(2,k+1) = (xn(2,k+1) + xn(1,k+1) * tan(pi/3)) / (1 + tan(pi/3)^2);
%             xn(1,k+1) = xn(2,k+1) * tan(pi/3);
%         end    
end
y(:,Nsample)=C*xn(:,Nsample)+H*de(:,Nsample)+v(:,Nsample);
Ub=G*(eye(NUI)-pinv(H)*H); 
S=[H C*Ub];
StableTest=eig(A);
biasTestGillijns2007=rank(C*G)-Nstate %if answer is 0 then Hsieh2006 is an Unbiased estimator for this system
biasTestHsieh2006=rank([C H])-Nmeas %if answer is 0 then Hsieh2006 is an Unbiased estimator for this system
biasTestHsieh2010=rank(S)-rank(H)-rank(Ub) %if answer is 0 then Hsieh2010 is an Unbiased estimator for this system
VarianceTestMajidiHsieh2014=eye(Nmeas)-S*pinv(S) %if answer is 0 then Hsieh2010 is Not a Minimum Variance estimator for this system 
% plot(y(1,:),'b');
% hold on;
% plot(x(1,:),'g');
% hold off;
% figure;
% plot(y(2,:),'b');
% hold on;
% plot(x(2,:),'g');
% hold off;