%x(k+1)=A*x(k)+B*u(k)+G*d(k)+w(k)
%y(k)=C*x(k)+H*d(k)+v(k)
clear all;
Nsample=1000; Nstate=2; Nmeas=4; NUI=3; NKI=1; NUP=6;
A=[0 1; -.7 1.5]; B=[1.8;2.3]; 

G=[1 0 0   ; 0  1 0];% F(2*2)
C=[2 1;1 0;1 2;2 1]; %H(3*3)

u=1.*randn(NKI,Nsample);
H = [1 0 0; 0 1 0;0 1 1;0.1 0 0];
Q=diag([0.001,0.001]);
R=[9 0 0 0;0 9 0 0;0 0 9 0; 0 0 0 9];

w=sqrt(Q)*randn(Nstate,Nsample);
v=sqrt(R)*randn(Nmeas,Nsample);
x=[1;1];
xn=[1;1];
de=zeros(NUI,Nsample);
y=zeros(Nmeas,Nsample);
de(1,200:500)=.5;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
de(2,500:700)=-.4; de(2,800:100)=0.2;
de(3,100:400)=-.0; de(3,800:100)=.0;
for k=1:Nsample-1
%    de(:,k)=(A-Ah)*xn(:,k);
   x(:,k+1)=A*xn(:,k)+B*u(:,k)+G*de(:,k);
   xn(:,k+1)=x(:,k+1)+w(:,k);
   y(:,k)=C*xn(:,k)+H*de(:,k)+v(:,k);
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