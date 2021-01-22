%x(k+1)=A*x(k)+B*u(k)+G*d(k)+w(k)
%y(k)=C*x(k)+H*d(k)+v(k)
clear all;
Nsample=1000; Nstate=4; Nmeas=4; NUI=4; NKI=1; NUP=6; T =3;tantheta=tan(pi/3);
A=[1 0 T 0 
    0 1 0 T
    0 0 1 0
    0 0 0 1]; %A(4*4)

B=0;
G=[0.1 0 0 0;0.1 1 0 0;0 0 1 -.01 ; 1 0 0 0.5];% F(4*4)
C=[2 1 0 0;0 1 0 0;0 0 1 0;0 0 0 1]; %H(3*3)

% matrix Dx = d
% D_c = [1 -tan(pi/3) 0 0 
%            0 0 1 -tan(pi/3)];
% d_c = [0;0];

 H=[1 0 0 0;0 1 0 0;0 0 1 0 ; 0 0 0 1]; %=[0 0;0 1]
u=1.*randn(NKI,Nsample);
% u=ones(NKI,Nsample);

Q=[0.0036 0 0 0;
   0  0.00342 0 0;
   0 0 0.0001 0 
   0 0 0 0.0001];
R=[0.51 0 0 0;
    0 0.26 0 0;
    0 0 .15 0;
    0 0 0 .1];
w=sqrt(Q)*randn(Nstate,Nsample);
v=sqrt(R)*randn(Nmeas,Nsample);
x=[1;1;1;1];
xn=[1;1;1;1];
de=zeros(NUI,Nsample);
y=zeros(Nmeas,Nsample);
de(1,200:500)=.3;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
de(2,300:500)=.1;  de(2,100:300)=0.2;
de(3,600:700)=-.4; de(3,800:900)=-.2;
de(4,100:700)=-.1; de(4,900:1000)=0.3;

for k=1:Nsample-1
%     de(:,k)=(A-Ah)*xn(:,k);
    x(:,k+1)=A*xn(:,k)+B*u(:,k)+G*de(:,k);
    xn(:,k+1)=x(:,k+1)+w(:,k);
    y(:,k)=C*xn(:,k)+H*de(:,k)+v(:,k);
    % Constrain the vehicle (i.e., the true state) to the straight road.
%     if abs(xn(1,k+1) - tantheta * xn(2,k+1)) > 2
%         xn(2,k+1) = (xn(2,k+1) + xn(1,k+1) * tantheta) / (1 + tantheta^2);
%         xn(1,k+1) = xn(2,k+1) * tantheta;
%     end
%     if abs(xn(3,k+1) - tantheta * xn(4,k+1)) > 0.2
%         xn(4,k+1) = (xn(4,k+1) + xn(3,k+1) * tantheta) / (1 + tantheta^2);
%         xn(3,k+1) = xn(4,k+1) * tantheta;
%     end
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