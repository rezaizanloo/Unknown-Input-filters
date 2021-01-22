
%x(k+1)=A*x(k)+B*u(k)+G*d(k)+w(k)
%y(k)=C*x(k)+H*d(k)+v(k)
clear all;
%% system 1
% Nsample=1000; Nstate=5; Nmeas=5; NUI=3; NKI=1; NUP=6;
% A=[0.5 2 0 0 0 ; 0 0.2 1 0 1;0 0 0.3 0 1; 0 0 0 0.7 1;0 0 0 0 0.1];
% G=[1 0 -.3;1 0 0;0 0 0;0 0 0;0 0 0];
% C=eye(Nmeas);
% H=[0 0 1;0 0 0;0 2 0; 0 0 0;1 0 0];
% Q=0.0001*[1 0 0 0 0;0 1  0.5 0 0; 0 0.5 1 0 0;0 0 0 1 0;0 0 0 0 1];
% R=0.01 * [1 0 0 0.5 0;0 1 0  0 0.3; 0 0 1 0 0; 0.5 0  0 1 0; 0 0.3 0 0 1];
% w=sqrt(Q)*randn(Nstate,Nsample);
% v=sqrt(R)*randn(Nmeas,Nsample);
% de=zeros(NUI,Nsample);
% y=zeros(Nmeas,Nsample);
% de(1,200:500)=.2;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
% de(2,500:700)=-.3; de(2,800:100)=0.2;
% de(3,500:700)=-.1; de(3,800:100)=-0.7;
% x=[1;1;1;1;1];
% xn=[1;1;1;1;1];
%% system 2
% Nsample=1000; Nstate=2; Nmeas=4; NUI=2; NKI=1; NUP=6;
% A=[0 1; -.7 1.5]; %A(2*2)
% G=[1 0;0 1];% F(2*2)
% C=[2 1;1 0;1 2;2 1]; %H(3*3)
% H=[1 0;0 1;0.1 0;1 .5]; %=[0 0;0 1]
% Q=[0.0036 0.00342;0.00342 0.03249];
% R=[0.51 0 0 0;0 0.26 0 0;0 0 .15 0; 0 0 0 .1];
% w=sqrt(Q)*randn(Nstate,Nsample);
% v=sqrt(R)*randn(Nmeas,Nsample);
% de=zeros(NUI,Nsample);
% y=zeros(Nmeas,Nsample);
% de(1,200:500)=.5;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
% de(2,500:700)=-.4; %de(2,700:900)=-.4;
% x=[1;1];
% xn=[1;1];
% % D_c = [1 -tan(pi/3); 1 0.1];
% % a_c = [0;0];
%% system 3
    % Nsample=1000; Nstate=4; Nmeas=2; NUI=2; NKI=1; NUP=6; T =3;tantheta=tan(pi/3);
    % A=[1 0 T 0 ; 0 1 0 T;  0 0 1 0; 0 0 0 1]; %A(4*4)
    % G= [1 0;1 0 ;0 0 ;0 0 ];% F(4*4)
    % C=[1 0 0 0;0 1 0 0]; %H(4*4)
    %  H=[1 0; 0 0]; %
    % Q=[4 0 0 0;0  4 0 0; 0 0 1 0 ;0 0 0 1];
    % R=[900 0 ;0 900 ];
    % w=sqrt(Q)*randn(Nstate,Nsample);
    % v=sqrt(R)*randn(Nmeas,Nsample);
    % de=zeros(NUI,Nsample);
    % y=zeros(Nmeas,Nsample);
    % de(1,200:500)=.3;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
    % de(2,300:500)=.1;  de(2,100:300)=0.2;
    % % de(3,600:700)=-.4; de(3,800:900)=-.2;
    % 
    % x=[1;1;1;1];
    % xn=[1;1;1;1];
    % % matrix Dx = d
    % % D_c = [1 -tan(pi/3) 0 0 ; 0 0 1 -tan(pi/3)];
    % % d_c = [0;0];
%% system 4
% Nsample=1000; Nstate=3; Nmeas=3; NUI=2;
% % original system
% A=[0.9944 -0.1203 -0.4302; 0.0017 0.9902 -0.0747;0 0.8187 0];
% G=[1 0;0 1;0.5 0];
% C=[1 0 0;0 1 0;0 0 1];
% H=[1 0 ;0 1; 0 0.5];
% Q=[0.01 0 0;0 0.01 0; 0 0 0.0001];
% R=[1 0 0;0 1 0; 0 0 1];
%     w=sqrt(Q)*randn(Nstate,Nsample);
%     v=sqrt(R)*randn(Nmeas,Nsample);
% x=zeros(Nstate,Nsample);x0=[1;1;1];
% y=zeros(Nmeas,Nsample);
% de(1,200:500)=.2;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
% de(2,500:700)=-.3; de(2,800:100)=0.2;
% 
% x=[1;1;1];
% xn=[1;1;1];

%% system 5
Nsample=1000; Nstate=2; Nmeas=2; NUI=2; NKI=1; NUP=6;
A=[-0.0005  -0.0084; 0.0517  0.8069]; 
B=[0.1815;1.7902]; 
G=[0.0129 0    ; -1.2504  0 ];% F(2*2)
C=[1 0;0 1]; %H(3*3)
H = [0 0 ; 0 1];
Q=[0.0036 0.0342;0.0342 0.3249];
R=[0.01 0;0 0.16 ];
w=sqrt(Q)*randn(Nstate,Nsample);
v=sqrt(R)*randn(Nmeas,Nsample);
de=zeros(NUI,Nsample);
y=zeros(Nmeas,Nsample);
de(1,200:500)=.5;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
de(2,500:700)=-.4; de(2,800:100)=0.2;
x=[1;1];
xn=[1;1];

%% system 6
% Nsample=1000; Nstate=4; Nmeas=4; NUI=4; NKI=1; NUP=6; T =1/2;tantheta=tan(pi/3);
% A=[1 0 T 0 ; 0 1 0 T;  0 0 1 0; 0 0 0 1]; %A(4*4)
% G= [1 0 0 0;0 1 0 0 ;0 0 1 0 ;0 0 0 1 ];% F(4*4)
% C=[1/2 0 0 0;0 1 0.5 0 ;0 0.5 1 0;0 0 0 0]; %H(4*4)
% H=[1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1 ]; %
% Q=0.001* [0.2 0 0 0;0  0.4 0.2 0; 0 0.2 0.1 0 ;0 0 0 .1];
% R=0.01 * [9 0 0 0 ;0 9 0 0; 0 0 1 0;0 0 0 1 ];
% w=sqrt(Q)*randn(Nstate,Nsample);
% v=sqrt(R)*randn(Nmeas,Nsample);
% de=zeros(NUI,Nsample);
% y=zeros(Nmeas,Nsample);
% de(1,200:500)=.3;  de(1,700:1000)=0.0;%de(1,1:Nsample)=1; %
% de(2,300:500)=.1;  de(2,100:300)=0.2;
% de(3,600:700)=-.4; de(3,800:900)=-1.9;
% de(4,100:700)=-5; de(4,900:1000)=0.3;
% x=[1;1;1;1];
% xn=[1;1;1;1];
teta = pi/3;
  D =  [1 -1 0 1
             1 0 -1 1];
    a =  [1;0 ];
for k=2:Nsample
%    de(:,k)=(A-Ah)*xn(:,k);
 xn(:,k)=A*xn(:,k-1)+G*de(:,k-1)+w(:,k-1);
%            x(:,k) = A * x(:,k) +  Qsqrt*randn(size(temp));
         % Constrain the vehicle (i.e., the true state) to the straight road.
%    if abs(x(1,k) - tan(teta) * x(2,k)) > 2
%       x(2,k) = (x(2,k) + x(1,k) * tan(teta)) / (1 + tan(teta)^2);
%       x(1,k) = x(2,k) * tan(teta);
%    end
%    if abs(x(3,k) - tan(teta) * x(4,k)) > 0.2
%       x(4,k) = (x(4,k) + x(3,k) * tan(teta)) / (1 + tan(teta)^2);
%       x(3,k) = x(4,k) * tan(teta);
%    end
   y(:,k-1)=C*xn(:,k-1)+H*de(:,k-1)+v(:,k-1);
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