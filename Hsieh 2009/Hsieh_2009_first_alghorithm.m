%C->H
clear  d uv xp Pxd Pdx Pd Phi xu Px t0 mu difU e tb tbs tbe con J X ;
close all;
A=ones(Nstate,Nstate);
xp=zeros(Nstate,Nsample); Pd(NUI,NUI,Nsample)=1;
xu=zeros(Nstate,Nsample);Px(Nstate,Nstate,Nsample)=1;
Px_pre(Nstate,Nstate,Nsample)=1;
Pxd(Nstate,NUI,Nsample)=1; Pdx (NUI,Nstate,Nsample)=1;
d = zeros(NUI,Nsample); 
Ss(Nmeas+Nstate,Nmeas,Nsample) = 1;
Phi(NUI,NUI,Nsample) = 1;
% Ft=C*G; alfa=[1 0 1 0 ; 0 1 0  1 ];
pu(:,:,1)=.1^2*eye( Nstate);
mu = 0.001;
S = [H  C*G*(eye(NUI) - pinv(H) * H)];
Te = [zeros(Nstate,NUI)  G*(eye(NUI) - pinv(H) * H)];
for k=2:1:Nsample
%% Estimation of unknown input
    Rt(:,:,k)=C*Px_pre(:,:,k)*C'+R;
    Ss(:,:,k) = pinv(S' * inv(Rt(:,:,k)) * S) * S' * inv(Rt(:,:,k));
    M(:,:,k)=[pinv(H) * H  zeros(NUI)] * Ss(:,:,k);
    Phi(:,:,k) = eye(NUI) - M(:,:,k) * H;
%  r =   rank(Rt(:,:,k));
%  r
   assump1 =  M(:,:,k) * H  % must  be equal to I
    d(:,k)=M(:,:,k)*(y(:,k)-C*xp(:,k));
    Pd(:,:,k) = M(:,:,k) * Rt(:,:,k) *M(:,:,k)' ;
%     xp(:,k)=xp(:,k)+G*d(:,k);
%% measurment update
K = Px_pre(:,:,k) * C' * inv(Rt(:,:,k));
L =  K + (Te - K * S)*Ss(:,:,k);
assump2 = L * H % must  be equal to zeros
xu(:,k) = xp(:,k) + L * (y(:,k) - C* xp(:,k) ) ;
khi = Rt(:,:,k)*L' - C * Px_pre(:,:,k);
Px(:,:,k) = Px_pre(:,:,k) - L* Rt(:,:,k)*L' + (L*khi)'; 
Pdx(:,:,k) = M(:,:,k) * khi;
Pxd(:,:,k) = Pdx(:,:,k)';

    %% time update
    xp(:,k+1) = A*xu(:,k) + B*u(k) + G*d(:,k);
    Px_pre(:,:,k+1)=[A G]*[Px(:,:,k) Pxd(:,:,k);Pdx(:,:,k) Pd(:,:,k)]*[A.' ; G.']+Q;

end
dim = 1;
temp1 = (xn(:, 1:end) - xu(:,1:end))';
RMSE_gilli = sqrt(sum(temp1.^2)/Nsample)

eGillijns2007=(x(:,1:end)-xu(:,1:end)).^2;
RMSEGillijns2007=sqrt(sum(eGillijns2007.')/(Nsample))
hold on;
title('Gillijns2007 Measurement:blue , Prediction:red , State1:green')
plot(y(1,:),'b');
plot(xu(1,:),'r');
plot(x(1,:),'g');
hold off;
figure;
hold on;
title('Gillijns2007 Measurement:blue , Prediction:red , State2:green')
plot(y(2,:),'b');
plot(xu(2,:),'r');
plot(x(2,:),'g');
hold off;