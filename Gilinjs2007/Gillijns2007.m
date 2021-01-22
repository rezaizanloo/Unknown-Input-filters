%C->H
clear  d uv xp pp xu pu t0 mu difU e tb tbs tbe con J X ;
close all;
Examples
A=ones(Nstate,Nstate);
xp=zeros(Nstate,Nsample); pp(Nstate,Nstate,Nsample)=0;
xu=zeros(Nstate,Nsample); pu(Nstate,Nstate,Nsample)=0;
Ft=C*G; alfa=[0.1 .01; 1 .1 ];
pu(:,:,1)=.1^2*eye( Nstate);
mu = 0.01;
D_c = eye(size(D,1));
for k=2:1:Nsample
   
    %% time update
    xp(:,k)=A*xu(:,k-1);
    pp(:,:,k)=A*pu(:,:,k-1)*A.'+Q;
    %% Estimation of unknown input
    Rt(:,:,k)=C*pp(:,:,k)*C.'+R;
    M(:,:,k)=inv(Ft.'*inv(Rt(:,:,k))*Ft)*Ft.'*inv(Rt(:,:,k));
    M(:,:,k) * C * G
    d(:,k-1)=M(:,:,k)*(y(:,k)-C*xp(:,k));
    xps(:,k)=xp(:,k)+G*d(:,k-1);
     %% measurment update
    pps(:,:,k)=(eye(Nstate)-G*M(:,:,k)*C)*pp(:,:,k)*(eye(Nstate)-G*M(:,:,k)*C).'+G*M(:,:,k)*R*(G*M(:,:,k)).';
    Ss(:,:,k)=-1.*G*M(:,:,k)*R;
%     sigma = 0.000001;
%     Rts(:,:,k)=(eye(Nmeas)-C*G*M(:,:,k))*Rt(:,:,k)*(eye(Nmeas)-C*G*M(:,:,k)).' + sigma * eye(Nstate);
    Rts(:,:,k)=(eye(Nmeas)-C*G*M(:,:,k))*Rt(:,:,k)*(eye(Nmeas)-C*G*M(:,:,k)).' ;
    Gain(:,:,k)=(pps(:,:,k)*C.'+Ss(:,:,k))*alfa.'*inv(alfa*Rts(:,:,k)*alfa.')*alfa;
    xu(:,k) = xps(:,k)+Gain(:,:,k)*(y(:,k)-C*xps(:,k));
    pu(:,:,k) = pps(:,:,k)-Gain(:,:,k)*(pps(:,:,k)*C.'+Ss(:,:,k)).';
end
dim = 1;
% temp1 = (xn(:, 1:end) - xu(:,1:end))';
% RMSE_gilli = sqrt(sum(temp1.^2)/Nsample)

eGillijns2007=(xn(:,1:end)-xu(:,1:end));
RMSEGillijns2007=sqrt(sum(eGillijns2007.').^2/(Nsample))
hold on;
title('Gillijns2007 Measurement:blue , Prediction:red , State1:green')
plot(y(1,:),'b');
plot(xu(1,:),'r');
plot(xn(1,:),'g');
hold off;
figure;
hold on;
title('Gillijns2007 Measurement:blue , Prediction:red , State2:green')
plot(y(2,:),'b');
plot(xu(2,:),'r');
plot(xn(2,:),'g');
hold off;