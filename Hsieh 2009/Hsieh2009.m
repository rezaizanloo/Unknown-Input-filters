xpg=zeros(Nstate,1); Pxg=zeros(Nstate,Nstate); n=length(y);
for k=1:n

gama=[zeros(Nstate,NUI) G*(eye(NUI)-pinv(H)*H)]; %(34)
S=[H C*G*(eye(NUI)-pinv(H)*H)]; %(Nstate,2*NUI)(39)
Rt(:,:,k)=C*Pxg(:,:,k)*C.'+R; %(3)
Sstar(:,:,k)=pinv(S.'*inv(Rt(:,:,k))*S)*S.'*inv(Rt(:,:,k));%(2*NUI,Nstate)
Mstar(:,:,k)=[pinv(H)*H zeros(NUI)]*Sstar(:,:,k); % (4)
dp(:,k)=Mstar(:,:,k)*(y(k)-C*xpg(:,k)); %(17)
Pd(:,:,k)=Mstar(:,:,k)*Rt(:,:,k)*Mstar(:,:,k).'; %(29)

K(:,:,k)=Pxg(:,:,k)*C.'*inv(Rt(:,:,k)); %(7)
Lstar(:,:,k)=K(:,:,k)+(gama-K(:,:,k)*S)*Sstar(:,:,k);%(34)
xp(:,k)=xpg(:,k)+Lstar(:,:,k)*(y(k)-C*xpg(:,k)); %(18)
Fi(:,:,k)=Rt(:,:,k)*Lstar(:,:,k).'-C*Pxg(:,:,k); %(32)
Px(:,:,k)=Pxg(:,:,k)-Lstar(:,:,k)*Rt(:,:,k)*(Lstar(:,:,k))'+Lstar(:,:,k)*Fi(:,:,k)+(Lstar(:,:,k)*Fi(:,:,k)).'; %(30)
Pdx(:,:,k)=Mstar(:,:,k)*Fi(:,:,k); %(10)
Pxd(:,:,k)=Pdx(:,:,k).'; %(10)

xpg(:,k+1)=A*xp(:,k)+G*dp(:,k)+B*u(k); %(11)
Pxg(:,:,k+1)=[A G]*[Px(:,:,k) Pxd(:,:,k);Pdx(:,:,k) Pd(:,:,k)]*[A.';G.']+Q; % (Nstate,Nstate)(12)
eHsieh2009(:,k)=(x(:,k)-xpg(:,k)).^2; % Error 

end
rmseHsieh2009=sqrt(sum(eHsieh2009));
close all
hold on;
title('Hsieh2009 Measurement:blue , Prediction:red , State1:green')
plot(y(1,:),'b');
plot(xpg(1,:),'r');
plot(x(1,:),'g');
hold off;
figure;
hold on;
title('Hsieh2009 Measurement:blue , Prediction:red , State2:green')
plot(y(2,:),'b');
plot(xpg(2,:),'r');
plot(x(2,:),'g');
hold off;
% figure;
% hold on;
% title('Hsieh2009 Measurement:blue , Prediction:red , State3:green')
% plot(y(3,:),'b');
% plot(xpg(3,:),'r');
% plot(x(3,:),'g');
% hold off;