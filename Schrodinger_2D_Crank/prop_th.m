clear all 
clc
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 1D
% 
% dx=0.01;
% dt=0.00001;
% t=0:dt:0.05;
% x=0:dx:5;
% 
% 
% sig=0.3;
% % hbar=6.62607*10^-34;
% % m=9.10938*10^-31;
% hbar=1;
% m=1;
% k=20;
% x_0=2;
% 
% Psy=zeros(length(t),length(x));
% norm=zeros(1,length(t));
% 
% figure()
% for i=1:length(t)
%     
%     tet=atan(2*hbar*t(i)/(m*sig^2));
%     delx=(sig/2)*sqrt(1+(4*hbar^2*t(i)^2)/(m^2*sig^4));
%     
%     Psy(i,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
%         exp(1i*k*(x-k*t(i)/2)).*exp(-(x-x_0-k*t(i)).^2 / (sig^2+2*1i*hbar*t(i)/m));
%     
%     norm(i)=trapeze(abs(Psy(i,:)).^2,x(1),x(end),length(x)-1);
%     
%     plot(x,abs(Psy(i,:)).^2)
%     title(sprintf('t=%d , norm=%d',t(i),norm(i)))
%     hold off
%     
%     pause(0.01)
%     
% end
% 
% 
% figure()
% plot(t,norm)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D

hbar=6.62607*10^-34;
m=9.10938*10^-31;

dx=1e-11;
dy=1e-11;
dt=1e-19;
x=0:dx:1e-9;
y=0:dy:1e-9;
t=0:dt:1e-17;


sig_x=1e-10;                 % Largeur initiale du packet d'onde
sig_y=1e-10;
lamda=1e-10;    % Accelerer du repos avec 50KeV
kp=2*(pi/lamda)*[1 0];

Psy_x=zeros(length(t),length(y));
Psy_y=zeros(length(t),length(x));
Psy=zeros(length(t),length(y),length(x));

norm=zeros(1,length(t));

figure()
for j=1:length(t)
    tet_x=atan(2*hbar*t(j)/(m*sig_x^2));
    delx=(sig_x/2)*sqrt(1+(4*hbar^2*t(j)^2)/(m^2*sig_x^4));

    Psy_x(j,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet_x/2).*...
    exp(1i*kp(1)*(x-kp(1)*hbar*t(j)/(2*m))).*exp(-(x-x(end)/4-hbar*kp(1)*t(j)/m).^2 / (sig_x^2+2*1i*hbar*t(j)/m));

    tet_y=atan(2*hbar*t(j)/(m*sig_y^2));
    dely=(sig_y/2)*sqrt(1+(4*hbar^2*t(j)^2)/(m^2*sig_y^4));

    Psy_y(j,:)=(1/sqrt(sqrt(2*pi)*dely))*exp(-1i*tet_y/2).*...
    exp(1i*kp(2)*(y-kp(2)*hbar*t(j)/(2*m))).*exp(-(y-y(end)/2-hbar*kp(2)*t(j)/m).^2 / (sig_y^2+2*1i*hbar*t(j)/m));

    Psy(j,:,:)=Psy_y(j,:,:)'*Psy_x(j,:,:);
    
    norm(j)=trapeze_2D(abs(squeeze(Psy(j,:,:))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    
    surf(x,y,squeeze(abs(Psy(j,:,:)).^2))
    title(sprintf('norm=%.10f',norm(j)))
    view(0,90)
    hold off
   
    pause(0.01)
    
end

