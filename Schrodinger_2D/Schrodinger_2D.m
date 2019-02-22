clear all;
clc;
close all;

global dx dy dt; 
format long

%% Discretisation
% x=linspace(0,2,1001);
% dx=x(2)-x(1);
% y=linspace(0,2,1001);
% dy=y(2)-y(1);
% t=linspace(0,0.3,1001);
% dt=t(2)-t(1);
dx=0.025;
dy=0.005;
dt=2.5*10e-5;
x=0:dx:1;
y=0:dy:1.5;
t=0:dt:0.002;


%% Parameter of wave packet
sig_x=0.03;
sig_y=0.06;
k=5;
k_y=1;
x_0=1;
y_0=1;

%% Wave packet initialisation
Psy=zeros(length(t),length(x),length(y));
% for k=1:length(x)
%     for j=1:length(y)
%         Psy(1,k,j)=exp(-((x(k)-x_0).^2)./(2.*sig_x.^2)).*exp(-((y(j)-y_0).^2)./(2.*sig_x.^2)).*exp(1i.*(k.*x(k)+k_y.*y(j)));
%     end
% end

Psy(1,:,:) = exp(-((x(:)-x_0).^2)./(2.*sig_x.^2)).*exp(1i.*(k.*x(:)))...
    *(exp(-((y(:)-y_0).^2)./(2.*sig_x.^2)).*exp(1i.*(k.*y(:))))';

% Normalisation
norm=sum(sum(abs(Psy(1,:,:))))*dx*dy;
Psy(1,:,:)=Psy(1,:,:)./norm;

%% Potentiel
pot=zeros(length(x),length(y));
% sig_g=0.01;
% pot = -1000 * 1i * exp( -((x-9).^2)/(2*sig_g^2) ) - 10000 * 1i * exp( -((x-9.8).^2)/(2*sig_g^2) );
% pot = -10000 * exp( -((x-7).^2)/(2*sig_g^2) );
% pot=(x-5).^2;

tic
for j=2:length(t)
    Psy(j,:,:)=run_kutt_4_2D(squeeze(Psy(j-1,:,:)) , pot);
    clc;
    j
end
toc

figure
for k = 1 : length(t)
    norm=sum(sum(abs(Psy(k,:,:))))*dx*dy;
%     plot(x, real(Psy(k,:,:)),'b');
%     hold on
%     plot(x,imag(Psy(k,:,:)),'r');
%     hold on
    surf(x,y,squeeze(abs(Psy(k,:,:))))
    hold off
%     ylim([0 1.5]);
%     plot(x,pot)
%     hold off
    title( sprintf('t = %.5f , Norme= %.6f', t(k), norm));
    pause(0.00001);
end