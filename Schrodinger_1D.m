clear all;
clc;
close all;

global dx dt; 
format long

%% Discretisation
x=linspace(0,10,1001);
dx=x(2)-x(1);
t=linspace(0,0.3,10001);
dt=t(2)-t(1);

%% Parameter of wave packet
sig=0.3;
k=100;
x_0=5;

%% Wave packet initialisation
psy=exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);
% psy=ifft(psy);

% Normalisation
norm=sum(abs(psy))*dx;
psy=psy./norm;

%% Initialisation of propagation matrix
Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

%% Potentiel
% pot=zeros(1,length(x));
sig_g=0.01;
% pot = -1000 * 1i * exp( -((x-9).^2)/(2*sig_g^2) ) - 10000 * 1i * exp( -((x-9.8).^2)/(2*sig_g^2) );
pot = -10000 * exp( -((x-7).^2)/(2*sig_g^2) );
% pot=(x-5).^2;

tic
for j=2:length(t)
    Psy(j,:)=run_kutt_2( x , Psy(j-1,:) , pot);
end
toc

figure
for k = 1 : length(t)
    norm=sum(abs(Psy(k,:)))*dx;
    plot(x, real(Psy(k,:)),'b');
    hold on
    plot(x,imag(Psy(k,:)),'r');
    hold on
    plot(x,abs(Psy(k,:)))
    hold on
    ylim([0 1.5]);
    plot(x,pot)
    hold off
    title( sprintf('t = %.5f , Norme= %.6f', t(k), norm));
    pause(0.00001);
end
% plot(x,abs(psy))