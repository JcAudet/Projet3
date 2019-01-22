clear all;
clc;
close all;

global dx dt; 
format long

%% Discretisation
x=linspace(0,10,500);
dx=x(2)-x(1);
t=linspace(0,1,5000);
dt=t(2)-t(1);

%% Parameter of wave packet
sig=0.3;
k=10;
x_0=5;

%% Wave packet initialisation
psy=exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);
% psy=ifft(psy);

% Normalisation
norm=sum(abs(psy));
psy=psy./norm;

%% Initialisation of propagation matrix
Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

%% Potentiel
pot=zeros(1,length(x));
% sig_g=0.01;
% pot=1000*exp(-((x-7.5).^2)/(2*sig_g^2));
% pot=(x-5).^2;


for j=2:length(t)
    Psy(j,:)=run_kutt_4( x , Psy(j-1,:) , pot);
end

figure
for k = 1 : length(t)
    norm=sum(abs(Psy(k,:)));
%     plot(x, real(Psy(k,:)),'b');
%     hold on
%     plot(x,imag(Psy(k,:)),'r');
%     hold on
    plot(x,abs(Psy(k,:)))
    hold on
    ylim([0 0.1]);
    plot(x,pot)
    hold off
    title( sprintf('t = %.5f , Norme= %.6f', t(k), norm));
    pause(0.00001);
end
