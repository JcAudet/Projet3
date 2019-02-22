clear all;
clc;
close all;

global dx dt A C; 
format long

%% Discretisation
% x=linspace(-1,1,1001);
% dx=x(2)-x(1);
% t=linspace(0,1e-4,1001);
% dt=t(2)-t(1);
x_i=0;
x_f=10;
dx=0.01;
dt=0.0001;
t=0:dt:0.03;
x=x_i:dx:x_f;

%% Stabilite
Cs=dt/dx^2

%% Parameter of wave packet
% sig=1e-10;
% k=2*pi/(0.05e-10);
sig=0.1;
k=20;
x_0=3;

%% Wave packet initialisation
psy=exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);
psy(1)=0;
psy(end)=0;

% Normalisation
psy_abs=abs(psy);
tic
norm=trapeze(psy_abs,x_i,x_f,length(psy_abs)-1);
toc
psy=psy./norm;

%% Initialisation of propagation matrix
Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

%% Potentiel
pot=zeros(1,length(x));
% sig_g=0.01;
% pot = -1000 * 1i * exp( -((x-9).^2)/(2*sig_g^2) ) - 10000 * 1i * exp( -((x-9.8).^2)/(2*sig_g^2) );
% pot = -10000 * exp( -((x-7).^2)/(2*sig_g^2) );
% pot=(x-5).^2;
% 
tic
% [C,A]=def_crank(x);
for j=2:length(t)
%     Psy(j,:)=crank_nicholson( x , Psy(j-1,:));
    Psy(j,:)=run_kutt_4( x , Psy(j-1,:), pot);
    Psy(j,1)=0;
    Psy(j,end)=0;
    j
end
toc

% tic
% figure
% % [C,A]=def_crank(x);
% for j=2:length(t)
% %     Psy(j,:)=crank_nicholson( x , Psy(j-1,:));
%     Psy(j,:)=run_kutt_4( x , Psy(j-1,:), pot);
%     j
%     
%     norm=sum(abs(Psy(j,:)))*dx;
%     plot(x, abs(Psy(j,:)),'b');
% %     ylim([0 1.5]);
%     title( sprintf('t = %.5f , Norme= %.6f', t(j), norm));
% %     pause(0.1);
% end
% toc

figure
for k = 1 : length(t)
    norm=trapeze(abs(Psy(k,:)),x_i,x_f,length(Psy(k,:))-1);
    plot(x, real(Psy(k,:)),'b');
    hold on
    plot(x,imag(Psy(k,:)),'r');
    hold on
    plot(x,abs(Psy(k,:)))
    hold on
    ylim([-10 10]);
    plot(x,pot)
    hold off
    title( sprintf('t = %.5f , Norme= %.6f', t(k), norm));
    pause(0.0000001);
end
