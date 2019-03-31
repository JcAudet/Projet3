clear all;
clc;
close all;

global A C; 
global C1 C2 C3 C4;
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation
x_i=0;
x_f=5;
dx=0.01;
dt=0.0001;
t=0:dt:0.3;
x=x_i:dx:x_f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter of wave packet
sig=0.3;
k=20;
x_0=2.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation
psy=exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);
psy(1)=0;
psy(end)=0;
norm=zeros(1,length(t));
% Normalisation
psy_abs=abs(psy);
norm(1)=trapeze(psy_abs,x_i,x_f,length(psy_abs)-1);
psy=psy./norm(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse Freq
[f,Phi,STD]=fft_wp(x,psy);
f_max=f(find(Phi==max(Phi)))+4*STD;
c_max=2*pi*f_max/k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Critere Stabilite

C=dt/c_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of propagation matrix
Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

%% Conditions fronti�res absorbantes (ABC's)
k1=20; % k1=k2=k donne le meilleur r�sultat ?
k2=20;
h_barre=1;
m=1;
%h_barre=1.054571800*10^(-34); vraie unit�s
%m=9.109*10^(-31);
%
%q1=h_barre*k1/m;
%q2=h_barre*k2/m;
q1=2*k1; % meilleurs q selon l'article
q2=2*k2;

%
c1=2/(q1+q2);
c2=q1*q2/(2*(q1+q2));
V=zeros(1,1000);
 
% Pour p=2;
C1=-1i/(2*(dx)) - 1i*c1/(2*dt) + (c1*V(1)-c2)/4;
C2=1i/(2*(dx)) - 1i*c1/(2*dt) + (c1*V(1)-c2)/4;
C3=1i/(2*(dx)) - 1i*c1/(2*dt) - (c1*V(1)-c2)/4;
C4=-1i/(2*(dx)) - 1i*c1/(2*dt) - (c1*V(1)-c2)/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potentiel
 V=zeros(1,length(x));
 sig_g=0.01;
% pot = -1000 * 1i * exp( -((x-9).^2)/(2*sig_g^2) ) - 10000 * 1i * exp( -((x-9.8).^2)/(2*sig_g^2) );
% V = - 100000000 * 1i * exp( -((x-4).^2)/(2*sig_g^2) );
% pot=(x-5).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition

a=1+1i*dt/(dx^2)+1i*dt*V/2;
b=ones(1,length(x)-1)*(-1i*dt/(2*dx^2));
g=1-1i*dt/(dx^2)-1i*dt*V/2; 

C=sparse( full( gallery('tridiag',-b,g,-b)));
A=sparse( full(gallery('tridiag',b,a,b)) );

%% Crank nich
tic
figure()
for k=1:length(t)
    
    Psy(k+1,:)=crank_nicholson_ABC(Psy(k,:));
    norm(k)=trapeze(abs(Psy(k,:)),x_i,x_f,length(Psy(k,:))-1);
    
    plot(x, real(Psy(k,:)),'b');
    hold on
    plot(x,imag(Psy(k,:)),'r');
    hold on
    plot(x,abs(Psy(k,:)))
    hold on
    ylim([-10 10]);
    plot(x,V)
    hold off
    title( sprintf('t = %.5f , Norme= %.6f', t(k), norm(k)));
    pause(0.01);
end
toc

% %% Runge Kutta O4
% tic
% for j=2:length(t)
%     Psy(j,:)=run_kutt_4( x , Psy(j-1,:), pot);
%     j
% end
% toc
% 
% figure
% for k = 1 : length(t)
%     norm(k)=trapeze(abs(Psy(k,:)),x_i,x_f,length(Psy(k,:))-1);
%     plot(x, real(Psy(k,:)),'b');
%     hold on
%     plot(x,imag(Psy(k,:)),'r');
%     hold on
%     plot(x,abs(Psy(k,:)))
%     hold on
%     ylim([-10 10]);
%     plot(x,pot)
%     hold off
%     title( sprintf('t = %.5f , Norme= %.6f', t(k), norm(k)));
%     pause(0.0000001);
% end


figure 
plot(t,norm)
title('Norm evolution')
xlabel('Time')
ylabel('Norm')