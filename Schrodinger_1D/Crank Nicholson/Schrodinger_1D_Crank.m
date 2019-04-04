clear all;
clc;
close all;

global dx C A; 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation (OK)

dx=0.01;
dt=0.0001;
t=0:dt:0.06;
x=0:dx:5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter of wave packet (OK, voir valeurs reelles)

sig=0.3;
k=50;
x_0=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation (OK)

[psy,norme]=wp_ini(x,sig,k,x_0);

norm=zeros(1,length(t));
norm(1)=norme;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse Freq et Critere stabilite

[f,Phi,STD]=fft_wp(x,psy);
f_max=f(find(Phi==max(Phi)))+4*STD;
c_max=2*pi*f_max/k;
C=dt/c_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of propagation matrix

Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potentiel

% V=zeros(1,length(x));
sig_g=0.05;
V = 1000 * exp( -((x-3.5).^2)/(2*sig_g^2) );
% V=zeros(1,length(x));
% V(350:355)=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crank nich
% Definition

a=1+1i*dt/(2*dx^2)+1i*V*dt/2;
b_C=ones(1,length(x)-1)*(1i*dt/(4*dx^2));
b_A=ones(1,length(x)-3)*(1i*dt/(4*dx^2));
g=1-1i*dt/(2*dx^2)-1i*V*dt/2;

C=sparse( full( gallery('tridiag',b_C,g,b_C)));
C=C(2:end-1,:);
C(1,1)=C(1,1)*2; C(end,end)=C(end,end)*2;

A=sparse( full(gallery('tridiag',-b_A,a(2:end-1),-b_A)) );

% Creation video
vid = VideoWriter('Schrod_Crank','MPEG-4');
vid.FrameRate = 120;
open(vid)

% Calcul
tic
gcf=figure();
for k=1:length(t)-1
    
    D = C * transpose(Psy(k,:));
    Psy(k+1,2:end-1)= mldivide(A,D);
    norm(k+1)=trapeze(abs(Psy(k+1,:)).^2,x(1),x(end),length(Psy(k+1,:))-1);
    
    plot(x, real(Psy(k,:)),'b');
    hold on
    plot(x,imag(Psy(k,:)),'r');
    hold on
    plot(x,abs(Psy(k,:)))
    hold on
    ylim([-10 10]);
    plot(x,V)
    hold off
    title( sprintf('t = %.5f , Norme= %.10f', t(k), norm(k)));

    F=getframe(gcf);
    writeVideo(vid,F);
    
end
toc

close(vid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure 
hold on
plot(t,norm)
title('Norm evolution')
xlabel('Time')
ylabel('Norm')


