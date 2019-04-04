clear all;
clc;
close all;

global dx; 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

dx=0.01;
dt=0.00005;
t=0:dt:0.07;
x=0:dx:5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation 

% Parametres
sig=0.3;
k=40;
x_0=2;

% Creation
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
%% Potentiel
% V=zeros(1,length(x));
sig_g=0.05;
V = 1000 * exp( -((x-3.5).^2)/(2*sig_g^2) );
% V=zeros(1,length(x));
% V(320:350)=1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runge Kutta O4

% Initialisation of propagation matrix
Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

% Creation video
vid = VideoWriter('Schrod_RK4','MPEG-4');
vid.FrameRate = 120;
open(vid)

tic
gcf=figure();
for k=1:length(t)-1
    
    Psy(k+1,:)=run_kutt_4(dt,dx,Psy(k,:),V);
    norm(k+1)=trapeze((abs(Psy(k,:))).^2,x(1),x(end),length(Psy(k,:))-1);
    
    plot(x, real(Psy(k,:)),'b');
    hold on
    plot(x, imag(Psy(k,:)),'r');
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

figure()
plot(t,norm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%