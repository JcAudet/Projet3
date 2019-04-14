clear all;
clc;
close all;

global dx; 
format long

hbar=1;
m=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

dx=0.01;
dt=0.00001;
t=0:dt:0.05;
x=0:dx:5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation 

% Parametres
sig=0.3;
k=20;
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
V=zeros(1,length(x));
% sig_g=0.05;
% V = 1000 * exp( -((x-3.5).^2)/(2*sig_g^2) );
% V=zeros(1,length(x));
% V(320:350)=1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runge Kutta O4

% Initialisation of propagation matrix
Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

Psy_th=zeros(length(t),length(x));
err=zeros(1,length(t));

tet=atan(2*hbar*t(1)/(m*sig^2));
delx=(sig/2)*sqrt(1+(4*hbar^2*t(1)^2)/(m^2*sig^4));
Psy_th(1,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
        exp(1i*k*(x-k*t(1)/2)).*exp(-(x-x_0-k*t(1)).^2 / (sig^2+2*1i*hbar*t(1)/m));

err(1)=trapeze(abs(abs(Psy(1,:))-abs(Psy_th(1,:))),x(1),x(end),length(x)-1);

% Creation video
% vid = VideoWriter('Schrod_RK4','MPEG-4');
% vid.FrameRate = 120;
% open(vid)

tic
% gcf=figure();
for i=1:length(t)-1
    
    Psy(i+1,:)=run_kutt_4(dt,dx,Psy(i,:),V);
    norm(i+1)=trapeze((abs(Psy(i,:))).^2,x(1),x(end),length(Psy(i,:))-1);
    
    tet=atan(2*hbar*t(i)/(m*sig^2));
    delx=(sig/2)*sqrt(1+(4*hbar^2*t(i)^2)/(m^2*sig^4));
    
    Psy_th(i+1,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
        exp(1i*k*(x-k*t(i)/2)).*exp(-(x-x_0-k*t(i)).^2 / (sig^2+2*1i*hbar*t(i)/m));
    
    err(i+1)=trapeze(abs(Psy(i,:)-Psy_th(i,:)),x(1),x(end),length(x)-1);
    
    
    plot(x, real(Psy(i,:)),'b');
    hold on
    plot(x, imag(Psy(i,:)),'r');
    hold on
    plot(x,abs(Psy(i,:)))
    hold on
    plot(x,abs(Psy_th(i,:)))
    hold on
    ylim([-10 10]);
    plot(x,V)
    hold off
    title( sprintf('t = %.5f , Norme= %.10f', t(i), norm(i)));
    pause(0.001)
    
%     F=getframe(gcf);
%     writeVideo(vid,F);
    
end
toc

figure()
plot(x,abs(Psy(1,:)))
hold on
plot(x,abs(Psy_th(1,:)))
hold on
plot(x,abs(Psy(end,:)))
hold on
plot(x,abs(Psy_th(end,:)))
hold on
plot(x,abs(Psy(floor(length(t)/2),:)))
hold on
plot(x,abs(Psy_th(floor(length(t)/2),:)))

% close(vid)

figure()
plot(t,norm-1)
hold on
plot(t,err)
legend('Norm-1','Norm-Norm theo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%