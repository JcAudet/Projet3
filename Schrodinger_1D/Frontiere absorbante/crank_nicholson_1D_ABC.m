clear all;
clc;
close all;

global A C; 
global C1 C2 C3 C4;
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

dx=0.01;
dt=0.0001;
t=0:dt:0.03;
x=0:dx:4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter of wave packet
sig=0.3;
k=50;
x_0=2.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation

norm=zeros(1,length(t));
Psy=zeros(length(t),length(x));

[Psy(1,:),norm(1)]=wp_ini(x,sig,k,x_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse Freq
% [f,Phi,STD]=fft_wp(x,psy);
% f_max=f(find(Phi==max(Phi)))+4*STD;
% c_max=2*pi*f_max/k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Critere Stabilite

% C=dt/c_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of propagation matrix

%% Conditions frontières absorbantes (ABC's)
k1=k; % k1=k2=k donne le meilleur résultat ?
k2=k;
h_barre=1;
m=1;
%h_barre=1.054571800*10^(-34); vraie unités
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
%  sig_g=0.01;
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

% Frontieres absorbante
A(1,1)=C1;      % Conditions frontières absorbantes
A(1,2)=C2;      % Valeurs à changer dans les matrices à résoudre
A(end,end)=C1;
A(end-1,end)=C2;

C(1,1)=C3;
C(1,2)=C4;
C(end,end)=C3;
C(end-1,end)=C4;

% Creation video
vid = VideoWriter('Schrod_Crank_abso','MPEG-4');
vid.FrameRate = 60;
open(vid)

%% Crank nich
tic
gcf=figure();
for k=1:length(t)-1
    
    b= C * transpose(Psy(k,:));
    Psy(k+1,:) = A \ b;
    
    norm(k)=trapeze(abs(Psy(k,:)).^2,x(1),x(end),length(Psy(k,:))-1);
    
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

figure 
plot(t,norm)
title('Norm evolution')
xlabel('Time')
ylabel('Norm')