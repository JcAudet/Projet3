clear all;
clc;
close all;

global dx C A sig x_0 k hbar m; 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation (OK)

hbar=1;
m=1;

dx=0.01;
dt=0.00001;
t=0:dt:0.05;
x=0:dx:5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation 

% Parametres
sig=0.2;
k=50;
x_0=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation (OK)

[psy1,norme1]=wp_ini(x,sig,k,x_0);
[psy2,norme2]=wp_ini(x,sig*2,k,x_0);
[psy3,norme3]=wp_ini(x,sig/2,k,x_0);

norm1=zeros(1,length(t));
norm1(1)=norme1;
norm2=zeros(1,length(t));
norm2(1)=norme2;
norm3=zeros(1,length(t));
norm3(1)=norme3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of propagation matrix

Psy1=zeros(length(t),length(x));
Psy1(1,:)=psy1;
Psy2=zeros(length(t),length(x));
Psy2(1,:)=psy2;
Psy3=zeros(length(t),length(x));
Psy3(1,:)=psy3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potentiel

V=zeros(1,length(x));

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


vid = VideoWriter('Heisenberg_1D_RK4','MPEG-4');
vid.FrameRate = 150;
vid.Quality = 100;
open(vid)

% Calcul
tic
gcf=figure();
for j=1:length(t)-1
    
    tic
    
    Psy1(j+1,:)=run_kutt_4(dt,dx,Psy1(j,:),V);
    Psy2(j+1,:)=run_kutt_4(dt,dx,Psy2(j,:),V);
    Psy3(j+1,:)=run_kutt_4(dt,dx,Psy3(j,:),V);
        
    plot(x,abs(Psy1(1,:)).^2,'b--')
    hold on
    plot(x,abs(Psy2(1,:)).^2,'r--')
    hold on
    plot(x,abs(Psy3(1,:)).^2,'g--')
    
    if j>floor(length(t)/2)
            plot(x,abs(Psy1(floor(length(t)/2),:)).^2,'b--')
            hold on
            plot(x,abs(Psy2(floor(length(t)/2),:)).^2,'r--')
            hold on
            plot(x,abs(Psy3(floor(length(t)/2),:)).^2,'g--')
    end
    
    plot(x,abs(Psy1(j,:)).^2,'b','linewidth',1)
    hold on
    plot(x,abs(Psy2(j,:)).^2,'r','linewidth',1)
    hold on
    plot(x,abs(Psy3(j,:)).^2,'g','linewidth',1)
    
    ylabel('|\Psi|^2')
    xlabel('Position (u.a)')
    
    hold off
    
    F=getframe(gcf);
    writeVideo(vid,F);
    
    toc           
end
toc

close(vid)