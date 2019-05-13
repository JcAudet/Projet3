clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global ly hbar m kp sig x0 y0

hbar=1.0545718*10^-34;
m=9.10938*10^-31;

sig=8e-11;           
lamda=5e-11;    
kp=2*(pi/lamda)*[1,0];

DX=1e-12


for dx=[DX DX*2 DX*4 DX*8 DX*16 DX*32 DX*64 DX*128 DX*256]
dy=dx;
dt=1e-20;
xf=2e-9;
yf=1.4e-9;
tf=1e-16;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/3;
y0=y(end)/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vid 1

[psy,norm] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

Psy=[];
Psy(:,1)=psy(:);
Psy_mat=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul des facteurs

v_mat=zeros(length(y),length(x));
V=v_mat(:);

b = 1 + 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) + 1i*dt*V/(2*hbar);
f = 1 - 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) - 1i*dt*V/(2*hbar);

c=-1i*hbar*dt*(1/dx^2)/(4*m);
d=-1i*hbar*dt*(1/dy^2)/(4*m);

g=1i*hbar*dt*(1/dx^2)/(4*m);
k=1i*hbar*dt*(1/dy^2)/(4*m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creation de M et M2, v et v2

diag=ones(1,length(V)-1)*d;
diag1=ones(1,length(V)-1)*d;
diag(ly:ly:length(V)-1)=0;
diag1(ly:ly:length(V)-1)=0;

diag2=ones(1,length(V)-1)*k;
diag3=ones(1,length(V)-1)*k;
diag2(ly:ly:length(V)-1)=0;
diag3(ly:ly:length(V)-1)=0;

M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
    [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
    ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c diag diag1]);
whos M
end


figure(1)
plot(MemoireETdiscretisation(:,1),MemoireETdiscretisation(:,2));
title('Memoire utilisee par la matrice M en fonction de la discretisation');
ylabel('Memoire [bytes]');
xlabel('Discretisation dx [m]');

P_x=polyfit(log10(MemoireETdiscretisation(:,1)),log10(MemoireETdiscretisation(:,2)),1)


figure(2)
loglog(MemoireETdiscretisation(:,1),MemoireETdiscretisation(:,2),'b','linewidth',1.5);
hold on
loglog(MemoireETdiscretisation(:,1),10.^(P_x(1)*log10(MemoireETdiscretisation(:,1)) +P_x(2)),'r--','LineWidth',1.5)

title('Memoire utilisee par la matrice M en fonction de la discretisation');
ylabel('Memoire [bytes]');
xlabel('Discretisation dx [m]');
legend('Demande en mémoire expérimental','Régression linéaire (y=-1.9877-15.3041)')
