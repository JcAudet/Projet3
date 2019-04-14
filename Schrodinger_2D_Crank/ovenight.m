clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global ly hbar m kp sig x0 y0


hbar=1.0545718*10^-34;
m=9.10938*10^-31;

sig=8e-11;           
lamda=5e-11;    
kp=2*(pi/lamda)*[1,0];

Imax=10;
Jmax=10;

dx=0.6e-11;
dy=0.6e-11;
dt=1e-20;
xf=1.8e-9;
yf=1.3e-9;
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

v_mat=barr(x,y,1000,x(end)/2,3e-11,y(end)/2,10e-11,5e-11,'Carre');
% v_mat=barr(x,y,1000,x(end)/2,3e-11,y(end)/2,20e-11,15e-11,'Carre');
% v_mat=barr_simple(x,y,1000,x(end)/2,3e-11,y(end)/2,40e-11,'Carre');
% v_mat=barr_simple_inv(x,y,1000,x(end)/2,3e-11,y(end)/2,10e-11,'Carre');
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
M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
    [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
    ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g diag2 diag3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul

vid1 = VideoWriter('Schrod_Crank_Diffra_1','MPEG-4');
vid1.FrameRate = 150;
vid1.Quality = 100;
open(vid1)

tic
gcf=figure();
colormap(gcf,'gray')
xlim([x(1) x(end)])
ylim([y(1) y(end)])
for  j = 1 : length(t)-1

    tic
    b=M2*Psy;
    Psy = mldivide(M,b);
    Psy_mat=vec2mat(Psy,length(y))';
    toc
    
    surf(x,y,abs(Psy_mat).^2,'edgecolor','none');
    hold on
    surf(x,y,max(max(abs(Psy_mat).^2))*v_mat/7000,'edgecolor','none')
    hold off
    view(0,90);
    daspect([1 1 1])
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    
    F=getframe(gcf);
    writeVideo(vid1,F);
    
end
toc
close(vid1)

clear M M2;
close all;

