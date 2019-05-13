clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global ly hbar m kp sig x0 y0

hbar=1.0545718*10^-34;
m=9.10938*10^-31;

sig=8e-11;           
lamda=5.5e-11;    
kp=2*(pi/lamda)*[1,0];

dx=5e-12;
dy=5e-12;
dt=1e-21;
xf=1e-9;
yf=0.7e-9;
tf=5e-17;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/3;
y0=y(end)/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vid 1

[psy,norm] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

Psy=psy(:);
Psy_mat=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul des facteurs

% Definition potentiel
pot=1.5*10^(-15);
L_mur=1.1e-11;
D_fente=0;
L_fente=0;

v_mat=barr(x,y,pot,x(end)/2,L_mur,y(end)/2,D_fente,L_fente,'Carre');
V=v_mat(:);

% Tunnel 1: 
% v_mat=barr(x,y,3*10^(-16),x(end)/2,3e-11,y(end)/2,10e-11,0,'Carre');
% Tunnel 2:
% v_mat=barr(x,y,5*10^(-16),x(end)/2,3e-11,y(end)/2,10e-11,0,'Carre');
% 1 fente (1):
%  v_mat=barr_simple(x,y,1000,x(end)/2,3e-11,y(end)/2,30e-11,'Carre');
% 1 fente (2):
%  v_mat=barr_simple(x,y,1000,x(end)/2,3e-11,y(end)/2,8e-11,'Carre');
% 2 fentes(2):
%   v_mat=barr(x,y,1000,x(end)/2,3e-11,y(end)/2,8e-11,3e-11,'Carre');
%    v_mat=barr(x,y,1000,x(end)/2,3e-11,y(end)/2,20e-11,15e-11,'Carre');
% Inverse
%  v_mat=barr_simple_inv(x,y,1000,x(end)/2,3e-11,y(end)/2,10e-11,'Carre');
% Demi mur
% barr2=barr_simple(x,y,1000,x(end)/2,3e-11,y(1),y(end),'Carre');

[M,M2]=MM2(V,dx,dy,dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul

vid1 = VideoWriter('Effet tunnel','MPEG-4');
vid1.FrameRate = 150;
vid1.Quality = 100;
open(vid1)

graph_resultat=zeros(length(y),length(t));

tic
gcf=figure();
colormap(gcf,'hot')
for  j = 1 : length(t)-1

    tic
    b=M2*Psy;
    Psy = mldivide(M,b);
    Psy_mat=vec2mat(Psy,length(y))';
    
    surf(x,y,abs(Psy_mat).^2,'edgecolor','none');
    hold on
    surf(x,y,max(max(abs(Psy_mat).^2)).*v_mat./7000,'edgecolor','none')
    hold off
    view(0,90);
    daspect([1 1 1])
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    
    graph_resultat(:,j)=abs(Psy_mat(:,floor(2*length(x)/3),1)).^2;

    F=getframe(gcf);
    writeVideo(vid1,F);
    
    toc
    
end
toc
close(vid1)

clear M M2;

I_0=1;
%dimensions
L=y(end);
H=(2/3-0.5)*x(end);
thetamax=atan(L/(2*H));
theta=-thetamax:0.01:thetamax;

%a : distance qui separe les fentes
%b : epaisseur des fentes
alpha=(pi*D_fente/lamda).*sin(theta);
beta=(pi*L_fente/lamda).*sin(theta);
I=I_0.*(cos(alpha).*sin(beta)./beta).^2;

%% Experience

maxb=max(graph_resultat);
i=find(maxb==max(maxb));
F=graph_resultat(:,i)/max(graph_resultat(:,i));

y=-y(end)/2:dy:y(end)/2;

figure
plot(y,F,'r')
hold on
plot(H*tan(theta),I,'b');
title(sprintf('Résulats expérimentaux vs théoriques (e=%e, d=%e, a=%e)',L_mur,D_fente,L_fente));
xlabel('Position sur l''écran (m)')
ylabel('Intensité normalisée')
legend('Patron de diffraction','Fraunhofer')

