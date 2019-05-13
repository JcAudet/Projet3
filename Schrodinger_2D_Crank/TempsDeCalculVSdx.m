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

dx=10e-12;
dy=10e-12;
time=zeros(1,8);
indice=1;

DX=1e-12


for dx=[DX DX*2 DX*4 DX*8 DX*16 DX*32 DX*64 DX*128 DX*256]
tic
dy=dx;
dt=1e-17;
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

Psy=psy(:);
Psy_mat=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul des facteurs

% Definition potentiel
pot=1000;
L_mur=3e-11;
D_fente=10e-11;
L_fente=5e-11;

v_mat=barr(x,y,pot,x(end)/2,L_mur,y(end)/2,D_fente,L_fente,'Carre');
V=v_mat(:);

[M,M2]=MM2(V,dx,dy,dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul

vid1 = VideoWriter('2 Fentes (fentes larges)','MPEG-4');
vid1.FrameRate = 150;
vid1.Quality = 100;
open(vid1)

graph_resultat=zeros(length(y),length(t));

%tic
gcf=figure();
colormap(gcf,'gray')
for  j = 1 : length(t)-1

    %tic
    b=M2*Psy;
    Psy = mldivide(M,b);
    Psy_mat=vec2mat(Psy,length(y))';
    
    surf(x,y,abs(Psy_mat).^2,'edgecolor','none');
    hold on
    surf(x,y,max(max(abs(Psy_mat).^2))*v_mat/7000,'edgecolor','none')
    hold off
    view(0,90);
    daspect([1 1 1])
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    
    graph_resultat(:,j)=abs(Psy_mat(:,floor(2*length(x)/3),1)).^2;

    F=getframe(gcf);
    writeVideo(vid1,F);
    
    %toc
    
end
%toc
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
time(indice)=toc
indice=indice+1;
end

vec=[DX DX*2 DX*4 DX*8 DX*16 DX*32 DX*64 DX*128 DX*256];

figure(1)
plot(vec,time);
title('Temps de calcul total en fonction de la discretisation dans le temps');
ylabel('Temps de Calcul[s]');
xlabel('Discretisation dx [m]');

P_x=polyfit(log10(vec(1:length(time)-4)),log10(time(1:length(time)-4)),1)

figure(2)
loglog(vec(1:length(time)-4),time(1:length(time)-4),'b','linewidth',1.5);
hold on
loglog(vec(1:length(time)-4),10.^(P_x(1)*log10(vec(1:length(time)-4)) +P_x(2)),'r','LineWidth',1.5)
title('Temps de calcul total en fonction de la discretisation dans le temps');
ylabel('Temps de Calcul (s)');
xlabel('Discretisation dx (m)');
legend('Temps de calcul expériemental','Régression linéaire (y=-0.87434x-14.2722)')
