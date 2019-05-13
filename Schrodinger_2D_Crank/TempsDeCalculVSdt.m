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
time=zeros(1,9);
indice=1;
for dt=[1e-18 2e-18 4e-18 8e-18 16e-18 32e-18 64e-18 128e-18 256e-18]
    tic
%dt=1e-17;
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

vec=[1e-18 2e-18 4e-18 8e-18 16e-18 32e-18 64e-18 128e-18];

figure(1)
plot([1e-18 2e-18 4e-18 8e-18 16e-18 32e-18 64e-18 128e-18],time(1:length(t)-1));
title('Temps de calcul total en fonction de la discretisation dans le temps');
ylabel('Temps de Calcul[s]');
xlabel('Discretisation dt [s]');

P_t=polyfit(log10([1e-18 2e-18 4e-18 8e-18 16e-18 32e-18 64e-18 128e-18]),log10(time(1:length(time)-1)),1)


figure(2)
loglog(vec,time(1:length(time)-1),'b','linewidth',1.5);
hold on
loglog(vec,10.^(P_t(1)*log10(vec) +P_t(2)),'r','LineWidth',1.5)
title('Temps de calcul total en fonction de la discretisation dans le temps');
ylabel('Temps de Calcul (s)');
xlabel('Discretisation dt (s)');
legend('Temps de calcul expériemental','Régression linéaire (y=-0.87434x-14.2722)')






