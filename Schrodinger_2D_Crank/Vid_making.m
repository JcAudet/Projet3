clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project          : Diffraction électronique
%% Function name    : Fichier principal de résolution
%
%% Author           : Jean-Christophe G.-A., Thomas R., Nouffel S.
%% Date             : 5/1/2019
%
%% Description      : Propagation d'un paquet d'onde selon les parmètres
%%                  définies dans la section "Définitions"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Définitions

profile on

global ly hbar m

% Constantes physiques
hbar=1.0545718*10^-34;
m=9.10938*10^-31;

% Paramètre du paquet d'onde
sig=8e-11;           
lamda=5e-11;    
kp=2*(pi/lamda)*[1,0];

% Discretisation spatiale et temporel
dx=5e-12;
dy=5e-12;
dt=1e-21;
xf=3e-11;
yf=3e-11;
tf=4e-17;

% Creation du domaine
x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/3;
y0=y(end)/2;

% Definition potentiel
pot=10000;
L_mur=1.1e-11;
D_fente=20e-11;
L_fente=3e-11;

% Initialisation Écrant
graph_resultat=zeros(length(y),length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation du paquet d'onde et mémoire

[psy,norm] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

Psy=psy(:);
Psy_mat=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation des matrices du système

% Création du potentiel
v_mat=barr(x,y,pot,x(end)/2,L_mur,y(end)/2,D_fente,L_fente,'Carre');
% v_mat=barr_simple(x,y,pot,x(end)/2,L_mur,y(end)/2,L_fente,'Carre');
V=v_mat(:);

% Creation des matrices
[M,M2]=MM2(V,dx,dy,dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation Video

vid1 = VideoWriter('Simple fente','MPEG-4');
vid1.FrameRate = 150;
vid1.Quality = 100;
open(vid1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propagation

gcf=figure();
colormap(gcf,'hot')

for  j = 1 : length(t)-1

    % Resolution du système linéraire
    b=M2*Psy;
    Psy = mldivide(M,b);
    Psy_mat=vec2mat(Psy,length(y))';
    
    % Affichage
    surf(x,y,abs(Psy_mat).^2,'edgecolor','none');
    hold on
    surf(x,y,max(max(abs(Psy_mat).^2)).*v_mat./7000,'edgecolor','none')
    hold off
    view(0,90);
    daspect([1 1 1])
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    
    % Prise de donnée de l'écrant
    graph_resultat(:,j)=abs(Psy_mat(:,floor(2*length(x)/3),1)).^2;

    % Création vidéo
    F=getframe(gcf);
    writeVideo(vid1,F);

    
end
close(vid1)

profile viewer


clear M M2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul du modère de Fraunhofer pour la configuration

L=y(end);
H=(2/3-0.5)*x(end);
thetamax=atan(L/(2*H));
theta=-thetamax:0.01:thetamax;

alpha=(pi*D_fente/lamda).*sin(theta);
beta=(pi*L_fente/lamda).*sin(theta);
I=(cos(alpha).*sin(beta)./beta).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Traitement des résultats de l'écran

maxb=max(graph_resultat);
i=find(maxb==max(maxb));
F=graph_resultat(:,i)/max(graph_resultat(:,i));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Affichage de la comparaison

figure
plot(-y(end)/2:dy:y(end)/2,F,'r')
hold on
plot(H*tan(theta),I,'b');
title(sprintf('Résulats expérimentaux vs théoriques (e=%e, d=%e, a=%e)',L_mur,D_fente,L_fente));
xlabel('Position sur l''écran (m)')
ylabel('Intensité normalisée')
legend('Patron de diffraction','Fraunhofer')

