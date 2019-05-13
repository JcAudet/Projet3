%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project          : Diffraction électronique
% Function name    : Effet_Tunnel_Graphiques.m
%
% Author           : Jean-Christophe G.-A., Thomas R., Nouffel S.
% Date             : 5/03/2019
%
% Description      : Graphique théorique vs experimental pour l'effet
% tunnel 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paramètres de l'expérience modifiables

a=100;                                   % Nombre d'iterations (on roule 'a' fois l'effet tunnel)
V_tunnel=linspace(10^-18,0.05*10^(-14),a);   % 'a' différentes valeurs de potentiel de la barrière
% V_tunnel=0.5*10^(-16);
L_mur=1.5e-11;                          % Largeur mur de potentiel
E=10^(-16);                             % Energie du paquet d'onde

%% Discretisation

global ly hbar m kp sig x0 y0

dx=5e-12;
dy=5e-12;
dt=5e-21;
xf=1e-9;
yf=2e-10;
tf=5e-17;
 
x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/3;
y0=y(end)/2;

%% Constantes et paramètres du paquet d'onde
hbar=1.0545718*10^-34;
m=9.10938*10^-31;
sig=8e-11;                                  % Équart type du paquet d'onde         
lambda=sqrt((6.62607*10^-34)^2/(2*m*E));    % Longeur d'onde, dépendante de l'énergie
kp=2*(pi/lambda)*[1,0];                     % Vecteur d'onde

% Initialisation de nos vecteurs de mesure
tunnel_transmise=zeros(1,a);
tunnel_reflechie=zeros(1,a);

for w=1:a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vid 1

[psy,norm] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

Psy=psy(:);
Psy_mat=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul des facteurs

% Definition potentiel
pot=V_tunnel(w);
D_fente=0;
L_fente=0;

v_mat=barr(x,y,pot,x(end)/2,L_mur,y(end)/2,D_fente,L_fente,'Carre');

V=v_mat(:);
[M,M2]=MM2(V,dx,dy,dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul
tic
for  j = 1 : length(t)-1

    tic
    b=M2*Psy;
    Psy = mldivide(M,b);
    Psy_mat=vec2mat(Psy,length(y))';
    toc   
end

% Calcul norme transimse et réfléchie
norm_transmise=trapeze_2D(abs(Psy_mat(:,round(length(x)/2+1):1:length(x))).^2,x(round(length(x)/2)),x(end),y(1),y(end),round(length(x)/2)-2,length(y)-1);
norm_reflechie=trapeze_2D(abs(Psy_mat(:,round(1:length(x)/2))).^2,x(1),x(round(length(x)/2)),y(1),y(end),round(length(x)/2)-2,length(y)-1);
norm=norm_transmise+norm_reflechie;

tunnel_transmise(w)=norm_transmise;
tunnel_reflechie(w)=norm_reflechie;
end
toc

%% Graphiques théorique et experimentaux
T=zeros(1,a);
R=zeros(1,a);

for i=1:1:a
V=V_tunnel(i);                    %Potentiel de la barrière 
lambda=L_mur*sqrt(2*m*V/hbar^2);
epsilon=E/V;
    if epsilon > 1
     T1=(1+(1/(4*epsilon*(epsilon-1)))*(sin(lambda*sqrt(epsilon-1)))^2)^(-1);
     R1=(1+(4*epsilon*(epsilon-1))/(sin(lambda*sqrt(epsilon-1)))^2)^(-1);
    else 
     T1=(1+1/(4*epsilon*(1-epsilon))*(sinh(lambda*sqrt(1-epsilon)))^2)^(-1);
     R1=(T1/(4*epsilon*(1-epsilon))*(sinh(lambda*sqrt(1-epsilon)))^2);
    end
T(i)=T1;
R(i)=R1;
end

figure
hold on
plot(V_tunnel,T,'blue')
plot(V_tunnel,R,'red')
plot(V_tunnel,tunnel_transmise,'--b')
plot(V_tunnel,tunnel_reflechie,'--r')

legend('Coeff de transmission théorique','Coeff de réflexion théoriques','Coeff de transmission experimentaux','Coeff de réflexion experimentaux')
xlabel('Potentiel de la barrière (eV)')
ylabel('Coefficients')
title('Coefficients de réflexion et transmission pour E=50KeV')

hold off



