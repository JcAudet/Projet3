clc;
clear all;
close all;

global dx dt

%En utilisant la methode Crank-Nicholson presentée dans les diapos de
%PHS3903, voir les diapos pour la resolutions matricielle

N=1000;%nombre de point sur la grille en x
%le vecteur d'onde contient N-2 point. -2 parce que les deux extremite
%gauche et droite sont fixee, c'est les conditions frontiere.

D = 1i; %complex(0,1)*(6*10^-34)/(4*pi*9.109*10^-31);
%D dans le systeme dunite du document schrod dynamique : on prend hbar=1 et
%m=1

%% Discretisation
x=linspace(-5,5,N);
dx=x(2)-x(1);
t=linspace(0,1,N);
dt=10^-5;


%% Initialisation et normalisation
Psy_i=Psy_0(x);
Norm=sum(abs(Psy_i))*dx;
Psy_i=Psy_i/Norm;
Psy=zeros(N,length(x));
Psy(1,:)=Psy_i;
% PsyInitial(1)=0;
% PsyInitial(end)=0;

%%ECRIRE LES DEUX MATRICES, ATTENTION AU DIMENSIO
%Ecriture des parametre pour les matrices
gamma=1-D*dt/(dx^2);
alpha=1+D*dt/(dx^2);
beta=-0.5*D*dt/(dx^2);
%% Construction de C et c
C=full(gallery('tridiag',length(x),-beta,gamma,-beta));

c=zeros(1,length(x));
c(1)=-2*beta*Psy(1);
c(end)=-2*beta*Psy(end);

%Le systeme a resoudre est de la forme A*(Psy_t+deltat)=b, ou b=C*(Psy_t)+c

A=full(gallery('tridiag',length(x),beta,alpha,beta));

figure
for i=1:N
    b=C*transpose(Psy(i,:));% + transpose(c);
    Psy(i+1,:) = A\b;
    plot(x,Psy(i,:));
    axis([-4 4 -2 2])
    pause(0.01)
end
% 
% figure
% for k = 1 : length(t)
% %     norm=sum(abs(Psy(k,:)))*dx;
% %     plot(x, real(Psy(k,:)),'b');
% %     hold on
% %     plot(x,imag(Psy(k,:)),'r');
% %     hold on
%     plot(x,abs(Psy(k,:)))
% %     hold on
%     ylim([0 1.5]);
% %     plot(x,pot)
%     hold off
% %     title( sprintf('t = %.5f , Norme= %.6f', t(k), norm));
% %     pause(0.0000001);
% end



