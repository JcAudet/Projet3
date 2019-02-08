clc;clear all;

%En utilisant la methode Crank-Nicholson presentée dans les diapos de
%PHS3903, voir les diapos pour la resolutions matricielle

N=100;%nombre de point sur la grille en x
%le vecteur d'onde contient N-2 point. -2 parce que les deux extremite
%gauche et droite sont fixee, c'est les conditions frontiere.

D = complex(0,1)*(6*10^-34)/(4*pi*9.109*10^-31);
%D dans le systeme dunite du document schrod dynamique : on prend hbar=1 et
%m=1
D=complex(0,1)/2;



deltax=0.01;
deltat=10^-5;


debut=-5;
fin=5;;
x=debut+deltax:deltax:fin-deltax;
%conditions frontieres
Psy_debut=0;
Psy_fin=0;

N=length(x)+2;

PsyInitial=Psy_0(x);

%%ECRIRE LES DEUX MATRICES ATTENTION AU DIMENSIO
%Ecriture des parametre pour les matrices
gamma=1-D*deltat/(deltax^2);
alpha=1+D*deltat/(deltax^2);
beta=-0.5*D*deltat/(deltax^2);
%Construction de C et c
C=full(gallery('tridiag',length(x),-beta,gamma,-beta));
c(1)=-2*beta*Psy_debut;
c(length(x))=-2*beta*Psy_fin;

%Le systeme a resoudre est de la forme A*(Psy_t+deltat)=b, ou b=C*(Psy_t)+c

A=full(gallery('tridiag',length(x),-beta,alpha,-beta));

Psy=transpose(PsyInitial);
modulecarre=transpose(Psy)*conj(Psy)*deltax;
Psy=Psy/modulecarre^0.5;
cell{1}=Psy;

k=1;
m=1000;
for i=1:m
   %c(1)=-2*beta*Psy(1);
   %c(length(x))=-2*beta*Psy(end);
    b=C*Psy + transpose(c);
    Psy_t_deltat = A\b;
    modulecarre=transpose(Psy_t_deltat)*conj(Psy_t_deltat)*deltax;
    Psy_t_deltat=Psy_t_deltat/modulecarre^0.5;
    cell{1+i}=Psy_t_deltat;
    Psy=Psy_t_deltat;
    plot(x,Psy_t_deltat);
    axis([-4 4 -2 2])
    pause(0.01)
    
    %loading
    loading=100*k/(m);
    disp(loading);
    k=k+1;
end
% 
% n=1;
% 
%  while n<length(cell)
%  
%  plot(x,(cell{n}));
%  axis([-4 4 -2 2])
%  pause(0.05);  
%  n=n+1;
%  end
% plot(x,cell{end})



