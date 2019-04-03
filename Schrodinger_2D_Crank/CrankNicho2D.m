clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

hbar=6.62607*10^-34;
m=9.10938*10^-31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global dx dy dt Imax Jmax

dx=10^-14;
dy=10^-14;
dt=10^-25;

Imax=150;                % Indices d'espaces max
Jmax=150;
Kmax=1000;               % Indice de temps max

x=0:dx:(Imax-1)*dx;
y=0:dy:(Jmax-1)*dy;

x_0=x(end)/3;y_0=y(end)/2;    % Centre en (x,y)=(5,5)
sig_x=1e-13;
sig_y=1e-13;% Proprietse initiales
lamda=5.4847*10^-12;    % Accelerer du repos avec 50KeV
% kp=2*(pi/lamda)*[1,1];
kp=1e15*[1,0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONDITIONS INITIALES
% La condition initiale est une fonciton Psy(x,y). Pour representer une 
% fonction de 2 variables dans un vecteur, il faut adopter une convention.

[psy,norm] = wp_ini_2D(x,y,sig_x,sig_y,kp,x_0,y_0);

Psy=[];
Psy(:,1)=psy(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul des facteurs

v_mat=( barr(x,y,1000,7e-13,0.2e-13,y(end)/2,1e-13,0.4e-13,'Carre') )';
V=v_mat(:);

b = 1 + 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) + 1i*dt*V/(2*hbar);
f = 1 - 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) - 1i*dt*V/(2*hbar);

c=-1i*hbar*dt*(1/dx^2)/(4*m);
a=c;
d=-1i*hbar*dt*(1/dy^2)/(4*m);
e=d;

g=1i*hbar*dt*(1/dx^2)/(4*m);
h=g;
k=1i*hbar*dt*(1/dy^2)/(4*m);
p=k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creation de M et M2, v et v2

%M=sparse(ones(length(b_vect),1),ones(length(b_vect),1),b_vect);
M=diag(b);
v=zeros(length(b),1);

%M2=sparse(ones(length(b_vect),1),ones(length(b_vect),1),f_vect);
M2=diag(f);
v2=zeros(length(f),1);

tic

for j = 1 : Jmax
    for i = 1 : Imax 

        ind=indexeur(i,j);
        
        if i > 1
            M(ind,indexeur(i-1,j)) = c;
            M2(ind,indexeur(i-1,j)) = h;
        else
            v(ind)=v(ind)+c*Psy_haut(j);
            v2(ind)=v2(ind)+g*Psy_haut(j);
        end
        
        if i < Imax 
            M(ind,indexeur(i+1,j)) = a;
            M2(ind,indexeur(i+1,j)) = g;
        else
            v(ind)=v(ind)+a*Psy_bas(j);
            v2(ind)=v2(ind)+h*Psy_bas(j);
        end
        
        if j < Jmax 
            M(ind,indexeur(i,j+1)) = d;
            M2(ind,indexeur(i,j+1)) = k;
        else
            v(ind)=v(ind)+d*Psy_droite(i);
            v2(ind)=v2(ind)+k*Psy_droite(i);
        end
        
        if j > 1
            M(ind,indexeur(i,j-1)) = e;
            M2(ind,indexeur(i,j-1)) = p;
        else
            v(ind)=v(ind)+e*Psy_gauche(i);
            v2(ind)=v2(ind)+p*Psy_gauche(i);
        end
    end
end
toc

tic
M=sparse(M);
M2=sparse(M2);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul

Psy_mat=zeros(Jmax,Imax,Kmax);
Psy_mat(:,:,1)=psy;

tic
figure()
for k = 2 : Kmax
    b=M2*Psy(:,k-1)+v2-v;
    Psy(:,k) = mldivide(M,b);
    
    Psy_mat(:,:,k)=vec2mat(Psy(:,k),Imax);
    norm=trapeze_2D(abs(Psy_mat(:,:,k)),x(1),x(end),y(1),y(end),Imax-1,Jmax-1);
    surf(x,y,abs(Psy_mat(:,:,k)));
    hold on
    surf(x,y,v_mat)
    hold off
    title(sprintf('Norme: %.10f',norm))
    view(0,90);
    pause(0.02);
end
toc


