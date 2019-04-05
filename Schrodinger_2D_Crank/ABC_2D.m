clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

% hbar=6.62607*10^-34;
% m=9.10938*10^-31;
hbar=1;
m=5*10^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global dx dy dt Imax Jmax

dx=10^-2;
dy=10^-2;
dt=1.75*10^-4;

Imax=110;                % Indices d'espaces max
Jmax=110;
Kmax=100;               % Indice de temps max

x=0:dx:(Imax-1)*dx;
y=0:dy:(Jmax-1)*dy;
t=0:dt:(Kmax-1)*dt;

x_0=x(end)/3;y_0=y(end)/2;    % Centre en (x,y)=(5,5)
sig_x=x(end)/15;
sig_y=y(end)/15;% Proprietse initiales
lamda=1;    % Accelerer du repos avec 50KeV
% kp=2*(pi/lamda)*[1,1];
kp=117*[-1,0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONDITIONS INITIALES
% La condition initiale est une fonciton Psy(x,y). Pour representer une 
% fonction de 2 variables dans un vecteur, il faut adopter une convention.

norm=zeros(1,length(t));
[psy,norm(1)] = wp_ini_2D(x,y,sig_x,sig_y,kp,x_0,y_0);

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

%% ABC's 2D
kx=kp(1);
ky=kp(2);
h_barre=1;

% ax=h_barre*kx/m;
% ay=h_barre*ky/m;
% 
% h1=m*(3*ax);
% h2=m^2*(ax^3)*(3/ax);
% h3=m^3*ax^3;
%
% Conditions frontières
% C1= -(-1i*h_barre^3/dx^3 - h_barre^2*h1/(2*dx));
% C2= 2*1i*h_barre^3/dx^3 + 2*m*h_barre^2/(dx*dt) + h_barre^2*h1/dx^2 - 1i*m*h_barre*h1/dt + 1i*h_barre/(2*dx) + h3/4;
% C3= -(-1i*h_barre^3 - h_barre^2*h1/(2*dx^2));
% C4= -(1i*h_barre^3/dx^3 - h_barre^2*h1/(2*dx^2));
% C5= -2*1i*h_barre^3/dx^3 - 2*m*h_barre^2/(dx*dt) + h_barre^2*h1/dx^2 - 1i*m*h_barre*h1/dt + 1i*h_barre*h2/(2*dx) + h3/4;
% C6= -(1i*h_barre^3/dx^3 - h_barre^2*h1/(2*dx^2));
% C7=-(-2*m*h_barre^2/(dx*dt) + 1i*m*h_barre*h1/dt + 1i*h_barre*h2/(2*dx) + h3/4);
% C8= -(2*m*h_barre^2/(dx*dt) + 1i*m*h_barre*h1/dt + 1i*h_barre*h2/(2*dx) + h3/4);

A1=1i*h_barre^3/(2*m)*(1/dx^3);
A2=1i*h_barre^3/(2*m)*(-2/dx^3);
A3=1i*h_barre^3/(2*m)*(1/dx^3);
A4=1i*h_barre^3/(2*m)*(-1/dx^3);
A5=1i*h_barre^3/(2*m)*(2/dx^3);
A6=1i*h_barre^3/(2*m)*(-1/dx^3);

B1=-h_barre^2*(1/(dx*dt));
B2=-h_barre^2*(-1/(dx*dt));
B3=-h_barre^2*(-1/(dx*dt));
B4=-h_barre^2*(1/(dx*dt));

D1=3*h_barre^3*kx/(2*m)*(1/(2*dx^2));
D2=3*h_barre^3*kx/(2*m)*(-2/(2*dx^2));
D3=3*h_barre^3*kx/(2*m)*(1/(2*dx^2));
D4=3*h_barre^3*kx/(2*m)*(1/(2*dx^2));
D5=3*h_barre^3*kx/(2*m)*(-2/(2*dx^2));
D6=3*h_barre^3*kx/(2*m)*(1/(2*dx^2));

E1=3*1i*h_barre^2*kx*(1/(2*dt));
E2=3*1i*h_barre^2*kx*(-1/(2*dt));
E3=3*1i*h_barre^2*kx*(1/(2*dt));
E4=3*1i*h_barre^2*kx*(-1/(2*dt));

F1=3*1i*h_barre^3*kx^3/(2*m)*(1/(2*dx));
F2=3*1i*h_barre^3*kx^3/(2*m)*(-1/(2*dx));
F3=3*1i*h_barre^3*kx^3/(2*m)*(1/(2*dx));
F4=3*1i*h_barre^3*kx^3/(2*m)*(-1/(2*dx));

G1=h_barre^3*kx^3/(2*m)*(1/4);
G2=h_barre^3*kx^3/(2*m)*(1/4);
G3=h_barre^3*kx^3/(2*m)*(1/4);
G4=h_barre^3*kx^3/(2*m)*(1/4);

C1=(A1+D1);
C2=(A2+B1+D2+E1+F1+G1);
C3=(A3+D3);
C4=(A4+D4);
C5=(A5+B3+D5+E3+F2+G3);
C6=(A6+D6);
C7=(B2+E2+F3+G2);
C8=(B4+E4+F4+G4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul

Psy_mat=zeros(Jmax,Imax,Kmax);
Psy_mat(:,:,1)=psy;

tic
figure()
for k = 2 : Kmax
    b=M2*Psy(:,k-1)+v2-v;
    Psy(:,k) = mldivide(M,b);
    
    for j=1:Imax
            if j==1
               Psy(Imax,k)=(C1*Psy(Imax*2,k-1)+C7*Psy(Imax,k-1))/C2;
               Psy(Imax-1,k)=(C4*Psy(Imax*2-1,k-1)+C8*Psy(Imax-1,k-1))/C5; 
            elseif j==Imax
               Psy(Imax^2,k)=(C3*Psy(Imax*(j-1),k-1)+C7*Psy(Imax*j,k-1))/C2;
               Psy(Imax^2-1,k)=(C6*Psy((Imax*(j-1))-1,k-1)+C8*Psy(Imax*(j)-1,k-1))/C5;   
            else
            Psy(Imax*j,k)=(C1*Psy(Imax*(j+1),k-1)+C3*Psy(Imax*(j-1),k-1)+C7*Psy(Imax*j,k-1))/C2;
            Psy((Imax*j)-1,k)=(C4*Psy(Imax*(j+1)-1,k-1)+C6*Psy((Imax*(j-1))-1,k-1)+C8*Psy(Imax*(j)-1,k-1))/C5;
            end
    end
    
    Psy_mat(:,:,k)=vec2mat(Psy(:,k),Imax);
    norm(k)=trapeze_2D(abs(Psy_mat(:,:,k)).^2,x(1),x(end),y(1),y(end),Imax-1,Jmax-1);
    surf(x,y,abs(Psy_mat(:,:,k)).^2);
    hold on
    surf(x,y,v_mat)
    hold off
    title(sprintf('Norme: %.10f',norm(k)))
    view(0,90);
    pause(0.02);
end
toc