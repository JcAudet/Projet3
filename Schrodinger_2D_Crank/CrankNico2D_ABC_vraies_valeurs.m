clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

hbar=6.62607*10^-34*(1/(2*pi));
m=9.10938*10^-31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global dx dy dt Imax Jmax

dx=10^-11;     % On veut que les dx et dt soit au environ 10x la largeur du packet d'onde
dy=10^-11;
dt=10^-18;

Imax=150;                % Indices d'espaces max
Jmax=150;
Kmax=250;               % Indice de temps max

x=0:dx:(Imax-1)*dx;
y=0:dy:(Jmax-1)*dy;
t=0:dt:(Kmax-1)*dt;

x_0=x(end)/4;y_0=y(end)/2;    % Centre en (x,y)=(5,5)

sig_x=1e-10;                 % Largeur initiale du packet d'onde
sig_y=1e-10;

lamda=5.4847*10^-12;    % Accelerer du repos avec 50KeV
kp=2*(pi/lamda)*[1,0];
% kp=1e17*[1,0];


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

v_mat=( barr(x,y,1000,x(length(x)/2),4e-11,y(length(y)/2),15e-11,10e-11,'Carre') )';
% v_mat=sparse( ( barr_simple(x,y,1000,x(length(x)/2),1e-14,y(length(y)/2),8e-14,'Carre') )' );
V=v_mat(:);

b = 1 + 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) + 1i*dt*V/(2*hbar);
f = 1 - 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) - 1i*dt*V/(2*hbar);

c=-1i*hbar*dt*(1/dx^2)/(4*m);
d=-1i*hbar*dt*(1/dy^2)/(4*m);

g=1i*hbar*dt*(1/dx^2)/(4*m);
k=1i*hbar*dt*(1/dy^2)/(4*m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creation de M et M2, v et v2

%M=sparse(ones(length(b_vect),1),ones(length(b_vect),1),b_vect);
M=diag(b);

%M2=sparse(ones(length(b_vect),1),ones(length(b_vect),1),f_vect);
M2=diag(f);

tic
for j = 1 : Jmax
    for i = 1 : Imax 

        ind=indexeur(i,j);
        
        if i > 1
            M(ind,indexeur(i-1,j)) = c;
            M2(ind,indexeur(i-1,j)) = g;
        end
        
        if i < Imax 
            M(ind,indexeur(i+1,j)) = c;
            M2(ind,indexeur(i+1,j)) = g;
        end
        
        if j < Jmax 
            M(ind,indexeur(i,j+1)) = d;
            M2(ind,indexeur(i,j+1)) = k;
        end
        
        if j > 1
            M(ind,indexeur(i,j-1)) = d;
            M2(ind,indexeur(i,j-1)) = k;
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
h_barre=hbar;

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

% Creation video
vid = VideoWriter('Schrod_Crank_Diffra','MPEG-4');
vid.FrameRate = 60;
open(vid)

tic
gcf=figure();
for k = 2 : Kmax
    b=M2*Psy(:,k-1);
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
  
    
    subplot(1,5,[1 2 3])
    Psy_mat(:,:,k)=vec2mat(Psy(:,k),Imax);
    norm(k)=trapeze_2D(abs(Psy_mat(:,:,k)).^2,x(1),x(end),y(1),y(end),Imax-1,Jmax-1);
    surf(x,y,abs(Psy_mat(:,:,k)).^2,'edgecolor','none');
    hold on
    surf(x,y,2e19*v_mat'/1000)
    hold off
    title(sprintf('Temps = %e  Norme: %.10f',t(k),norm(k)))
    view(-60,20);
    
    subplot(1,5,[4 5])
    plot(abs(Psy_mat(:,floor(2*length(x)/3),k)).^2,y)
    xlim([-2e19 2e19])
    
    F=getframe(gcf);
    writeVideo(vid,F);
    
end
toc

close(vid)