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
Kmax=250;               % Indice de temps max

x=0:dx:(Imax-1)*dx;
y=0:dy:(Jmax-1)*dy;
t=0:dt:(Kmax-1)*dt;

x_0=x(end)/4;y_0=y(end)/2;    % Centre en (x,y)=(5,5)
sig_x=1e-13;
sig_y=1e-13;% Proprietse initiales
lamda=5.4847*10^-12;    % Accelerer du repos avec 50KeV
% kp=2*(pi/lamda)*[1,1];
kp=1e16*[1,0];


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

v_mat=( barr(x,y,1000,x(length(x)/2),2e-14,y(length(y)/2),12e-14,10e-14,'Carre') )';
% v_mat=sparse( ( barr_simple(x,y,1000,x(length(x)/2),1e-14,y(length(y)/2),8e-14,'Carre') )' );
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

for j = 1 : length(y)
    for i = 1 : length(x)

        ind=indexeur(i,j);
        
        if i > 1
            M(ind,indexeur(i-1,j)) = c;
            M2(ind,indexeur(i-1,j)) = h;
        else
            v(ind)=v(ind)+c*Psy_haut(j);
            v2(ind)=v2(ind)+g*Psy_haut(j);
        end
        
        if i < length(x)
            M(ind,indexeur(i+1,j)) = a;
            M2(ind,indexeur(i+1,j)) = g;
        else
            v(ind)=v(ind)+a*Psy_bas(j);
            v2(ind)=v2(ind)+h*Psy_bas(j);
        end
        
        if j < length(y)
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

Psy_mat=zeros(length(y),length(x),length(t));
Psy_mat(:,:,1)=psy;

% Creation video
vid = VideoWriter('Schrod_Crank_Diffra','MPEG-4');
vid.FrameRate = 60;
open(vid)

tic
gcf=figure();
for k = 2 : length(t)
    b=M2*Psy(:,k-1)+v2-v;
    Psy(:,k) = mldivide(M,b);
    
    subplot(1,5,[1 2 3])
    Psy_mat(:,:,k)=vec2mat(Psy(:,k),length(x));
    norm(k)=trapeze_2D(abs(Psy_mat(:,:,k)).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    surf(x,y,abs(Psy_mat(:,:,k)).^2,'edgecolor','none');
    hold off
    title(sprintf('Temps = %e  Norme: %.10f',t(k),norm(k)))
    view(0,90);
    
    subplot(1,5,[4 5])
    plot(abs(Psy_mat(:,floor(2*length(x)/3),k)).^2,y)
    xlim([-2e25 2e25])
    
    F=getframe(gcf);
    writeVideo(vid,F);
    
end
toc

close(vid)
