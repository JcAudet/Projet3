clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

hbar=6.62607*10^-34;
m=9.10938*10^-31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global lx ly

dx=1e-11;     % On veut que les dx et dt soit au environ 10x la largeur du packet d'onde
dy=1e-11;
dt=1e-21;

x=0:dx:1e-9;
y=0:dy:1e-9;
t=0:dt:1e-16;

lx=length(x);
ly=length(y);

x_0=x(end)/4;
y_0=y(end)/2;    % Centre en (x,y)=(5,5)
sig_x=1e-10;                 % Largeur initiale du packet d'onde
sig_y=1e-10;
lamda=5e-12;    % Accelerer du repos avec 50KeV
kp=2*(pi/lamda)*[1 0];


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

v_mat=zeros(length(y),length(x));
% v_mat=( barr(x,y,1000,x(floor(length(x)/2)),4e-11,y(floor(length(y)/2)),15e-11,10e-11,'Carre') )';
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
for j = 1 : length(y)
    for i = 1 : length(x)

        ind=indexeur(i,j);
        
        if i > 1
            M(ind,indexeur(i-1,j)) = c;
            M2(ind,indexeur(i-1,j)) = g;
        end
        
        if i < length(x)
            M(ind,indexeur(i+1,j)) = c;
            M2(ind,indexeur(i+1,j)) = g;
        end
        
        if j < length(y)
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
for k = 1 : length(t)-1
    b=M2*Psy(:,k);
    Psy(:,k+1) = mldivide(M,b);
    
    Psy_mat(:,:,k+1)=vec2mat(Psy(:,k+1),length(x))';
    norm(k+1)=trapeze_2D(abs(Psy_mat(:,:,k+1)).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    
    subplot(1,5,[1 2 3])    
    surf(x,y,abs(Psy_mat(:,:,k)).^2,'edgecolor','none');
    hold off
    title(sprintf('Temps = %e  Norme: %.10f',t(k),norm(k)))
    view(0,90);
    
    subplot(1,5,[4 5])
    plot(abs(Psy_mat(:,floor(2*length(x)/3),k)).^2,y,'b')
    xlim([0 5e19])
    
    F=getframe(gcf);
    writeVideo(vid,F);
    
end
toc

close(vid)
