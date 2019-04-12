clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

hbar=6.62607*10^-34;
m=9.10938*10^-31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global ly

dx=0.68e-11;     % On veut que les dx et dt soit au environ 10x la largeur du packet d'onde
dy=0.68e-11;
dt=5e-19;

x=0:dx:1e-9;
y=0:dy:1e-9;
t=0:dt:5e-17;

lx=length(x);
ly=length(y);

x_0=x(end)/4;
y_0=y(end)/2;    % Centre en (x,y)=(5,5)
sig_x=1e-10;                 % Largeur initiale du packet d'onde
sig_y=1e-10;
lamda=5e-11;    % Accelerer du repos avec 50KeV
kp=2*(pi/lamda)*[1 0];


z=zeros(length(x),length(y));
z(:,floor(2*length(y)/3))=5*10^(19);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONDITIONS INITIALES
% La condition initiale est une fonciton Psy(x,y). Pour representer une 
% fonction de 2 variables dans un vecteur, il faut adopter une convention.

norm=zeros(1,length(t));
[psy,norm(1)] = wp_ini_2D(x,y,sig_x,sig_y,kp,x_0,y_0);

Psy=[];
Psy(:,1)=psy(:);

% Psy_mat=zeros(length(y),length(x),length(t));
Psy_mat=zeros(length(y),length(x),1);
Psy_mat(:,:,1)=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul des facteurs

% v_mat=zeros(length(y),length(x));
v_mat=( barr(x,y,1000,x(floor(length(x)/2)),4e-11,y(floor(length(y)/2)),15e-11,10e-11,'Carre') );
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

        ind=indexeur(j,i);
        
        if i > 1
            M(ind,indexeur(j,i-1)) = c;
            M2(ind,indexeur(j,i-1)) = g;
        end
        
        if i < length(x)
            M(ind,indexeur(j,i+1)) = c;
            M2(ind,indexeur(j,i+1)) = g;
        end
        
        if j < length(y)
            M(ind,indexeur(j+1,i)) = d;
            M2(ind,indexeur(j+1,i)) = k;
        end
        
        if j > 1
            M(ind,indexeur(j-1,i)) = d;
            M2(ind,indexeur(j-1,i)) = k;
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

% Creation video
vid = VideoWriter('Schrod_Crank_Diffra','MPEG-4');
vid.FrameRate = 60;
%Qualite = 4k
open(vid)

tic
gcf=figure();
for k = 1 : length(t)-1
    b=M2*Psy(:,k);
    Psy(:,k+1) = mldivide(M,b);
        
    %Psy_mat(:,:,k+1)=vec2mat(Psy(:,k+1),length(y))';
    Psy_mat(:,:,1)=vec2mat(Psy(:,k+1),length(y))';
    norm(k+1)=trapeze_2D(abs(Psy_mat(:,:,1)).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
   
colormap('gray')

% map = [0 0 0
%     0.80 0.80 0.80
%     0.81 0.81 0.81
%     0.82 0.82 0.82
%     0.83 0.83 0.83
%     0.84 0.84 0.84
%     0.85 0.85 0.85
%     0.86 0.86 0.86
%     0.87 0.87 0.87
%     0.88 0.88 0.88
%     0.89 0.89 0.89
%     0.90 0.90 0.90
%     0.91 0.91 0.91
%     0.92 0.92 0.92
%     0.93 0.93 0.93
%     0.94 0.94 0.94
%     0.95 0.95 0.95
%     0.96 0.96 0.96
%     0.97 0.97 0.97
%     0.98 0.98 0.98
%     0.99 0.99 0.99
%     1 1 1];
% colormap(map)

subplot(1,5,[1 2 3])      
surf(x,y,abs(Psy_mat(:,:,1)).^2,'edgecolor','none');
% hold on
% surf(x,y,z,'edgecolor','none')
% %plot3([floor(2*length(y)/3), floor(2*length(y)/3)],[y(1) y(end)],[0 0])
% hold off
title(sprintf('Temps = %e  Norme: %.10f',t(k),norm(k)))
view(0,90);
 
subplot(1,5,[4 5])
plot(abs(Psy_mat(:,floor(2*length(x)/3),1)).^2,y,'b')
xlim([0 5e19])
    
F=getframe(gcf);
writeVideo(vid,F);   
end
toc

close(vid)
close(vid)