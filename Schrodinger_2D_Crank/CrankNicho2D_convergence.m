clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

hbar=6.62607*10^-34;
m=9.10938*10^-31;

sig_x=1e-13;
sig_y=1e-13;% Proprietse initiales
lamda=5.4847*10^-12;    % Accelerer du repos avec 50KeV
% kp=2*(pi/lamda)*[1,1];
kp=1e16*[1,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global dx dy dt Imax Jmax

dx=10^-14;
dy=dx;
dt=10^-25;

Imax=50;                % Indices d'espaces max
Jmax=Imax;
Kmax=50;               % Indice de temps max

xyi=0;
xyf=(Imax-1)*dx;
x1=xyi:dx:xyf;
y1=x1;
t=0:dt:(Kmax-1)*dt;


Dxy=[dx dx/2];
norm_x=zeros(length(Dxy),length(t));
Psy_mem=zeros(length(Dxy),Kmax,Jmax,Imax);

for w=1:length(Dxy)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initialisation
    
    x=0:Dxy(w):xyf;
    y=x;
   
    [psy,norm_x(w,1)] = wp_ini_2D(x,y,sig_x,sig_y,kp,x(end)/4,y(end)/2);

    Psy=zeros(length(t),length(y)*length(x));
    Psy(1,:)=psy(:);
    Psy_mem(w,1,:,:)=psy(1:2^(w-1):length(psy),1:2^(w-1):length(psy));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat = zeros(length(y),length(x));
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

    M=diag(b);
    v=zeros(length(b),1);

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul
    
    
    for k = 1 : Kmax-1

        b=M2*Psy(k,:)'+v2-v;
        Psy(k+1,:) = mldivide(M,b);
        psy=vec2mat(Psy(k+1,:),length(x));
        Psy_mem(w,k+1,:,:)=psy(1:2^(w-1):length(psy),1:2^(w-1):length(psy));
    end

    w
end


err_x=zeros(length(Dxy)-1,length(t),length(0:Dxy(1):5));
err_x_int = zeros(length(Dx)-1,length(t));

% Calcul erreur
for i=1:length(Dx)-1
    
    err_x(i,:,:)=abs(Psy_mem_x(i,:,:)-Psy_mem_x(i+1,:,:));
    
    for j = 1:length(t)
        
        err_x_int(i,j)=trapeze(err_x(i,j,:),0,5,length(err_x(i,j,:))-1);
        
    end
end


figure()
hold on
for i=1:length(Dx)-1
    plot(t,err_x_int(i,:))
end
legend(sprintf('dx = %f',Dx(1)),sprintf('dx = %f',Dx(2)),sprintf('dx = %f',Dx(3)),sprintf('dx = %f',Dx(4)),sprintf('dx = %f',Dx(5)))
title('Erreur en fonction du temps et de la discretisation de l''espace')



