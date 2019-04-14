clear all
clc
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global ly hbar m kp sig x0 y0

T=zeros(2,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DX

dx=0.94e-11;
Dxy=[dx dx*2 dx*4 dx*8 dx*16];

for w=1:length(Dxy)
    
    tic
    
    hbar=6.62607*10^-34;
    m=9.10938*10^-31;   

    sig=8e-11;           
    lamda=1e-10;    
    kp=2*(pi/lamda)*[1,0];

    dx=Dxy(w);
    dy=Dxy(w);
    dt=1e-20;
    xf=2e-9;
    yf=1.4e-9;
    tf=2e-17;

    x=0:dx:xf;
    y=0:dx:yf; ly=length(y);
    t=0:dt:tf;
    x0=x(end)/3;
    y0=y(end)/2;
    
    norm=zeros(1,length(t));
    [psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

    Psy=[];
    Psy(:,1)=psy(:);
    Psy_mat=psy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat=zeros(length(y),length(x));
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

    M=sparse(M);
    M2=sparse(M2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

    for  j = 1 : length(t)-1

        b=M2*Psy;
        Psy = mldivide(M,b);

        Psy_mat=vec2mat(Psy,length(y))';

        norm(j+1)=trapeze_2D(abs(Psy_mat).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    end

    T(1,w)=toc()
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DT
clear all

dt=1e-20;
Dt=[dt dt*2 dt*4 dt*8 dt*16];

for w=1:length(Dx)
    
    tic
    
    hbar=6.62607*10^-34;
    m=9.10938*10^-31;   

    sig=8e-11;           
    lamda=1e-10;    
    kp=2*(pi/lamda)*[1,0];

    dx=0.94e-11;
    dy=0.94e-11;
    dt=Dt(w);
    xf=2e-9;
    yf=1.4e-9;
    tf=2e-17;

    x=0:dx:xf;
    y=0:dy:yf; ly=length(y);
    t=0:dt:tf;
    x0=x(end)/3;
    y0=y(end)/2;
    
    norm=zeros(1,length(t));
    [psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

    Psy=[];
    Psy(:,1)=psy(:);
    Psy_mat=psy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat=zeros(length(y),length(x));
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

    M=sparse(M);    
    M2=sparse(M2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

    for  j = 1 : length(t)-1

        b=M2*Psy;
        Psy = mldivide(M,b);

        Psy_mat=vec2mat(Psy,length(y))';

        norm(j+1)=trapeze_2D(abs(Psy_mat).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    end

    T(2,w)=toc()
end

