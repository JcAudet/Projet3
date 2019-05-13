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

dx=0.1e-11;
Dxy=[dx dx*2 dx*4 dx*8 dx*16 dx*32];

for w=1:length(Dxy)
    
    tic
    
    hbar=1.0545718*10^-34;
    m=9.10938*10^-31;

    sig=8e-11;           
    lamda=5e-11;    
    kp=2*(pi/lamda)*[1,0];

    dx=Dxy(w);
    dy=Dxy(w);
    dt=1e-20;
    xf=2e-9;
    yf=1.4e-9;
    tf=1e-19;

    x=0:dx:xf;
    y=0:dy:yf; ly=length(y);
    t=0:dt:tf;
    x0=x(end)/3;
    y0=y(end)/2;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% vid 1

    [psy,norm] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

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

    diag=ones(1,length(V)-1)*d;
    diag1=ones(1,length(V)-1)*d;
    diag(ly:ly:length(V)-1)=0;
    diag1(ly:ly:length(V)-1)=0;

    diag2=ones(1,length(V)-1)*k;
    diag3=ones(1,length(V)-1)*k;
    diag2(ly:ly:length(V)-1)=0;
    diag3(ly:ly:length(V)-1)=0;

    M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c diag diag1]);
    M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g diag2 diag3]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

    graph_resultat=zeros(length(y),length(t));
    for  j = 1 : length(t)-1

        b=M2*Psy;
        Psy = mldivide(M,b);
        Psy_mat=vec2mat(Psy,length(y))';
        graph_resultat(:,j)=abs(Psy_mat(:,floor(2*length(x)/3),1)).^2;

    end
    clear M M2;

    T(1,w)=toc()
end
% T = 1.0e+03 * [2.9902 0.1703 0.0263 0.0043 0.0007 0.0002]
       a =   1.152e-46  (-5.081e-46, 7.385e-46)
       b =      -4.118  (-4.314, -3.922)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DT
clear all

dt=1e-20;
Dt=[dt dt*2 dt*4 dt*8 dt*16 dt*32];

for w=1:length(Dt)
    
    tic
    
    hbar=1.0545718*10^-34;
    m=9.10938*10^-31;

    sig=8e-11;           
    lamda=5e-11;    
    kp=2*(pi/lamda)*[1,0];

    Imax=10;
    Jmax=10;

    dx=1e-11;
    dy=1e-11;
    dt=Dt(w);
    xf=2e-9;
    yf=1.4e-9;
    tf=1e-17;

    x=0:dx:xf;
    y=0:dy:yf; ly=length(y);
    t=0:dt:tf;
    x0=x(end)/3;
    y0=y(end)/2;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% vid 1

    [psy,norm] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

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

    diag=ones(1,length(V)-1)*d;
    diag1=ones(1,length(V)-1)*d;
    diag(ly:ly:length(V)-1)=0;
    diag1(ly:ly:length(V)-1)=0;

    diag2=ones(1,length(V)-1)*k;
    diag3=ones(1,length(V)-1)*k;
    diag2(ly:ly:length(V)-1)=0;
    diag3(ly:ly:length(V)-1)=0;

    M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c diag diag1]);
    M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g diag2 diag3]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

    graph_resultat=zeros(length(y),length(t));
    for  j = 1 : length(t)-1

        b=M2*Psy;
        Psy = mldivide(M,b);
        Psy_mat=vec2mat(Psy,length(y))';

    graph_resultat(:,j)=abs(Psy_mat(:,floor(2*length(x)/3),1)).^2;

    end
    clear M M2;

    T(2,w)=toc()
end

a =   5.766e-19 
b =       -1.03     

DT=linspace(DT(1),DT(end),100)
       
figure()
plot(Dt,T(2,:),'ro','LineWidth',1.5)
hold on
plot(DT,a.*DT.^b,'LineWidth',1.5)
xlabel('dt')
ylabel('Temps écoulé (s)')

% T =
% 
%   Columns 1 through 3
% 
%          0         0         0
%   224.6262  110.6219   53.9606
% 
%   Columns 4 through 5
% 
%          0         0
%    26.8821   11.6820
