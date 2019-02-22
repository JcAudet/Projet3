%J'ai modifier : 
%ligne 92 et 100 pour les C.F.
%J'ai ecrit M avec sparse. Et j'ai créée une nouvelle variable qui sappel
%exigence_memoire, elle contient  16*(nnz(L)+nnz(U)) comme dans les
%consignes
clc;clear all;
% Distribution de la température dans un appartement d'un immeuble de plusieurs étages

% Équation de transfert de chaleur:
% k*(d^2 T(x,y)/dx^2 + d^2 T(x,y)/dy^2)+S=0

% Conditions aux limites:

% (1) Condition convective (de Robin) à x=0 et à x=Lx (faces externes du mur):
% -k*dT(x=0,y)/dx=-h*(T-Ta)
% -k*dT(x=L,y)/dx=h*(T-Ta)
Ta=-10; %oC

% (2) Condition de Dirichlet sur le plafond et sur ??le plancher
% T(x, y=0 ou y=Ly)=Tp
Tp=-10; %oC

% Dimensions d'appartement
Lx=4; %[m]
Ly=2.4;  %[m]

% Parametres d'un mur d'isolation thermique
Lm=0.4; %m ; Épaisseur du mur en brique
km=0.85;%W/(m*K); La conductivité thermique de la brique
h=10; %W/(m^2*K); Coefficient de transfert thermique sur les surfaces extérieures du mur

% Paramètres de l'air qui remplit l'appartement
ka=0.024;

d_ar=[];tini_ar=[];tinv_ar=[];mem_ar=[];Tm_ar=[];
for fact=(1/12).*[16,8,4,2,1]
    d=0.1*fact; %Pas de discrétisation en [m]
    d_ar=[d_ar d]
    Nx=round(Lx/d)+1;
    Ny=round(Ly/d)+1;
    
    tic
    % Initialisation de la source de chaleur et de la conductivité thermique
    S=zeros(Ny,Nx); k=zeros(Ny,Nx);
    for i=1:Ny
        y=(i-1)*d;
        for j=1:Nx
            x=(j-1)*d;
            
            % Sourse volumique de chaleur q[W/m^3] d'épaisseur dL.
            % La source est intégrée dans la partie intérieure du mur a x=0 et
            % il occupe le tiers du mur dans la direction verticale
            dL=0.1;
            q=1e4;% W/m^3;
            if (x<=Lm)&&(y<=Ly/3)
                % À l'intérieur de l'élément chauffant
                S(i,j)=q*exp(-((x-Lm)/dL).^2);
            else
                % À l'extérieur de l'élément chauffant
                S(i,j)=0;
            end
            
            % l'espace de vie de l'appartement est délimitée par
            % les parois d'épaisseur Lm de tous les quatre côtés
            if (x<=Lm)||(x>=(Lx-Lm))||(y<=Lm)||(y>=(Ly-Lm))
                % À l'intérieur du mur
                k(i,j)=km;
            else
                % À l'intérieurde de l'appartement
                k(i,j)=ka;
            end
        end
    end
    
    M=sparse(Nx*Ny,Nx*Ny);
    b=zeros(Nx*Ny,1);
    
    for i=1:Ny
        for j=1:Nx
            % remplir la ligne pl de la matrice M
            pl=(i-1)*Nx+j;
            
            if ((i>1)&(i<Ny))&((j>1)&(j<Nx))
                % noeud qui est strictement à l'intérieur de la cellule de simulation
                pc=pl;M(pl,pc)=-4; % contribution de noeud (i,j)
                pc=(i-1)*Nx+j-1;M(pl,pc)=1; % contribution de noeud (i,j-1)
                pc=(i-1)*Nx+j+1;M(pl,pc)=1; % contribution de noeud (i,j+1)
                pc=(i-2)*Nx+j;M(pl,pc)=1; % contribution de noeud (i-1,j)
                pc=(i)*Nx+j;M(pl,pc)=1; % contribution de noeud (i+1,j)
                b(pl)=-d^2*S(i,j)/k(i,j);
            elseif (i==1)
                 %CODE MODIFIER : 
                pc=pl;M(pl,pc)=3+2*d*h/k(i,j); % contribution de noeud (1,j)
                pc=(i-1)*Nx+j+1;M(pl,pc)=-4; % contribution de noeud (2,j)
                pc=(i-1)*Nx+j+2;M(pl,pc)=1; % contribution de noeud (3,j)
                b(pl)=2*d*h*Ta/k(i,j);
                
            elseif (i==Ny)
               
                %CODE MODIFIER : 
                pc=pl;M(pl,pc)=3+2*d*h/k(i,j); % contribution de noeud (Ny,j)
                pc=(i-1)*Nx+j-1;M(pl,pc)=-4; % contribution de noeud (Ny-1,j)
                pc=(i-1)*Nx+j-2;M(pl,pc)=1; % contribution de noeud (Ny-2,j)
                b(pl)=2*d*h*Ta/k(i,j);
                
            elseif (j==1)
                % noeud à la surface externe du mur x=0
                pc=pl;M(pl,pc)=3+2*d*h/k(i,j); % contribution de noeud (i,1)
                pc=(i-1)*Nx+j+1;M(pl,pc)=-4; % contribution de noeud (i,2)
                pc=(i-1)*Nx+j+2;M(pl,pc)=1; % contribution de noeud (i,3)
                b(pl)=2*d*h*Ta/k(i,j);
            elseif (j==Nx)
                % noeud à la surface externe du mur x=Nx
                pc=pl;M(pl,pc)=3+2*d*h/k(i,j); % contribution de noeud (i,Nx)
                pc=(i-1)*Nx+j-1;M(pl,pc)=-4; % contribution de noeud (i,Nx-1)
                pc=(i-1)*Nx+j-2;M(pl,pc)=1; % contribution de noeud (i,Nx-2)
                b(pl)=2*d*h*Ta/k(i,j);
            else
                display('Erreur dans la définition de la matrice de coefficients');
            end
        end
    end
    tini_ar=[tini_ar toc]
    
    tic
    %T=M\b;
    [L,U]=lu(M);T=U\(L\b);
    tinv_ar=[tinv_ar toc]
    
    exigence_memoire=16*(nnz(L)+nnz(U));
    mem_ar=[mem_ar 16*(nnz(L)+nnz(U))]
    
    Tr=reshape(T,Nx,Ny)';
    
    Tm_ar=[Tm_ar Tr(round(Ly/d/2+1),round(Lx/d/2+1))]
    
end

figure(1)
h=pcolor((0:d:Lx),(0:d:Ly),S);set(h,'LineStyle','none')
colorbar
xlabel('x [m]'); ylabel('y [m]'); title('S(x,y) [W/m^3]')

figure(2)
h=pcolor((0:d:Lx),(0:d:Ly),k);set(h,'LineStyle','none')
colorbar
xlabel('x [m]'); ylabel('y [m]'); title('k(x,y) [W/(m^2\cdotK)]')

figure(3)
h=pcolor((0:d:Lx),(0:d:Ly),Tr);set(h,'LineStyle','none')
colormap(hot);
colorbar
xlabel('x [m]'); ylabel('y [m]'); title('T(x,y) [K]')

figure(4)
loglog(d_ar,mem_ar/1024^3,':o')
axis([1e-3 1 1e-4 128])
axis square
xlabel('\delta')
ylabel('mem (Gb)')

figure(5)
loglog(d_ar(2:end),abs(Tm_ar(2:end)-Tm_ar(1:end-1)),':o')
axis tight
xlabel('\delta')
ylabel('Err')

figure(6)
loglog(d_ar,tini_ar,':o')
axis([1e-3 1 1e-3 10000])
axis square
xlabel('\delta')
ylabel('t_{ini} (s)')

figure(7)
loglog(d_ar,tinv_ar,':o')
axis([1e-3 1 1e-3 10000])
axis square
xlabel('\delta')
ylabel('t_{inv} (s)')




%%QUESTION 1 : Err_precis et Tm_precis

Err_precis=min(abs(Tm_ar(2:end)-Tm_ar(1:end-1)));
pos=find(abs(Tm_ar(2:end)-Tm_ar(1:end-1))== Err_precis);
Tm_precis=Tm_ar(pos);

disp(Err_precis)
disp(Tm_precis)


%QUESTION 2 : A___ et P___  , ___ peut etre tinv, tini, err, mem
%%Moindre carre pour trouver les parametre de chacun des graph

nsx=transpose(log(mem_ar/1024^3));
nsf=transpose(log(d_ar));
C=[1 nsx(1);1 nsx(2);1 nsx(3);1 nsx(4);1 nsx(5)];
b=nsf;
XmoindreCarre=(inv(transpose(C)*C))*(transpose(C))*b;

Amem=10^XmoindreCarre(1);
Pmem=XmoindreCarre(2);



nsx=transpose(log(tini_ar));
nsf=transpose(log(d_ar));
C=[1 nsx(1);1 nsx(2);1 nsx(3);1 nsx(4);1 nsx(5)];
b=nsf;
XmoindreCarre=(inv(transpose(C)*C))*(transpose(C))*b;

Atini=10^XmoindreCarre(1);
Ptini=XmoindreCarre(2);


nsx=transpose(log(tinv_ar));
nsf=transpose(log(d_ar));
C=[1 nsx(1);1 nsx(2);1 nsx(3);1 nsx(4);1 nsx(5)];
b=nsf;
XmoindreCarre=(inv(transpose(C)*C))*(transpose(C))*b;

Atinv=10^XmoindreCarre(1);
Ptinv=XmoindreCarre(2);



nsx=transpose(log(abs(Tm_ar(2:end)-Tm_ar(1:end-1))));
nsf=transpose(log(d_ar(2:end)));
C=[1 nsx(1);1 nsx(2);1 nsx(3);1 nsx(4)];
b=nsf;
XmoindreCarre=(inv(transpose(C)*C))*(transpose(C))*b;

Aerr=10^XmoindreCarre(1);
Perr=XmoindreCarre(2);



%%QUESTION 3 avec 128GB on trouve delta, et ensuite on trouve le reste
delta_128GB = exp((log(128)-log(Amem))/Pmem);

tini_128GB = exp(log(Atini) + Ptini*log(delta_128GB));
tinv_128GB = exp(log(Atinv) + Ptinv*log(delta_128GB));
Err_128GB = exp(log(Aerr) + Perr*log(delta_128GB));




