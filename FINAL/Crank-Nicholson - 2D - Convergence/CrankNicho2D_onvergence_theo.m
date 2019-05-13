clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project          : Diffraction électronique
%% Function name    : Analyse de convergence de Crank-Nicholson 2D
%
%% Author           : Jean-Christophe G.-A., Thomas R., Nouffel S.
%% Date             : 5/1/2019
%
%% Description      : Obtention de l'ordre spatial et temporel de 
%%                  l'algorithme de Crank-Nicholson imprémenté
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Définitions

global ly hbar m sig kp x0 y0

% Constantes physiques
hbar=1.0545718*10^-34;
m=9.10938*10^-31;

% Paramètre du paquet d'onde
sig=8e-11;           
lamda=5e-11;    
kp=2*(pi/lamda)*[1,0];

% Discretisation spatiale et temporel
dx=5e-12;
dy=5e-12;
dt=1e-21;
xf=1.4e-9;
yf=0.8e-9;
tf=6e-17;

% Definition potentiel
pot=10000;
L_mur=1.1e-11;
D_fente=20e-11;
L_fente=3e-11;

% Initialisation Écrant
graph_resultat=zeros(length(y),length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse de convergence temporelle

%% Creation du domaine
x=0:dx:xf;
y=0:dy:yf; ly=length(y);
x0=x(end)/3;
y0=y(end)/2;

%% Creation des pas de temps à traité
Dt=[dt dt*2 dt*4 dt*8 dt*16 dt*32];
t1=0:Dt(end):tf;

%% Variables Mémoires
norm_t=zeros(length(Dt),length(t1));
norm_th_t=zeros(length(Dt),length(t1));
err_t=zeros(length(Dt),length(t1));

%% Propagation
for w=1:length(Dt)
    
    % Variables temporaires
    t=0:Dt(w):tf;
    norm=zeros(1,length(t));
    norm_th=zeros(1,length(t));
    err=zeros(1,length(t));
    
    % Initialisation du paquet d'onde et mémoire
    [psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);
    Psy=psy(:);
    
    % Initialisation des matrices du système
    v_mat=zeros(length(y),length(x));
    V=v_mat(:);
    [M,M2]=MM2(V,dx,dy,Dt(w));
    
    % Ménage de mémoire
    clear b c d diag diag2 diag3 f g k V v_mat
    
    % Initialisation du paquet d'onde théorique et calcul de la première
    % erreur
    [psy_th,norm_th(1)]=analy(x,y,t(1));
    err(1)=trapeze_2D(abs(squeeze(abs(psy_th))-squeeze(abs(psy))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    
    % Propagation
    for j = 1 : length(t)-1
        
        % Resolution du système linéraire
        b=M2*Psy;
        Psy = mldivide(M,b);
        psy=vec2mat(Psy,length(y))';
        
        % Calcul de la norme
        norm(j+1)=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        
        % Calcul de l'erreur
        [psy_th,norm_th(j+1)]=analy(x,y,t(j+1));
        err(j+1)=trapeze_2D(abs(squeeze(abs(psy_th))-abs(psy)),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);     
            
    end
    
    % Adaptation des données aux conteneurs mémoires
    norm_t(w,:)=norm(1:2^(length(Dt)-w):length(norm)); 
    norm_th_t(w,:)=norm_th(1:2^(length(Dt)-w):length(norm_th)); 
    err_t(w,:)=err(1:2^(length(Dt)-w):length(err));
    
end

%% Affichage Évolution norme numérique
figure()
hold on
for i=1:length(Dt)
    plot(t,norm_t(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dt(1)),sprintf('dx = dy = %d',Dt(2)),sprintf('dx = %d',Dt(3)))
title('Evolution de la norme de la resolution en fonction de DT')

%% Affichage Évolution norme théorique
figure()
hold on
for i=1:length(Dt)
    plot(t,norm_th_t(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dt(1)),sprintf('dx = dy = %d',Dt(2)),sprintf('dx = %d',Dt(3)))
title('Evolution de la norme de la solution analytique en fonction de DT')

%% Affichage Évolution erreur
color=['r' 'b' 'k' 'c' 'm'];
figure()
hold on
for i=1:length(Dt)
    plot(t,err_t(i,:),color(i))
end
legend(sprintf('\\Deltat=%.1e',Dt(1)),sprintf('\\Deltat=%.1e',Dt(2)),sprintf('\\Deltat=%.1e',Dt(3)),sprintf('\\Deltat=%.1e',Dt(4)),sprintf('\\Deltat=%.1e',Dt(5)),'Location','Best')


%% Trouver l'ordre de convergence (Naturel)
P_t=polyfit(log10(Dt),log10(err_t(:,end)'),1)

figure()
loglog(Dt,err_t(:,end),'ro')
hold on
loglog(Dt,10.^(P_t(1)*log10(Dt) +P_t(2)),'b','LineWidth',1.5)
xlabel('Log(dt)')
ylabel('log(Err)')
legend('Erreur après un temps fixe de propagation',sprintf('Régression y=%.3f*dt+%.3f',P_t(1),P_t(2)),'Location','Best')

%% Trouver l'ordre de convergence (Annulation Dx)
P_rigged=polyfit(log10(Dt(2:end))',log10(err_t(2:end,end)-err_t(1,end)),1)

figure()
loglog(Dt,err_t(:,end)-err_t(1,end),'ro')
hold on
loglog(Dt,10.^(P_rigged(1)*log10(Dt) +P_rigged(2)),'b','LineWidth',1.5)
xlabel('Log(dt)')
ylabel('log(Err)-log(Err_D_X)')
legend('Erreur après un temps fixe de propagation',sprintf('Régression y=%.3f*dt+%.3f',P_rigged(1),P_rigged(2)),'Location','Best')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse de convergence spatial

%% Creation du domaine
x=0:dx:xf;
y=0:dy:yf; ly=length(y);
x0=x(end)/3;
y0=y(end)/2;

%% Creation des pas d'espace à traiter
Dx=[dx dx*2 dx*4 dx*8];
x1=0:Dx(end):xf;

%% Variables Mémoires
norm_x=zeros(length(Dx),length(t));
norm_th_x=zeros(length(Dx),length(t));
err_x=zeros(length(Dx),length(t));

for w=1:length(Dx)

    
    % Variables temporaires
    x=0:Dx(w):xf;
    
    % Initialisation du paquet d'onde et mémoire
    [psy,norm_x(w,1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);
    Psy=psy(:);
    Psy_mat=psy;

    % Initialisation des matrices du système
    v_mat=zeros(length(y),length(x));
    V=v_mat(:);
    [M,M2]=MM2(V,Dx(w),dy,dt);
    
    % Initialisation du paquet d'onde théorique et calcul de la première
    % erreur
    [Psy_th,norm_th_x(w,1)]=analy(x,y,t(1));
    err_x(w,1)=trapeze_2D(abs(squeeze(abs(Psy_th))-squeeze(abs(Psy_mat))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    % Propagation
    for j = 1 : length(t)-1
        
        % Resolution du système linéraire
        b=M2*Psy;
        Psy = mldivide(M,b);
        Psy_mat=vec2mat(Psy,length(y))';

        % Calcul de la norme
        norm_x(w,j+1)=trapeze_2D(abs(squeeze(Psy_mat)).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
       
        % Calcul de l'erreur
        [Psy_th,norm_th_x(w,j+1)]=analy(x,y,t(j+1));
        err_x(w,j+1)=trapeze_2D(abs(squeeze(abs(Psy_th))-squeeze(abs(Psy_mat))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    end 
end

%% Affichage Évolution norme numérique
figure()
hold on
for i=1:length(Dx)
    plot(t,norm_x(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de la norme de la resolution en fonction de DX')

%% Affichage Évolution norme théorique
figure()
hold on
for i=1:length(Dx)
    plot(t,norm_th_x(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de la norme de la solution analyique en fonction de DX')

%% Affichage Évolution erreur
figure()
hold on
for i=1:length(Dx)
    plot(t,err_x(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de l''erreur en fonction de DX')

%% Trouver l'ordre de convergence (Naturel)
P_x=polyfit(log10(Dx),log10(err_x(:,end)'),1)

figure()
loglog(Dx,err_x(:,end),'ro')
hold on
loglog(Dx,10.^(P_x(1)*log10(Dx) +P_x(2)),'b','LineWidth',1.5)
legend('Erreur après un temps fixe de propagation',sprintf('Régression y=%.3f*dx+%.3f',P_x(1),P_x(2)),'Location','Best')
xlabel('Log(dx)')
ylabel('log(Erreur)')

