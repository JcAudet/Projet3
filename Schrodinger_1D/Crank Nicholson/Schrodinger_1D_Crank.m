clear all;
clc;
close all;

global dx C A sig x_0 k hbar m; 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation (OK)

hbar=1;
m=1;

dx=0.1;
dt=0.00001;
t=0:dt:0.05;
x=0:dx:5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation 

% Parametres
sig=0.3;
k=20;
x_0=2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation (OK)

[psy,norme]=wp_ini(x,sig,k,x_0);

norm=zeros(1,length(t));
norm(1)=norme;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of propagation matrix

Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potentiel

V=zeros(1,length(x));
% sig_g=0.05;
% V = 1000 * exp( -((x-3.5).^2)/(2*sig_g^2) );
% V=zeros(1,length(x));
% V(350:355)=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crank nich
% Definition

a=1+1i*dt/(2*dx^2)+1i*V*dt/2;
b_C=ones(1,length(x)-1)*(1i*dt/(4*dx^2));
b_A=ones(1,length(x)-3)*(1i*dt/(4*dx^2));
g=1-1i*dt/(2*dx^2)-1i*V*dt/2;

C=sparse( full( gallery('tridiag',b_C,g,b_C)));
C=C(2:end-1,:);
C(1,1)=C(1,1)*2; C(end,end)=C(end,end)*2;

A=sparse( full(gallery('tridiag',-b_A,a(2:end-1),-b_A)) );


% Calcul
tic
gcf=figure();
for j=1:length(t)-1
    
    D = C * transpose(Psy(j,:));
    Psy(j+1,2:end-1)= mldivide(A,D);
    
    plot(x,abs(Psy(j,:)).^2)
    hold off
    
    pause(0.0001)
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% figure()
% subplot(1,3,1)
% plot(x,abs(Psy(1,:).^2),'r','linewidth',1)
% hold on
% plot(x,abs(psy(1,:)).^2,'b--','linewidth',1)
% ylim([0 3]);
% ylabel('$ |\Psi|^2 $','Interpreter','latex')
% title(sprintf('Temps = %.2f seconde',t(1)))
% subplot(1,3,2)
% plot(x,abs(Psy(2500,:).^2),'r','linewidth',1)
% hold on
% plot(x,abs(psy(2500,:)).^2,'b--','linewidth',1)
% ylim([0 3]);
% xlabel('x(u.a.)')
% title(sprintf('Temps = %.2f seconde',t(2500)))
% subplot(1,3,3)
% plot(x,abs(Psy(4000,:).^2),'r','linewidth',1)
% hold on
% plot(x,abs(psy(4000,:)).^2,'b--','linewidth',1)
% ylim([0 3]);
% legend('Résolution numérique','Solution analytique')
% title(sprintf('Temps = %.2f seconde',t(4000)))



