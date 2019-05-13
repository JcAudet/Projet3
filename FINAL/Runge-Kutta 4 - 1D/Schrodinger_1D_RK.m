clear all;
clc;
close all;

global dx; 
format long

hbar=1;
m=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

dx=0.01;
dt=0.00001;
t=0:dt:0.05;
x=0:dx:5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation 

% Parametres
sig=0.3;
k=20;
x_0=2;

% Creation
[psy,norme]=wp_ini(x,sig,k,x_0);
norm=zeros(1,length(t));
norm(1)=norme;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potentiel

V=zeros(1,length(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runge Kutta O4

% Initialisation of propagation matrix
Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

Psy_th=zeros(length(t),length(x));
err=zeros(1,length(t));

tet=atan(2*hbar*t(1)/(m*sig^2));
delx=(sig/2)*sqrt(1+(4*hbar^2*t(1)^2)/(m^2*sig^4));
Psy_th(1,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
        exp(1i*k*(x-k*t(1)/2)).*exp(-(x-x_0-k*t(1)).^2 / (sig^2+2*1i*hbar*t(1)/m));

err(1)=trapeze(abs(abs(Psy(1,:))-abs(Psy_th(1,:))),x(1),x(end),length(x)-1);

tic
for i=1:length(t)-1
    
    Psy(i+1,:)=run_kutt_4(dt,dx,Psy(i,:),V);
    norm(i+1)=trapeze((abs(Psy(i,:))).^2,x(1),x(end),length(Psy(i,:))-1);
    
    tet=atan(2*hbar*t(i)/(m*sig^2));
    delx=(sig/2)*sqrt(1+(4*hbar^2*t(i)^2)/(m^2*sig^4));
    
    Psy_th(i+1,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
        exp(1i*k*(x-k*t(i)/2)).*exp(-(x-x_0-k*t(i)).^2 / (sig^2+2*1i*hbar*t(i)/m));
    
    err(i+1)=trapeze(abs(Psy(i,:)-Psy_th(i,:)),x(1),x(end),length(x)-1);
    
    
end
toc

figure()
subplot(1,3,1)
plot(x,abs(Psy(1,:).^2),'r','linewidth',1)
hold on
plot(x,abs(Psy_th(1,:)).^2,'b--','linewidth',1)
ylim([0 3]);
ylabel('$ |\Psi|^2 $','Interpreter','latex')
title(sprintf('Temps = %.2f seconde',t(1)))
subplot(1,3,2)
plot(x,abs(Psy(2500,:).^2),'r','linewidth',1)
hold on
plot(x,abs(Psy_th(2500,:)).^2,'b--','linewidth',1)
ylim([0 3]);
xlabel('x(u.a.)')
title(sprintf('Temps = %.2f seconde',t(2500)))
subplot(1,3,3)
plot(x,abs(Psy(4000,:).^2),'r','linewidth',1)
hold on
plot(x,abs(Psy_th(4000,:)).^2,'b--','linewidth',1)
ylim([0 3]);
legend('Résolution numérique','Solution analytique')
title(sprintf('Temps = %.2f seconde',t(4000)))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%