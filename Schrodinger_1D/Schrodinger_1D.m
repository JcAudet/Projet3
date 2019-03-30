clear all;
clc;
close all;

global dx; 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation (OK)

dx=0.01;
dt_0=0.0001;
dt_1=0.00005;
dt_2=0.000025;
t_0=0:dt_0:0.3;
t_1=0:dt_1:0.3;
t_2=0:dt_2:0.3;
x=0:dx:10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter of wave packet (OK, voir valeurs reelles)

sig=0.3;
k=20;
x_0=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation (OK)

psy=exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);
psy(1)=0;
psy(end)=0;

norm_0=zeros(1,length(t_1));
norm_1=zeros(1,length(t_1));
norm_2=zeros(1,length(t_1));

% Normalisation

psy_abs=abs(psy);
norme=trapeze(psy_abs,x(1),x(end),length(psy_abs)-1);
norm_0=norme;norm_1=norme;norm_2=norme;
psy=psy./norme;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse Freq

[f,Phi,STD]=fft_wp(x,psy);
f_max=f(find(Phi==max(Phi)))+4*STD;
c_max=2*pi*f_max/k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Critere Stabilite

C=dt_1/c_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of propagation matrix

Psy=zeros(length(t_1),length(x));Psy_0=zeros(length(t_1),length(x));Psy_2=zeros(length(t_1),length(x));
Psy(1,:)=psy;Psy_0(1,:)=psy;Psy_2(1,:)=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potentiel
V=zeros(1,length(x));
% sig_g=0.01;
% pot = -1000 * 1i * exp( -((x-9).^2)/(2*sig_g^2) ) - 10000 * 1i * exp( -((x-9.8).^2)/(2*sig_g^2) );
% V = - 100000000 * 1i * exp( -((x-4).^2)/(2*sig_g^2) );
% pot=(x-5).^2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Solution reelle
% 
% % Déplacement temporel
% 
% a=1;
% m=9.109*10^(-31);
% h_bar=6.62607004*10^(-34)/(2*pi);
% 
% norme_re=zeros(1,length(t));
% Psy_re=zeros(length(t),length(x));
% Psy_re(1,:)=abs(psy);
% 
% figure()
% for i=1:length(t)-1
%     
%     y=sqrt(1+4*h_bar^2*t(i+1)^2/(m^2*a^4));
%     
%     for j=1:length(x)
%         
%         Psy_re(i+1,j)= sqrt(2/(pi*a^2)) * 1/y * exp(-2/(a*y)^2 * ( x(j)-x_0-h_bar * k * t(i+1)/m )^2);
%     
%     end 
%     
%     plot(x,Psy_re(i+1,:))
%     norme_re(i)=trapeze(abs(psy),x(1),x(end),length(psy)-1);
%     title(sprintf('Norme: %.10f',norme_re(i)))
%     pause(0.002)
%     
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Crank nich
%
% % Definition
% 
% a=1+1i*dt/(dx^2)+1i*dt*V/2;
% b=ones(1,length(x)-1)*(-1i*dt/(2*dx^2));
% g=1-1i*dt/(dx^2)-1i*dt*V/2;
% 
% C=sparse( full( gallery('tridiag',-b,g,-b)));
% 
% A=sparse( full(gallery('tridiag',b,a,b)) );

% % Calcul
%
% tic
% figure()
% for k=1:length(t)
%     
%     Psy(k+1,:)=crank_nicholson(Psy(k,:));
%     norm(k)=trapeze(abs(Psy(k,:)),x(1),x(end),length(Psy(k,:))-1);
%     
%     plot(x, real(Psy(k,:)),'b');
%     hold on
%     plot(x,imag(Psy(k,:)),'r');
%     hold on
%     plot(x,abs(Psy(k,:)))
%     hold on
%     ylim([-10 10]);
%     plot(x,V)
%     hold off
%     title( sprintf('t = %.5f , Norme= %.6f', t(k), norm(k)));
%     pause(0.01);
% end
% toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runge Kutta O4

tic
% figure()
% for i=1:length(t_0)-1
%     Psy_0(i+1,:)=run_kutt_4(t_0,Psy_0(i,:),V);
%     norm_0(i+1)=trapeze(abs(Psy_0(i,:)),x(1),x(end),length(Psy_0(i,:))-1);
%     i
% end
for j=1:length(t_1)-1
    Psy(j+1,:)=run_kutt_4(t_1,Psy(j,:),V);
    norm_1(j+1)=trapeze(abs(Psy(j,:)),x(1),x(end),length(Psy(j,:))-1);
    j
end
% for k=1:length(t_2)-1
%     Psy_2(k+1,:)=run_kutt_4(t_2,Psy_2(k,:),V);
%     norm_2(k+1)=trapeze(abs(Psy_2(k,:)),x(1),x(end),length(Psy_2(k,:))-1);
%     k
% end
toc

% for k=1:length(t_0)
%     
%     subplot(2,3,1)
%     plot(x,abs(Psy_0(k,:)))
%     subplot(2,3,2)
%     plot(x,abs(Psy(k*2,:)))
%     subplot(2,3,3)
%     plot(x,abs(Psy_2(k*4,:)))
%     subplot(2,3,4)
%     plot(t_1(1:k),norm_0(1:k))
%     subplot(2,3,5)
%     plot(t_1(1:k*2),norm_1(1:k*2))
%     subplot(2,3,6)
%     plot(t_1(1:k*4),norm_2(1:k*4))
%     
%     pause(0.01)
%     
% end
%     plot(x, real(Psy(k,:)),'b');
%     hold on
%     plot(x, Psy_re(k,:),'r');
%     hold on
%     plot(x,abs(Psy(k,:)))
%     hold on
%     ylim([-10 10]);
%     plot(x,V)
%     hold off
%     title( sprintf('t = %.5f , Norme= %.6f', t(k), norm(k)));
%     pause(0.01);
% end
% toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure 
hold on
plot(t_0,norm_0)
plot(t_1,norm_1)
plot(t_2,norm_2)
title('Norm evolution')
xlabel('Time')
ylabel('Norm')