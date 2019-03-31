clear all;
clc;
close all;

global dx; 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation (OK)

dx=0.01;
dt=0.00005;
t=0:dt:0.3;
x=0:dx:5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter of wave packet (OK, voir valeurs reelles)

sig=0.3;
k=20;
x_0=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation (OK)

[psy,norme]=wp_ini(x,sig,k,x_0);

norm=zeros(1,length(t));
norm(1)=norme;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse Freq et Critere stabilite

[f,Phi,STD]=fft_wp(x,psy);
f_max=f(find(Phi==max(Phi)))+4*STD;
c_max=2*pi*f_max/k;
C=dt/c_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of propagation matrix

Psy=zeros(length(t),length(x));
Psy(1,:)=psy;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runge Kutta O4

tic
figure()
for k=1:length(t)-1
    
    Psy(k+1,:)=run_kutt_4(t,Psy(k,:),V);
    norm(k+1)=trapeze(abs(Psy(k,:)),x(1),x(end),length(Psy(k,:))-1);
    
    subplot(2,1,1)
    plot(x, real(Psy(k,:)),'b');
    hold on
    plot(x,abs(Psy(k,:)))
    hold on
    ylim([-10 10]);
    plot(x,V)
    hold off
    title( sprintf('t = %.5f , Norme= %.6f', t(k), norm(k)));
    pause(0.0001);
    
    subplot(2,1,2)
    plot(t(1:k),norm(1:k))
    title('Norm evolution')
    xlabel('Time')
    ylabel('Norm')
    xlim([0,0.3])
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%