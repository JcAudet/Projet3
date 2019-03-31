clear all;
clc;
close all;

format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation (OK)

dx_0=0.05;
dx=0.01;
dx_2=0.007;
dt_0=0.0001;
dt=0.00005;
dt_2=0.000025;
t_0=0:dt_0:0.085;
t=0:dt:0.085;
t_2=0:dt_2:0.085;
x_0=0:dx_0:5;
x=0:dx:5;
x_2=0:dx_2:5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter of wave packet (OK, voir valeurs reelles)

sig=0.3;
k=20;
x0=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave packet initialisation (OK)

[psy_0,norme_0]=wp_ini(x_0,sig,k,x0);
[psy,norme]=wp_ini(x,sig,k,x0);
[psy_2,norme_2]=wp_ini(x_2,sig,k,x0);

norm_t_0=zeros(1,length(t_0));
norm_t_1=zeros(1,length(t));
norm_t_2=zeros(1,length(t_2));

norm_x_0=zeros(1,length(t));
norm_x_1=zeros(1,length(t));
norm_x_2=zeros(1,length(t));

norm_t_0(1)=norme;norm_t_1(1)=norme;norm_t_2(1)=norme;
norm_x_0(1)=norme_0;norm_x_1(1)=norme;norm_x_2(1)=norme_2;

%% Initialisation of propagation matrix

Psy_t_0=zeros(length(t_0),length(x));
Psy_t_1=zeros(length(t),length(x));
Psy_t_2=zeros(length(t_2),length(x));
Psy_t_0(1,:)=psy;
Psy_t_1(1,:)=psy;
Psy_t_2(1,:)=psy;

Psy_x_0=zeros(length(t),length(x_0));
Psy_x_1=zeros(length(t),length(x));
Psy_x_2=zeros(length(t),length(x_2));
Psy_x_0(1,:)=psy_0;
Psy_x_1(1,:)=psy;
Psy_x_2(1,:)=psy_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potentiel
V=zeros(1,length(x));
V_0=zeros(1,length(x_0));
V_2=zeros(1,length(x_2));
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

% dt
tic
for i=1:length(t_0)-1
    Psy_t_0(i+1,:)=run_kutt_4(dt_0, dx,Psy_t_0(i,:),V);
    norm_t_0(i+1)=trapeze(abs(Psy_t_0(i,:)),x(1),x(end),length(Psy_t_0(i,:))-1);
    i
end
for j=1:length(t)-1
    Psy_t_1(j+1,:)=run_kutt_4(dt, dx, Psy_t_1(j,:),V);
    norm_t_1(j+1)=trapeze(abs(Psy_t_1(j,:)),x(1),x(end),length(Psy_t_1(j,:))-1);
    j
end
% figure()
for k=1:length(t_2)-1
    Psy_t_2(k+1,:)=run_kutt_4(dt_2, dx, Psy_t_2(k,:),V);
    norm_t_2(k+1)=trapeze(abs(Psy_t_2(k,:)),x(1),x(end),length(Psy_t_2(k,:))-1);
    
%     plot(x,Psy_t_2(k,:))
%     title( sprintf('t = %.5f , Norme= %.6f', t_2(k), norm_t_2(k)));
%     hold off
%     pause(0.001)
    k
end
toc

% dx
tic
% figure()
for ii=1:length(t)-1
    Psy_x_0(ii+1,:)=run_kutt_4(dt, dx_0, Psy_x_0(ii,:),V_0);
    norm_x_0(ii+1)=trapeze(abs(Psy_x_0(ii,:)),x_0(1),x_0(end),length(Psy_x_0(ii,:))-1);
    
%     plot(x,Psy_x_0(ii,:))
%     title( sprintf('t = %.5f , Norme= %.6f', t(k), norm_x_0(ii)));
%     hold off
%     pause(0.001)
    
    ii
end
for jj=1:length(t)-1
    Psy_x_1(jj+1,:)=run_kutt_4(dt, dx, Psy_x_1(jj,:),V);
    norm_x_1(jj+1)=trapeze(abs(Psy_x_1(jj,:)),x(1),x(end),length(Psy_x_1(jj,:))-1);
    jj
end
for kk=1:length(t)-1
    Psy_x_2(kk+1,:)=run_kutt_4(dt, dx_2, Psy_x_2(kk,:),V_2);
    norm_x_2(kk+1)=trapeze(abs(Psy_x_2(kk,:)),x_2(1),x_2(end),length(Psy_x_2(kk,:))-1);
    kk
end
toc


% for k=1:length(t_0)
%     
%     subplot(2,3,1)
%     plot(x_2,abs(Psy_t_0(k,:)))
%     subplot(2,3,2)
%     plot(x_2,abs(Psy_t_1(k*2,:)))
%     subplot(2,3,3)
%     plot(x_2,abs(Psy_t_2(k*4,:)))
%     subplot(2,3,4)
%     plot(t(1:k),norm_t_0(1:k))
%     subplot(2,3,5)
%     plot(t(1:k*2),norm_t_1(1:k*2))
%     subplot(2,3,6)
%     plot(t(1:k*4),norm_t_2(1:k*4))
%     
%     pause(0.01)
%     
% end

figure()

subplot(2,3,1)
plot(t_0,norm_t_0)
subplot(2,3,2)
plot(t,norm_t_1)
subplot(2,3,3)
plot(t_2,norm_t_2)
subplot(2,3,4)
plot(t,norm_x_0)
subplot(2,3,5)
plot(t,norm_x_1)
subplot(2,3,6)
plot(t,norm_x_2)

T=[dt_0 dt dt_2];
T_e=[norm_t_0(end)-1 norm_t_1(end)-1 norm_t_2(end)-1];
X=[dx_0 dx dx_2];
X_e=[norm_x_0(end)-1 norm_x_1(end)-1 norm_x_2(end)-1]

figure()
subplot(1,2,1)
hold on
loglog(T,T_e,'.b')
subplot(1,2,2)
hold on
loglog(X,X_e,'.b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

