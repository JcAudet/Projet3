clear all;
clc;
close all;

format long

hbar=1;
m=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametres
sig=0.3;
k=20;
x0=2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Theo
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Etude convergence DX
% dx=0.05;
% Dx=[dx dx/2 dx/4 dx/8 dx/16 dx/32];
% dt=0.000001;
% t=0:dt:0.01;
% 
% norm_x=zeros( length(Dx),length(t) );
% norm_th=zeros(length(Dx),length(t));
% err_x=zeros(length(Dx),length(t));
% % Psy_mem_x = zeros(length(Dx),length(t),length(0:dx:5));
% 
% for i=1:length(Dx)
%     
%     x=0:Dx(i):5;
%     V=zeros(1,length(x));
%     
%     Psy=zeros(length(t),length(x));
%     Psy_th=zeros(length(t),length(x));
%     [Psy(1,:),norm_x(i,1)]=wp_ini(x,sig,k,x_0);
%     
%     tet=atan(2*hbar*t(1)/(m*sig^2));
%     delx=(sig/2)*sqrt(1+(4*hbar^2*t(1)^2)/(m^2*sig^4));
%     
%     Psy_th(1,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
%         exp(1i*k*(x-k*t(1)/2)).*exp(-(x-x_0-k*t(1)).^2 / (sig^2+2*1i*hbar*t(1)/m));
%     
%     for j=1:length(t)-1
%         Psy(j+1,:)=run_kutt_4(dt, Dx(i), Psy(j,:),V);
%         
%         tet=atan(2*hbar*t(j)/(m*sig^2));
%         delx=(sig/2)*sqrt(1+(4*hbar^2*t(j)^2)/(m^2*sig^4));
%     
%         Psy_th(j+1,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
%         exp(1i*k*(x-k*t(j)/2)).*exp(-(x-x_0-k*t(j)).^2 / (sig^2+2*1i*hbar*t(j)/m));
%         
%         norm_x(i,j+1)=trapeze(abs(Psy(j+1,:)).^2,x(1),x(end),length(x)-1);
%         err_x(i,j+1)=trapeze(abs(Psy(j+1,:)-Psy_th(j+1,:)),x(1),x(end),length(x)-1);
%     
%     end    
%     i
% end
% 
% figure()
% hold on
% for i=1:length(Dx)
%     plot(t,err_x(i,:))
% end
% legend(sprintf('dx = %f',Dx(1)),sprintf('dx = %f',Dx(2)),sprintf('dx = %f',Dx(3)),sprintf('dx = %f',Dx(4)),sprintf('dx = %f',Dx(5)),sprintf('dx = %f',Dx(6)))
% title('Erreur en fonction du temps et de la discretisation de l''espace')
% 
% P=polyfit(log10(Dx),log10(err_x(:,end)'),1)
% 
% figure()
% loglog(Dx,err_x(:,end))
% hold on
% loglog(Dx,10^(P(1)*log10(Dx)+P(2)))
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Etude convergence DT
% 
% dx=0.05;
% dt=0.000005;
% Dt=[dt dt/2 dt/4 dt/8 dt/16 dt/32];
% x=0:dx:5;
% 
% norm_t = zeros(length(Dt),length(0:dt:0.01));
% err_t = zeros(length(Dt),length(0:dt:0.01));
% 
% for i=1:length(Dt)
%     
%     t=0:Dt(i):0.01;
%     V=zeros(1,length(x));
%     err=zeros(1,length(t));
%     
%     Psy=zeros(length(t),length(x));
%     Psy_th=zeros(length(t),length(x));
%     [Psy(1,:),norm_t(i,1)]=wp_ini(x,sig,k,x_0);
%     
%     tet=atan(2*hbar*t(1)/(m*sig^2));
%     delx=(sig/2)*sqrt(1+(4*hbar^2*t(1)^2)/(m^2*sig^4));
%     
%     Psy_th(1,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
%         exp(1i*k*(x-hbar*k*t(1)/(2*m))).*exp(-(x-x_0-hbar*k*t(1)/m).^2 / (sig^2+2*1i*hbar*t(1)/m));
%     
%     err(1)=trapeze(abs(Psy(1,:)-Psy_th(1,:)),x(1),x(end),length(x)-1);
%     
%     for j=1:length(t)-1
%         Psy(j+1,:)=run_kutt_4(Dt(i), dx, Psy(j,:),V);
%         
%         tet=atan(2*hbar*t(j)/(m*sig^2));
%         delx=(sig/2)*sqrt(1+(4*hbar^2*t(j)^2)/(m^2*sig^4));
%     
%         Psy_th(j+1,:)=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet/2).*...
%         exp(1i*k*(x-hbar*k*t(j)/(2*m))).*exp(-(x-x_0-hbar*k*t(j)/m).^2 / (sig^2+2*1i*hbar*t(j)/m));
%         
%         norm_t(i,j+1)=trapeze(abs(Psy(j+1,:)).^2,x(1),x(end),length(x)-1);
%         err(j+1)=trapeze(abs(Psy(j+1,:)-Psy_th(j+1,:)),x(1),x(end),length(x)-1);
%                
%     end   
%    
% %     figure()
% %     plot(x,abs(Psy(1,:)))
% %     hold on
% %     plot(x,abs(Psy_th(1,:)))
% %     hold on
% %     plot(x,abs(Psy(end,:)))
% %     hold on
% %     plot(x,abs(Psy_th(end,:)))
% %     hold on
% %     plot(x,abs(Psy(floor(length(t)/2),:)))
% %     hold on
% %     plot(x,abs(Psy_th(floor(length(t)/2),:)))
%     
%     err_t(i,:)=err(1:2^(i-1):length(t));
%     
%     i
% end
% 
% 
% 
% figure()
% hold on
% for i=1:length(Dt)
%     plot(0:dt:0.01,err_t(i,:))
% end
% legend(sprintf('dt = %f',Dt(1)),sprintf('dt = %f',Dt(2)),sprintf('dt = %f',Dt(3)),sprintf('dt = %f',Dt(4)),sprintf('dt = %f',Dt(5)),sprintf('dt = %f',Dt(6)))
% title('Erreur en fonction du temps et de la discretisation du temps')
% 
% 
% P_t=polyfit(log10(Dt(1:end-1)),log10(err_t_int(:,end)'),1)
% 
% figure()
% loglog(Dt(1:end-1),err_t_int(:,end),'ro')
% hold on
% loglog(Dt(1:end-1),10.^(P_t(1)*log10(Dt(1:end-1)) +P_t(2)),'b','LineWidth',1.5)
% xlabel('Log(dt)')
% ylabel('log(Err)')
% legend('Erreur après un temps fixe de propagation',sprintf('Régression y=%.3f*dt+%.3f',P_t(1),P_t(2)),'Location','Best')
% 



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Non theo
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Etude convergence DX
dx=0.1;
Dx=[dx dx/2 dx/4 dx/8 dx/16 dx/32];
dt=0.000005;
t=0:dt:0.05;

norm_x=zeros( length(Dx),length(t) );
Psy_mem_x = zeros(length(Dx),length(t),length(0:dx:5));

for i=1:length(Dx)
    
    x=0:Dx(i):5;
    
    V=zeros(1,length(x));
    
    Psy=zeros(length(t),length(x));
    [Psy(1,:),norm_x(i,1)]=wp_ini(x,sig,k,x0);
    Psy_mem_x(i,1,:)=Psy(1,1:2^(i-1):length(Psy(i,:)));
    
    for j=1:length(t)-1
        Psy(j+1,:)=run_kutt_4(dt, Dx(i), Psy(j,:),V);
        Psy_mem_x(i,j+1,:)=Psy(j+1,1:2^(i-1):length(Psy(j,:)));
        norm_x(i,j+1)=trapeze(abs(Psy(j,:)).^2,x(1),x(end),length(Psy(j,:))-1);
    
    end
    i
end

err_x=zeros(length(Dx)-1,length(t),length(0:Dx(1):5));
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


P_x=polyfit(log10(Dx(1:end-1)),log10(err_x_int(:,end)'),1)

figure()
loglog(Dx(1:end-1),err_x_int(:,end),'ro')
hold on
loglog(Dx(1:end-1),10.^(P_x(1)*log10(Dx(1:end-1)) +P_x(2)),'b','LineWidth',1.5)
legend('Erreur après un temps fixe de propagation',sprintf('Régression y=%.3f*dx+%.3f',P_x(1),P_x(2)),'Location','Best')
xlabel('Log(dx)')
ylabel('log(Erreur)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Etude convergence DT

dx=0.01;
dt=0.0001;
Dt=[dt dt/2 dt/4 dt/8 dt/16 dt/32];
x=0:dx:5;

norm_t = zeros(length(Dt),length(0:dt:0.05));
Psy_mem_t = zeros(length(Dt),length(0:dt:0.05),length(x));

for i=1:length(Dt)
    
    t=0:Dt(i):0.05;
    
    V=zeros(1,length(x));
    
    Psy=zeros(length(t),length(x));
    [Psy(1,:),norm_t(i,1)]=wp_ini(x,sig,k,x0);
    Psy_mem_t(i,1,:)=Psy(1,:);
    
    for j=1:length(t)-1
        Psy(j+1,:)=run_kutt_4(Dt(i),dx,Psy(j,:),V);
        norm_t(i,j+1)=trapeze(abs(Psy(j+1,:)).^2,x(1),x(end),length(Psy(j+1,:))-1);
    end
    
    Psy_mem_t(i,:,:)=Psy(1:2^(i-1):length(Psy(:,1)),:);
    
    i
end

err_t=zeros(length(Dt)-1,length(0:dt:0.05),length(x));
err_t_int = zeros(length(Dt)-1,length(0:dt:0.05));

% Calcul erreur
for i=1:length(Dt)-1
    
    err_t(i,:,:)=abs(Psy_mem_t(i,:,:)-Psy_mem_t(i+1,:,:));
    
    for j = 1:length(0:dt:0.05)
        
        err_t_int(i,j)=trapeze(err_t(i,j,:),0,5,length(err_t(i,j,:))-1);
        
    end
end


figure()
hold on
for i=1:length(Dt)-1
    plot(0:dt:0.05,err_t_int(i,:))
end
legend(sprintf('dt = %f',Dt(1)),sprintf('dt = %f',Dt(2)),sprintf('dt = %f',Dt(3)),sprintf('dt = %f',Dt(4)),sprintf('dt = %f',Dt(5)))
title('Erreur en fonction du temps et de la discretisation du temps')

P_t=polyfit(log10(Dt(1:end-1)),log10(err_t_int(:,end)'),1)

figure()
loglog(Dt(1:end-1),err_t_int(:,end),'ro')
hold on
loglog(Dt(1:end-1),10.^(P_t(1)*log10(Dt(1:end-1)) +P_t(2)),'b','LineWidth',1.5)
xlabel('Log(dt)')
ylabel('log(Erreur)')
legend('Erreur après un temps fixe de propagation',sprintf('Régression y=%.3f*dt+%.3f',P_t(1),P_t(2)),'Location','Best')


