clear all;
clc;
close all;

format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametres
sig=0.3;
k=20;
x0=2;

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


figure()
loglog(Dx(1:end-1),err_x_int(:,end))

P=polyfit(log(Dx(1:end-1)),log(err_x_int(:,end)'),1)

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


figure()
loglog(Dt(1:end-1),err_t_int(:,end))

P=polyfit(log(Dt(1:end-1)),log(err_t_int(:,end)'),1)

