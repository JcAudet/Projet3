clear all;
clc;
close all;

format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametres
sig=0.3;
k=50;
x0=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Etude convergence DX
Dx=[0.05 0.01 0.005 0.001];
dt=0.0001;
t=0:dt:0.05;

norm_x=zeros(length(Dx),length(t));

for i=1:length(Dx)
    
    x=0:Dx(i):5;
    
    V=zeros(1,length(x));
    
    Psy=zeros(length(t),length(x));
    [Psy(1,:),norm_x(i,1)]=wp_ini(x,sig,k,x0);
    
    for j=1:length(t)-1
        Psy(j+1,:)=run_kutt_4(dt, Dx(1), Psy(j,:),V);
        norm_x(i,j+1)=trapeze(abs(Psy(j,:)).^2,x(1),x(end),length(Psy(j,:))-1);
    end
    i
end

figure()
hold on
for i=1:length(Dx)
    plot(t,norm_x(i,:)-1)
end
legend(sprintf('dx = %f',Dx(1)),sprintf('dx = %f',Dx(2)),sprintf('dx = %f',Dx(3)),sprintf('dx = %f',Dx(4)))
title('Erreur en fonction du temps et de la discretisation de l''espace')

figure()
loglog(Dx,abs(norm_x(:,end)-1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Etude convergence DT

dx=0.01;
Dt=[0.0001 0.00005 0.00001 0.000005];
x=0:dx:5;

norm_t=zeros(length(Dt),length(0:Dt(end):0.05));

for i=1:length(Dt)
    
    t=0:Dt(i):0.05;
    
    V=zeros(1,length(x));
    
    Psy=zeros(length(t),length(x));
    [Psy(1,:),norm_t(i,1)]=wp_ini(x,sig,k,x0);
    
    for j=1:length(t)-1
        Psy(j+1,:)=run_kutt_4(Dt(i),dx,Psy(j,:),V);
        norm_t(i,j+1)=trapeze(abs(Psy(j,:)).^2,x(1),x(end),length(Psy(j,:))-1);
    end
    i
end

figure()
norme_t_end=zeros(1,length(Dt));
hold on
for i=1:length(Dt)
    norme=norm_t(i,:);
    norme=norme(norme>0.5);
    norme_t_end(i)=norme(end);
    plot(0:Dt(i):0.05,norme-1)
end
legend(sprintf('dt = %f',Dt(1)),sprintf('dt = %f',Dt(2)),sprintf('dt = %f',Dt(3)),sprintf('dt = %f',Dt(4)))
title('Erreur en fonction du temps et du timestep')

figure()
loglog(Dt,abs(norme_t_end-1))


