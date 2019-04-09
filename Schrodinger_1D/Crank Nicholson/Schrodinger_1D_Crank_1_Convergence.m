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
t=0:dt:0.01;

norm_x=zeros( length(Dx),length(t) );
Psy_mem_x = zeros(length(Dx),length(t),length(0:dx:5));

figure()
hold on
for i=1:length(Dx)
    x=0:Dx(i):5;
    
    V=zeros(1,length(x));
    
    Psy=zeros(length(t),length(x));
    [Psy(1,:),norm_x(i,1)]=wp_ini(x,sig,k,x0);
    Psy_mem_x(i,1,:)=Psy(1,1:2^(i-1):length(Psy(1,:)));
    
    a=1+1i*dt/(2*Dx(i)^2)+1i*V*dt/2;
    b_C=ones(1,length(x)-1)*(1i*dt/(4*Dx(i)^2));
    b_A=ones(1,length(x)-3)*(1i*dt/(4*Dx(i)^2));
    g=1-1i*dt/(2*Dx(i)^2)-1i*V*dt/2;

    C=sparse( full( gallery('tridiag',b_C,g,b_C)));
    C=C(2:end-1,:);
    C(1,1)=C(1,1)*2; C(end,end)=C(end,end)*2;

    A=sparse( full(gallery('tridiag',-b_A,a(2:end-1),-b_A)) );
    
    for j=1:length(t)-1
        
        D = C * transpose(Psy(j,:));
        Psy(j+1,2:end-1)= mldivide(A,D);
        
        Psy_mem_x(i,j+1,:)=Psy(j+1,1:2^(i-1):length(Psy(j,:)));
        norm_x(i,j+1)=trapeze(abs(Psy(j+1,:)).^2,x(1),x(end),length(Psy(j+1,:))-1);
    end
    
    plot(x,Psy(end,:))
    pause(0.05)
    
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

norm_t = zeros(length(Dt),length(0:dt:0.01));
Psy_mem_t = zeros(length(Dt),length(0:dt:0.01),length(x));

figure()
hold on
for i=1:length(Dt)
    
    t=0:Dt(i):0.01;
    V=zeros(1,length(x));
    
    Psy=zeros(length(t),length(x));
    [Psy(1,:),norm_t(i,1)]=wp_ini(x,sig,k,x0);
    Psy_mem_t(i,1,:)=Psy(1,:);
    
    a=1+1i*Dt(i)/(2*dx^2)+1i*V*Dt(i)/2;
    b_C=ones(1,length(x)-1)*(1i*Dt(i)/(4*dx^2));
    b_A=ones(1,length(x)-3)*(1i*Dt(i)/(4*dx^2));
    g=1-1i*Dt(i)/(2*dx^2)-1i*V*Dt(i)/2;

    C=sparse( full( gallery('tridiag',b_C,g,b_C)));
    C=C(2:end-1,:);
    C(1,1)=C(1,1)*2; C(end,end)=C(end,end)*2;

    A=sparse( full(gallery('tridiag',-b_A,a(2:end-1),-b_A)) );
    
    figure()
    for j=1:length(t)-1
        D = C * transpose(Psy(j,:));
        Psy(j+1,2:end-1)= mldivide(A,D);
        norm_t(i,j+1)=trapeze(abs(Psy(j+1,:)).^2,x(1),x(end),length(Psy(j+1,:))-1);
    end
    
    Psy_mem_t(i,:,:)=Psy(1:2^(i-1):length(Psy(:,1)),:);
    
    plot(x,Psy(end,:))
    pause(0.01)
    
    i
end
hold off

err_t=zeros(length(Dt)-1,length(0:dt:0.01),length(x));
err_t_int = zeros(length(Dt)-1,length(0:dt:0.01));

% Calcul erreur
for i=1:length(Dt)-1
    
    err_t(i,:,:)=abs(Psy_mem_t(i,:,:)-Psy_mem_t(i+1,:,:));
    
    for j = 1:length(0:Dt(1):0.01)
        
        err_t_int(i,j)=trapeze(err_t(i,j,:),0,5,length(err_t(i,j,:))-1);
        
    end
end


figure()
hold on
for i=1:length(Dt)-1
    plot(0:dt:0.01,err_t_int(i,:))
end
legend(sprintf('dt = %f',Dt(1)),sprintf('dt = %f',Dt(2)),sprintf('dt = %f',Dt(3)),sprintf('dt = %f',Dt(4)),sprintf('dt = %f',Dt(5)))
title('Erreur en fonction du temps et de la discretisation du temps')


figure()
loglog(Dt(1:end-1),err_t_int(:,end))

P=polyfit(log(Dt(1:end-1)),log(err_t_int(:,end)'),1)

