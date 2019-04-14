
clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

global ly hbar m kp x0 y0 sig

% hbar=6.62607*10^-34;
% m=9.10938*10^-31;
% 
% sig=8e-11;           
% lamda=1e-10;    
% kp=2*(pi/lamda)*[1,0];
% 
% dx=0.4e-11;
% dy=0.4e-11;
% dt=8e-21;
% xf=1e-9;
% yf=0.5e-9;
% tf=5e-18;
% 
% x=0:dx:xf;
% y=0:dy:yf; ly=length(y);
% t=0:dt:tf;
% x0=4*x(end)/10;
% y0=y(end)/2;


hbar=1.0545718*10^-34;
m=9.10938*10^-31;

sig=8e-11;           
lamda=1e-10;    
kp=2*(pi/lamda)*[1,0];

dx=0.8e-11;
dy=0.8e-11;
dt=1e-20;
xf=2e-9;
yf=1.4e-9;
tf=5e-18;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/3;
y0=y(end)/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Theorie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dx

Dx=[dx dx*2 dx*4];
x1=0:Dx(end):xf;

% Init matrices
norm_x=zeros(length(Dx),length(t));
norm_th_x=zeros(length(Dx),length(t));
err_x=zeros(length(Dx),length(t));

for w=1:length(Dx)

    x=0:Dx(w):xf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONDITIONS INITIALES
    
    [psy,norm_x(w,1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

    Psy=[];
    Psy(:,1)=psy(:);

    Psy_mat=zeros(length(y),length(x),length(t));
    Psy_mat(:,:,1)=psy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat=zeros(length(y),length(x));
    % v_mat=( barr(x,y,1000,x(floor(length(x)/2)),4e-11,y(floor(length(y)/2)),15e-11,10e-11,'Carre') );
    % v_mat=sparse( ( barr_simple(x,y,1000,x(length(x)/2),1e-14,y(length(y)/2),8e-14,'Carre') )' );
    V=v_mat(:);

    b = 1 + 1i*hbar*dt*(1/Dx(w)^2 + 1/dy^2)/(2*m) + 1i*dt*V/(2*hbar);
    f = 1 - 1i*hbar*dt*(1/Dx(w)^2 + 1/dy^2)/(2*m) - 1i*dt*V/(2*hbar);

    c=-1i*hbar*dt*(1/Dx(w)^2)/(4*m);
    d=-1i*hbar*dt*(1/dy^2)/(4*m);

    g=1i*hbar*dt*(1/Dx(w)^2)/(4*m);
    k=1i*hbar*dt*(1/dy^2)/(4*m);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Creation de M et M2, v et v2

    diag=ones(1,length(V)-1)*d;
    diag1=ones(1,length(V)-1)*d;
    diag(ly:ly:length(V)-1)=0;
    diag1(ly:ly:length(V)-1)=0;

    diag2=ones(1,length(V)-1)*k;
    diag3=ones(1,length(V)-1)*k;
    diag2(ly:ly:length(V)-1)=0;
    diag3(ly:ly:length(V)-1)=0;

    M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c diag diag1]);
    M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g diag2 diag3]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Init theorie
    
    Psy_x=zeros(1,length(x));
    Psy_y=zeros(1,length(y));
    Psy_th=zeros(length(t),length(y),length(x));
    
    [Psy_th(1,:,:),norm_th_x(w,1)]=analy(x,y,t(1));
    
    norm_th_x(w,1)=trapeze_2D(abs(squeeze(Psy_th(1,:,:))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    err_x(w,1)=trapeze_2D(abs(squeeze(abs(Psy_th(1,:,:)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

%     figure()
    tic
    for j = 1 : length(t)-1
        
        b=M2*Psy(:,j);
        Psy(:,j+1) = mldivide(M,b);
        
        psy=vec2mat(Psy(:,j+1),length(y))';
        Psy_mat(:,:,j+1)=psy;
        
        [Psy_th(j+1,:,:),norm_th_x(w,j+1)]=analy(x,y,t(j+1));
        
        norm_x(w,j+1)=trapeze_2D(abs(squeeze(Psy_mat(:,:,j+1))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        err_x(w,j+1)=trapeze_2D(abs(squeeze(abs(Psy_th(j+1,:,:)))-squeeze(abs(Psy_mat(:,:,j+1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
              
    end
    toc
    w    
end

figure()
hold on
for i=1:length(Dx)
    plot(t,norm_x(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de la norme de la resolution en fonction de DX')

figure()
hold on
for i=1:length(Dx)
    plot(t,norm_th_x(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de la norme de la solution analyique en fonction de DX')

figure()
hold on
for i=1:length(Dx)
    plot(t,err_x(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de l''erreur en fonction de DX')

P_x=polyfit(log(Dx),log(err_x(:,end)'),1)

figure()
plot(log(Dx),log(err_x(:,end)))
hold on
plot(log(Dx),P_x(1)*log(Dx)+P_x(2))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dt

dt=5e-19;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
x0=4*x(end)/10;
y0=y(end)/2;

Dt=[dt dt*2 dt*4];
t1=0:Dt(end):tf;

norm_t=zeros(length(Dt),length(t1));
norm_th_t=zeros(length(Dt),length(t1));
err_t=zeros(length(Dt),length(t1));


for w=1:length(Dt)
    
    t=0:Dt(w):tf;
    
    norm=zeros(1,length(t));
    norm_th=zeros(1,length(t));
    err=zeros(1,length(t));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONDITIONS INITIALES
    
    [psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

    Psy=[];
    Psy(:,1)=psy(:);

    Psy_mat=zeros(length(y),length(x),length(t));
    Psy_mat(:,:,1)=psy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat=zeros(length(y),length(x));
    % v_mat=( barr(x,y,1000,x(floor(length(x)/2)),4e-11,y(floor(length(y)/2)),15e-11,10e-11,'Carre') );
    % v_mat=sparse( ( barr_simple(x,y,1000,x(length(x)/2),1e-14,y(length(y)/2),8e-14,'Carre') )' );
    V=v_mat(:);

    b = 1 + 1i*hbar*Dt(w)*(1/dx^2 + 1/dy^2)/(2*m) + 1i*Dt(w)*V/(2*hbar);
    f = 1 - 1i*hbar*Dt(w)*(1/dx^2 + 1/dy^2)/(2*m) - 1i*Dt(w)*V/(2*hbar);

    c=-1i*hbar*Dt(w)*(1/dx^2)/(4*m);
    d=-1i*hbar*Dt(w)*(1/dy^2)/(4*m);

    g=1i*hbar*Dt(w)*(1/dx^2)/(4*m);
    k=1i*hbar*Dt(w)*(1/dy^2)/(4*m);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Creation de M et M2, v et v2

    diag=ones(1,length(V)-1)*d;
    diag1=ones(1,length(V)-1)*d;
    diag(ly:ly:length(V)-1)=0;
    diag1(ly:ly:length(V)-1)=0;

    diag2=ones(1,length(V)-1)*k;
    diag3=ones(1,length(V)-1)*k;
    diag2(ly:ly:length(V)-1)=0;
    diag3(ly:ly:length(V)-1)=0;

    M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c diag diag1]);
    M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g diag2 diag3]);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Init theorie

    Psy_th=zeros(length(t),length(y),length(x));
    
    [Psy_th(1,:,:),norm_th(1)]=analy(x,y,t(1));
    
    err(1)=trapeze_2D(abs(squeeze(abs(Psy_th(1,:,:)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

    tic
%     figure()
    for j = 1 : length(t)-1
        b=M2*Psy(:,j);
        Psy(:,j+1) = mldivide(M,b);

        Psy_mat(:,:,j+1)=vec2mat(Psy(:,j+1),length(y))';
        norm(j+1)=trapeze_2D(abs(Psy_mat(:,:,j+1)).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        
        [Psy_th(j+1,:,:),norm_th(j+1)]=analy(x,y,t(j+1));
        
        err(j+1)=trapeze_2D(abs(squeeze(abs(Psy_th(j+1,:,:)))-abs(Psy_mat(:,:,j+1))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

%         subplot(2,2,1)
%         surf(x,y,abs(Psy_mat(:,:,j)).^2,'edgecolor','none');
%         hold off
%         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm(j)))
%         view(0,90);
%         daspect([1 1 1])    
% 
%         subplot(2,2,2)
%         plot(t(1:j),err(1:j))
% 
%         subplot(2,2,3)
%         plot(x,abs(Psy_mat(floor(length(y)/2),:,j)).^2)
%         hold on
%         plot(x,abs(squeeze(Psy_th(j,floor(length(y)/2),:))).^2)
%         legend('CN','ANALY')
% 
%         hold off
% 
%         subplot(2,2,4)
%         plot(t(1:j),norm(1:j))
%         hold on
%         plot(t(1:j),norm_th(1:j))
%         legend('CN','ANALY')
% 
% 
%         hold off
%         pause(0.00001)
        
        
    end
    
    norm_t(w,:)=norm(1:2^(length(Dt)-w):length(norm)); 
    norm_th_t(w,:)=norm_th(1:2^(length(Dt)-w):length(norm_th)); 
    err_t(w,:)=err(1:2^(length(Dt)-w):length(err));
    
    toc
end

figure()
hold on
for i=1:length(Dt)
    plot(t,norm_t(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dt(1)),sprintf('dx = dy = %d',Dt(2)),sprintf('dx = %d',Dt(3)))
title('Evolution de la norme de la resolution en fonction de DT')

figure()
hold on
for i=1:length(Dt)
    plot(t,norm_th_t(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dt(1)),sprintf('dx = dy = %d',Dt(2)),sprintf('dx = %d',Dt(3)))
title('Evolution de la norme de la solution analytique en fonction de DT')


figure()
hold on
for i=1:length(Dt)
    plot(t,err_t(i,:))
end
legend(sprintf('dx = dy = %d',Dt(1)),sprintf('dx = dy = %d',Dt(2)),sprintf('dx = %d',Dt(3)))
title('Evolution de l''erreur en fonction de DT')

P_t=polyfit(log(Dt),log(err_t(:,end))',1)

figure()
plot(log(Dt),log(err_t(:,end)))
hold on
plot(log(Dt),P_t(1)*log(Dt)+P_t(2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lam

Lam=[lamda lamda*2 lamda*4 lamda*8];

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;

% Init matrices
norm_l=zeros(length(Lam),length(t));
norm_th_l=zeros(length(Lam),length(t));
err_l=zeros(length(Lam),length(t));

% M et M2
v_mat=zeros(length(y),length(x));
V=v_mat(:);

b = 1 + 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) + 1i*dt*V/(2*hbar);
f = 1 - 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) - 1i*dt*V/(2*hbar);

c=-1i*hbar*dt*(1/dx^2)/(4*m);
d=-1i*hbar*dt*(1/dy^2)/(4*m);

g=1i*hbar*dt*(1/dx^2)/(4*m);
k=1i*hbar*dt*(1/dy^2)/(4*m);

diag=ones(1,length(V)-1)*d;
diag1=ones(1,length(V)-1)*d;
diag(ly:ly:length(V)-1)=0;
diag1(ly:ly:length(V)-1)=0;

diag2=ones(1,length(V)-1)*k;
diag3=ones(1,length(V)-1)*k;
diag2(ly:ly:length(V)-1)=0;
diag3(ly:ly:length(V)-1)=0;

M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
    [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
    ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c diag diag1]);
M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
    [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
    ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g diag2 diag3]);

for w=1:length(Lam)

    kp=2*(pi/Lam(w))*[1,0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONDITIONS INITIALES
    
    [psy,norm_l(w,1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

    Psy=[];
    Psy(:,1)=psy(:);

    Psy_mat=zeros(length(y),length(x),length(t));
    Psy_mat(:,:,1)=psy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat=zeros(length(y),length(x));
    % v_mat=( barr(x,y,1000,x(floor(length(x)/2)),4e-11,y(floor(length(y)/2)),15e-11,10e-11,'Carre') );
    % v_mat=sparse( ( barr_simple(x,y,1000,x(length(x)/2),1e-14,y(length(y)/2),8e-14,'Carre') )' );
    V=v_mat(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Init theorie
    
    Psy_th=zeros(length(t),length(y),length(x));
    
    [Psy_th(1,:,:),norm_th_l(w,1)]=analy(x,y,t(1));
    
    norm_th_l(w,1)=trapeze_2D(abs(squeeze(Psy_th(1,:,:))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    err_l(w,1)=trapeze_2D(abs(squeeze(abs(Psy_th(1,:,:)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

%     figure()
    tic
    for j = 1 : length(t)-1
        
        b=M2*Psy(:,j);
        Psy(:,j+1) = mldivide(M,b);
        
        psy=vec2mat(Psy(:,j+1),length(y))';
        Psy_mat(:,:,j+1)=psy;
        
        [Psy_th(j+1,:,:),norm_th_l(w,j+1)]=analy(x,y,t(j+1));
        
        norm_l(w,j+1)=trapeze_2D(abs(squeeze(Psy_mat(:,:,j+1))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        err_l(w,j+1)=trapeze_2D(abs(squeeze(abs(Psy_th(j+1,:,:)))-squeeze(abs(Psy_mat(:,:,j+1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        
%         subplot(1,2,1)
%         surf(x,y,abs(squeeze(Psy_mat(j,:,:))).^2,'edgecolor','none');
%         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_x(w,j)))
%         view(0,90);
%         subplot(1,2,2)
%         surf(x,y,abs(squeeze(Psy_th(j,:,:))).^2,'edgecolor','none');
%         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_th_x(w,j)))
%         view(0,90);
%         hold off
%         
%         pause(0.01)
        
    end
    toc
    w    
end

figure()
hold on
for i=1:length(Dx)
    plot(t,norm_l(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de la norme de la resolution en fonction de la longueur d''onde')

figure()
hold on
for i=1:length(Dx)
    plot(t,norm_th_l(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de la norme de la solution analyique en fonction de la longueur d''onde')

figure()
hold on
for i=1:length(Lam)
    plot(t,err_l(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Lam(1)),sprintf('dx = dy = %d',Lam(2)),sprintf('dx = %d',Lam(3)),sprintf('dx = %d',Lam(4)))
title('Evolution de l''erreur avec differentes valeurs de longueur d''onde')

P_l=polyfit(log(Lam),log(err_l(:,end)'),1)

figure()
plot(log(Lam),log(err_l(:,end)))
hold on
plot(log(Lam),P_l(1)*log(Lam)+P_l(2))

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Sigma
% 
% Sig=[sig/2 sig sig*2];
% 
% x=0:dx:xf;
% y=0:dy:yf; ly=length(y);
% t=0:dt:tf;
% 
% % Init matrices
% norm_s=zeros(length(Sig),length(t));
% norm_th_s=zeros(length(Sig),length(t));
% err_s=zeros(length(Sig),length(t));
% 
% 
% for w=1:length(Sig)
% 
%     sig=Sig(w);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% CONDITIONS INITIALES
%     
%     [psy,norm_s(w,1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);
% 
%     Psy=[];
%     Psy(:,1)=psy(:);
% 
%     Psy_mat=zeros(length(y),length(x),length(t));
%     Psy_mat(:,:,1)=psy;
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Calcul des facteurs
% 
%     
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Init theorie
%     
%     Psy_th=zeros(length(t),length(y),length(x));
%     
%     [Psy_th(1,:,:),norm_th_s(w,1)]=analy(x,y,t(1));
%     
%     norm_th_s(w,1)=trapeze_2D(abs(squeeze(Psy_th(1,:,:))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%     err_s(w,1)=trapeze_2D(abs(squeeze(abs(Psy_th(1,:,:)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Calcul
% 
% %     figure()
%     tic
%     for j = 1 : length(t)-1
%         
%         b=M2*Psy(:,j);
%         Psy(:,j+1) = mldivide(M,b);
%         
%         psy=vec2mat(Psy(:,j+1),length(y))';
%         Psy_mat(:,:,j+1)=psy;
%         
%         [Psy_th(j+1,:,:),norm_th_s(w,j+1)]=analy(x,y,t(j+1));
%         
%         norm_s(w,j+1)=trapeze_2D(abs(squeeze(Psy_mat(:,:,j+1))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%         err_s(w,j+1)=trapeze_2D(abs(squeeze(abs(Psy_th(j+1,:,:)))-squeeze(abs(Psy_mat(:,:,j+1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%         
% %         subplot(1,2,1)
% %         surf(x,y,abs(squeeze(Psy_mat(j,:,:))).^2,'edgecolor','none');
% %         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_x(w,j)))
% %         view(0,90);
% %         subplot(1,2,2)
% %         surf(x,y,abs(squeeze(Psy_th(j,:,:))).^2,'edgecolor','none');
% %         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_th_x(w,j)))
% %         view(0,90);
% %         hold off
% %         
% %         pause(0.01)
%         
%     end
%     toc
%     w    
% end
% 
% figure()
% hold on
% for i=1:length(Sig)
%     plot(t,norm_x(i,:))
% end
% hold off
% legend(sprintf('Sigma = %d',Sig(1)),sprintf('Sigma = %d',Sig(2)),sprintf('Sigma = %d',Sig(3)))
% title('Evolution de la norme de la resolution en fonction de l''ecart type')
% 
% figure()
% hold on
% for i=1:length(Sig)
%     plot(t,norm_th_x(i,:))
% end
% hold off
% legend(sprintf('Sigma = %d',Sig(1)),sprintf('Sigma = %d',Sig(2)),sprintf('Sigma = %d',Sig(3)))
% title('Evolution de la norme de la solution analyique en fonction de l''ecart type')
% 
% figure()
% hold on
% for i=1:length(Sig)
%     plot(t,err_s(i,:))
% end
% hold off
% legend(sprintf('Sigma = %d',Lam(1)),sprintf('Sigma = %d',Lam(2)),sprintf('Sigma = %d',Lam(3)),sprintf('dx = %d',Lam(4)))
% title('Evolution de l''erreur avec differentes valeurs de longueur d''onde')
% 
% P_s=polyfit(log(Sig),log(err_s(:,end)'),1)
% 
% figure()
% plot(log(Sig),log(err_s(:,end)))
% hold on
% plot(log(Sig),P_s(1)*log(Sig)+P_s(2))




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% DY
% 
% Dy=[dy dy*2 dy*4];
% y1=0:Dy(end):yf;
% 
% x=0:dx:xf;
% y=0:dy:yf; ly=length(y);
% t=0:dt:tf;
% 
% % Init matrices
% norm_y=zeros(length(Dy),length(t));
% norm_th_y=zeros(length(Dy),length(t));
% err_y=zeros(length(Dy),length(t));
% 
% for w=1:length(Dy)
% 
%     y=0:Dy(w):yf;
%     
%     ly=length(y);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% CONDITIONS INITIALES
%     
%     [psy,norm_y(w,1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);
% 
%     Psy=[];
%     Psy(:,1)=psy(:);
% 
%     Psy_mat=zeros(length(y),length(x),length(t));
%     Psy_mat(:,:,1)=psy;
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Calcul des facteurs
% 
%     v_mat=zeros(length(y),length(x));
%     % v_mat=( barr(x,y,1000,x(floor(length(x)/2)),4e-11,y(floor(length(y)/2)),15e-11,10e-11,'Carre') );
%     % v_mat=sparse( ( barr_simple(x,y,1000,x(length(x)/2),1e-14,y(length(y)/2),8e-14,'Carre') )' );
%     V=v_mat(:);
% 
%     b = 1 + 1i*hbar*dt*(1/dx^2 + 1/Dy(w)^2)/(2*m) + 1i*dt*V/(2*hbar);
%     f = 1 - 1i*hbar*dt*(1/dx^2 + 1/Dy(w)^2)/(2*m) - 1i*dt*V/(2*hbar);
% 
%     c=-1i*hbar*dt*(1/dx^2)/(4*m);
%     d=-1i*hbar*dt*(1/Dy(w)^2)/(4*m);
% 
%     g=1i*hbar*dt*(1/dx^2)/(4*m);
%     k=1i*hbar*dt*(1/Dy(w)^2)/(4*m);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Creation de M et M2, v et v2
% 
%     %M=sparse(ones(length(b_vect),1),ones(length(b_vect),1),b_vect);
%     M=diag(b);
% 
%     %M2=sparse(ones(length(b_vect),1),ones(length(b_vect),1),f_vect);
%     M2=diag(f);
% 
%     tic
%     for j = 1 : length(y)
%         for i = 1 : length(x)
% 
%             ind=indexeur(j,i);
% 
%             if i > 1
%                 M(ind,indexeur(j,i-1)) = c;
%                 M2(ind,indexeur(j,i-1)) = g;
%             end
% 
%             if i < length(x)
%                 M(ind,indexeur(j,i+1)) = c;
%                 M2(ind,indexeur(j,i+1)) = g;
%             end
% 
%             if j < length(y)
%                 M(ind,indexeur(j+1,i)) = d;
%                 M2(ind,indexeur(j+1,i)) = k;
%             end
% 
%             if j > 1
%                 M(ind,indexeur(j-1,i)) = d;
%                 M2(ind,indexeur(j-1,i)) = k;
%             end
%         end
%     end
%     toc
% 
%     tic
%     M=sparse(M);
%     M2=sparse(M2);
%     toc
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Init theorie
%     
%     Psy_x=zeros(1,length(x));
%     Psy_y=zeros(1,length(y));
%     Psy_th=zeros(length(t),length(y),length(x));
%     
%     [Psy_th(1,:,:),norm_th_y(w,1)]=analy(x,y,t(1));
%     
%     norm_th_y(w,1)=trapeze_2D(abs(squeeze(Psy_th(1,:,:))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%     err_y(w,1)=trapeze_2D(abs(squeeze(abs(Psy_th(1,:,:)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Calcul
% 
% %     figure()
%     tic
%     for j = 1 : length(t)-1
%         
%         b=M2*Psy(:,j);
%         Psy(:,j+1) = mldivide(M,b);
%         
%         psy=vec2mat(Psy(:,j+1),length(y))';
%         Psy_mat(:,:,j+1)=psy;
%         
%         [Psy_th(j+1,:,:),norm_th_y(w,j+1)]=analy(x,y,t(j+1));
%         
%         norm_y(w,j+1)=trapeze_2D(abs(squeeze(Psy_mat(:,:,j+1))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%         err_y(w,j+1)=trapeze_2D(abs(squeeze(abs(Psy_th(j+1,:,:)))-squeeze(abs(Psy_mat(:,:,j+1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%         
% %         subplot(1,2,1)
% %         surf(x,y,abs(squeeze(Psy_mat(j,:,:))).^2,'edgecolor','none');
% %         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_x(w,j)))
% %         view(0,90);
% %         subplot(1,2,2)
% %         surf(x,y,abs(squeeze(Psy_th(j,:,:))).^2,'edgecolor','none');
% %         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_th_x(w,j)))
% %         view(0,90);
% %         hold off
% %         
% %         pause(0.01)
%         
%     end
%     toc  
%     w    
% end
% 
% figure()
% hold on
% for i=1:length(Dy)
%     plot(t,norm_x(i,:))
% end
% hold off
% legend(sprintf('dy = %d',Dy(1)),sprintf('dy = %d',Dy(2)),sprintf('dy = %d',Dy(3)))
% title('Evolution de la norme de la resolution en fonction de DY')
% 
% figure()
% hold on
% for i=1:length(Dy)
%     plot(t,norm_th_y(i,:))
% end
% hold off
% legend(sprintf('dy = %d',Dy(1)),sprintf('dy = %d',Dy(2)),sprintf('dy = %d',Dy(3)))
% title('Evolution de la norme de la solution analyique en fonction de DY')
% 
% figure()
% hold on
% for i=1:length(Dy)
%     loglog(t,err_y(i,:))
% end
% hold off
% legend(sprintf('dy = %d',Dy(1)),sprintf('dy = %d',Dy(2)),sprintf('dy = %d',Dy(3)))
% title('Evolution de l''erreur en fonction de DY')
% 
% P_y=polyfit(log(Dy),log(err_y(:,end)'),1)
% 
% figure()
% plot(log(Dy),log(err_y(:,end)))
% hold on
% plot(log(Dy),P_y(1)*log(Dy)+P_y(2))


