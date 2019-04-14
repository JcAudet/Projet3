clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

global ly hbar m kp sig_x sig_y x0 y0

hbar=1.0545718*10^-34;
m=9.10938*10^-31;

sig_x=8e-11;                 % Largeur initiale du packet d'onde
sig_y=8e-11;
% lamda=6e-12;
lamda=2e-10;    % Accelerer du repos avec 50KeV
kp=2*(pi/lamda)*[1,0];

dx=0.63e-11;
dy=1e-11;
dt=1e-21;
xf=2e-9;
yf=1e-9;
tf=1e-18;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/3;
y0=y(end)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dx

Dx=[dx*8 dx*4 dx*2 dx];

x1=0:Dx(1):xf;

norm_x=zeros(length(Dx),length(t));
Psy_mem_x=zeros(length(Dx),length(t),length(y),length(x1));

for w=1:length(Dx)

    x=0:Dx(w):xf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONDITIONS INITIALES
    
    [psy,norm_x(w,1)] = wp_ini_2D(x,y,sig_x,sig_y,kp,x(end)/3,y(end)/2);

    Psy=[];
    Psy(:,1)=psy(:);
    Psy_mem_x(w,1,:,:)=psy(:,1:2^(w-1):length(x));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat=zeros(length(y),length(x));
    V=v_mat(:);

    b = 1 + 1i*hbar*dt*(1/Dx(w)^2 + 1/dy^2)/(2*m) + 1i*dt*V/(2*hbar);
    f = 1 - 1i*hbar*dt*(1/Dx(w)^2 + 1/dy^2)/(2*m) - 1i*dt*V/(2*hbar);

    c=-1i*hbar*dt*(1/Dx(w)^2)/(4*m);
    d=-1i*hbar*dt*(1/dy^2)/(4*m);

    g=1i*hbar*dt*(1/Dx(w)^2)/(4*m);
    k=1i*hbar*dt*(1/dy^2)/(4*m);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Creation de M et M2, v et v2

    M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c ones(1,length(V)-1)*d ones(1,length(V)-1)*d]);
    M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
        [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
        ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g ones(1,length(V)-1)*k ones(1,length(V)-1)*k]);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

%     figure()
    tic
    for k = 1 : length(t)-1
        
        b=M2*Psy(:,k);
        Psy(:,k+1) = mldivide(M,b);

        psy=vec2mat(Psy(:,k+1),length(y))';
        
        Psy_mem_x(w,k+1,:,:)=psy(:,1:2^(w-1):length(x));
        
        norm_x(w,k+1)=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

        surf(x,y,abs(psy).^2,'edgecolor','none');
        hold off
        title(sprintf('Temps = %e  Norme: %.10f',t(k),norm_x(w,k)))
        view(0,90);
        
        pause(0.01)
        
    end
    toc
    w    
end

err_x = zeros(length(Dx)-1,length(t),length(y),length(x1));
err_x_int = zeros(length(Dx)-1,length(t));

% Calcul erreur
for i=1:length(Dx)-1
    
    err_x(i,:,:,:)=abs( Psy_mem_x(i,:,:,:)-Psy_mem_x(i+1,:,:,:) );
    
    for j = 1:length(t)
        err_x_int(i,j)=trapeze_2D(squeeze(err_x(i,j,:,:)),x1(1),x1(end),y(1),y(end),length(x1)-1,length(y)-1);
    end
    
end

figure()
hold on
for i=1:length(Dx)
    plot(t,norm_x(i,:))
end
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)),sprintf('dx = %d',Dx(4)))
title('Evolution de la norme en fonction de DX')
ylabel('Norme')
xlabel('Time')

figure()
hold on
for i=1:length(Dx)-1
    plot(t,err_x_int(i,:))
end
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)),sprintf('dx = %d',Dx(4)))
title('Evolution de l''erreur en fonction de DX')
ylabel('Erreur')
xlabel('Time')

P_x=polyfit(log(Dx(1:end-1)),log(err_x_int(:,end)'),1)

figure()
loglog(Dx(1:end-1),err_x_int(:,end))
hold on
loglog(Dx(1:end-1),P_x(1)*Dx(1:end-1)+P_x(2))

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% DY
% 
% dx=0.56e-11;
% dy=0.56e-11;
% dt=5e-21;
% 
% x=0:dx:1e-9;
% y1=0:dy*8:1e-9;
% t=0:dt:1e-17;
% 
% Dy=[dy dy*2 dy*4 dy*8];
% norm_y=zeros(length(Dy),length(t));
% Psy_mem_y=zeros(length(Dy),length(t),length(y1),length(x));
% 
% for w=1:length(Dy)
% 
%     
%     y=0:Dy(w):1e-9;
%     
%     ly=length(y);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% CONDITIONS INITIALES
%     
%     [psy,norm_y(w,1)] = wp_ini_2D(x,y,sig_x,sig_y,kp,x(end)/3,y(end)/2);
% 
%     Psy=[];
%     Psy(:,1)=psy(:);
%     Psy_mem_y(w,1,:,:)=psy(1:2^(length(Dy)-w):length(y),:);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Calcul des facteurs
% 
%     v_mat=zeros(length(y),length(x));
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
%            if i > 1
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
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Calcul
%     
%     tic
%     for k = 1 : length(t)-1
%         b=M2*Psy(:,k);
%         Psy(:,k+1) = mldivide(M,b);
% 
%         psy=vec2mat(Psy(:,k+1),length(y))';
%         Psy_mem_y(w,k+1,:,:)=psy(1:2^(length(Dy)-w):length(y),:);
%         
%         norm_y(w,k+1)=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%         
%     end
%     toc
% end
% 
% figure()
% hold on
% for i=1:length(Dy)
%     surf(x,y,abs(squeeze(Psy_mem_y(i,end,:,:))),'edgecolor','none')
%     view(0,90)
% end
% 
% 
% err_y=zeros(length(Dy)-1,length(t),length(y1),length(x));
% err_y_int = zeros(length(Dy)-1,length(t));
% 
% % Calcul erreur
% for i=1:length(Dy)-1
%     
%     err_y(i,:,:,:)=abs(Psy_mem_y(i,:,:,:)-Psy_mem_y(i+1,:,:,:));
%     
%     for j = 1:length(t)
%         
%         err_y_int(i,j)=trapeze_2D(squeeze(err_y(i,j,:,:)),x(1),x(end),y1(1),y1(end),length(x)-1,length(y1)-1);
%         
%     end
% end
% 
% 
% figure()
% hold on
% for i=1:length(Dy)-1
%     plot(t,err_y_int(i,:))
% end
% legend(sprintf('dx = dy = %d',Dy(1)),sprintf('dx = dy = %d',Dy(2)),sprintf('dx = %d',Dy(3)))
% title('Erreur en fonction du temps et de la discretisation de l''espace')
% 
% %% LOG10
% P_y=polyfit(log10(Dy(1:end-1)),log10(err_y_int(:,end)'),1)
% 
% figure()
% loglog(Dy(1:end-1),err_y_int(:,end))
% hold on
% loglog(Dy(1:end-1),P_y(1)*Dy(1:end-1)+P_y(2))
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Dt
% 
% dx=0.56e-11;
% dy=0.56e-11;
% dt=5e-21;
% 
% x=0:dx:1e-9;
% y=0:dy:1e-9;
% t1=0:dt*8:1e-17;
% 
% ly=length(y);
% 
% Dt=[dt dt*2 dt*4 dt*8];
% 
% Psy_mem_t=zeros(length(Dt),length(t1),length(y),length(x));
% 
% for w=1:length(Dt)
%     
%     t=0:Dt(w):1e-17;
%     
%     norm_t=zeros(1,length(t));
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% CONDITIONS INITIALES
%     
%     [psy,norm_t(1)] = wp_ini_2D(x,y,sig_x,sig_y,kp,x(end)/3,y(end)/2);
% 
%     Psy=[];
%     Psy(:,1)=psy(:);
%     Psy_mat=zeros(length(t),length(y),length(x));
%     Psy_mat(1,:,:)=psy(:,:);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Calcul des facteurs
% 
%     v_mat=zeros(length(y),length(x));
%     V=v_mat(:);
% 
%     b = 1 + 1i*hbar*Dt(w)*(1/dx^2 + 1/dy^2)/(2*m) + 1i*Dt(w)*V/(2*hbar);
%     f = 1 - 1i*hbar*Dt(w)*(1/dx^2 + 1/dy^2)/(2*m) - 1i*Dt(w)*V/(2*hbar);
% 
%     c=-1i*hbar*Dt(w)*(1/dx^2)/(4*m);
%     d=-1i*hbar*Dt(w)*(1/dy^2)/(4*m);
% 
%     g=1i*hbar*Dt(w)*(1/dx^2)/(4*m);
%     k=1i*hbar*Dt(w)*(1/dy^2)/(4*m);
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
%            if i > 1
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
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Calcul
%     
%     tic
%     for k = 1 : length(t)-1
%         b=M2*Psy(:,k);
%         Psy(:,k+1) = mldivide(M,b);
% 
%         psy=vec2mat(Psy(:,k+1),length(y))';
%         Psy_mat(k+1,:,:)=psy;
%         
%         norm_t(w,k+1)=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%         
%     end
%     
%     Psy_mem_t(w,:,:,:)=Psy_mat(1:2^(length(Dt)-w):length(t),:,:);
%     
%     toc
% end
% 
% figure()
% hold on
% for i=1:length(Dt)
%     surf(x,y,abs(squeeze(Psy_mem_t(i,end,:,:))),'edgecolor','none')
%     view(0,90)
% end
% 
% 
% err_t=zeros(length(Dt)-1,length(t1),length(y),length(x));
% err_t_int = zeros(length(Dt)-1,length(t1));
% 
% % Calcul erreur
% for i=1:length(Dt)-1
%     
%     err_t(i,:,:,:)=abs(Psy_mem_t(i,:,:,:)-Psy_mem_t(i+1,:,:,:));
%     
%     for j = 1:length(t)
%         
%         err_t_int(i,j)=trapeze_2D(squeeze(err_t(i,j,:,:)),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%         
%     end
% end
% 
% 
% figure()
% hold on
% for i=1:length(Dt)-1
%     plot(t,err_t_int(i,:))
% end
% legend(sprintf('dx = dy = %d',Dt(1)),sprintf('dx = dy = %d',Dt(2)),sprintf('dx = %d',Dt(3)))
% title('Erreur en fonction du temps et de la discretisation de l''espace')
% 
% P_t=polyfit(log10(Dt(1:end-1)),log(err_t_int10(:,end)'),1)
% 
% figure()
% loglog(Dt(1:end-1),err_t_int(:,end))
% hold on
% loglog(Dt(1:end-1),P_t(1)*Dt(1:end-1)+P_t(2))


