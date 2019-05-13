clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparaison vitesse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DT = 1e-20, DX = 2.5e-12 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes

global ly hbar m kp x0 y0 sig

hbar=1.0545718*10^-34;
m=9.10938*10^-31;

sig=8e-11;           
lamda=5e-11;    
kp=2*(pi/lamda)*[1,0];

dx=5e-12;
dy=5e-12;
dt=1e-21;
xf=1e-9;
yf=0.5e-9;
tf=1.5e-17;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/4;
y0=y(end)/2;

norm=zeros(1,length(t));
norm_th=zeros(1,length(t));
err=zeros(1,length(t));
mem1=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONDITIONS INITIALES

[psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

Psy=psy(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul des facteurs

v_mat=zeros(length(y),length(x));
V=v_mat(:);

[M,M2]=MM2(V,dx,dy,dt);

clear b c d diag diag2 diag3 f g k V v_mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init theorie

[psy_th,norm_th(1)]=analy(x,y,t(1));

err(1)=trapeze_2D(abs(squeeze(abs(psy_th))-squeeze(abs(psy))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul

vid = VideoWriter('dx=2e-12,dt=1e-20','MPEG-4');
vid.FrameRate = 150;
vid.Quality = 100;
open(vid)

gcf=figure();
colormap(gcf,'hot')
for j = 1 : length(t)-1
    tic
    b=M2*Psy;
    Psy = mldivide(M,b);
    
    psy=vec2mat(Psy,length(y))';
    norm(j+1)=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    
    [psy_th,norm_th(j+1)]=analy(x,y,t(j+1));
    
    err(j+1)=trapeze_2D(abs(squeeze(abs(psy_th))-abs(psy)),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    
    surf(x,y,abs(squeeze(abs(psy_th))-abs(psy)),'edgecolor','none');
    hold off
    view(0,90);
    daspect([1 1 1])
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    colorbar
    caxis([0 max(max(abs(psy)))]);
    
    F=getframe(gcf);
    writeVideo(vid,F);
    
    if mod(j,100)==0
        mem1=cat(3,mem1,abs(squeeze(abs(psy_th))-abs(psy)));
    end
    
    toc
end

close(vid)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% DT = 1e-20, DX = 5e-12
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Constantes
% 
% global ly hbar m kp x0 y0 sig
% 
% hbar=1.0545718*10^-34;
% m=9.10938*10^-31;
% 
% sig=8e-11;           
% lamda=5e-11;    
% kp=2*(pi/lamda)*[1,0];
% 
% dx=5e-12;
% dy=5e-12;
% dt=1e-20;
% xf=1.2e-9;
% yf=0.6e-9;
% tf=5e-17;
% 
% x=0:dx:xf;
% y=0:dy:yf; ly=length(y);
% t=0:dt:tf;
% x0=x(end)/4;
% y0=y(end)/2;
% 
% norm=zeros(1,length(t));
% norm_th=zeros(1,length(t));
% err=zeros(1,length(t));
% mem2=[];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% CONDITIONS INITIALES
% 
% [psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);
% 
% Psy=psy(:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Calcul des facteurs
% 
% v_mat=zeros(length(y),length(x));
% V=v_mat(:);
% 
% [M,M2]=MM2(V,dx,dy,dt);
% 
% clear b c d diag diag2 diag3 f g k V v_mat
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Init theorie
% 
% [psy_th,norm_th(1)]=analy(x,y,t(1));
% 
% err(1)=trapeze_2D(abs(squeeze(abs(psy_th))-squeeze(abs(psy))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Calcul
% 
% vid = VideoWriter('dx=5e-12,dt=1e-20','MPEG-4');
% vid.FrameRate = 150;
% vid.Quality = 100;
% open(vid)
% 
% gcf=figure();
% colormap(gcf,'hot')
% for j = 1 : length(t)-1
%     tic
%     b=M2*Psy;
%     Psy = mldivide(M,b);
%     
%     psy=vec2mat(Psy,length(y))';
%     norm(j+1)=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%     
%     [psy_th,norm_th(j+1)]=analy(x,y,t(j+1));
%     
%     err(j+1)=trapeze_2D(abs(squeeze(abs(psy_th))-abs(psy)),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%     
%     surf(x,y,abs(squeeze(abs(psy_th))-abs(psy)),'edgecolor','none');
%     hold off
%     view(0,90);
%     daspect([1 1 1])
%     xlim([x(1) x(end)])
%     ylim([y(1) y(end)])
%     colorbar
%     caxis([0 max(max(abs(psy)))]);
%     
%     F=getframe(gcf);
%     writeVideo(vid,F);
%     
%     if mod(j,100)==0
%         mem2=cat(3,mem2,abs(squeeze(abs(psy_th))-abs(psy)));
%     end
%     
%     toc
% end
% 
% close(vid)
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% DT = 5e-21, DX = 10e-12
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Constantes
% 
% global ly hbar m kp x0 y0 sig
% 
% hbar=1.0545718*10^-34;
% m=9.10938*10^-31;
% 
% sig=8e-11;           
% lamda=5e-11;    
% kp=2*(pi/lamda)*[1,0];
% 
% dx=10e-12;
% dy=10e-12;
% dt=1e-20;
% xf=1.2e-9;
% yf=0.6e-9;
% tf=5e-17;
% 
% x=0:dx:xf;
% y=0:dy:yf; ly=length(y);
% t=0:dt:tf;
% x0=x(end)/4;
% y0=y(end)/2;
% 
% norm=zeros(1,length(t));
% norm_th=zeros(1,length(t));
% err=zeros(1,length(t));
% mem3=[];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% CONDITIONS INITIALES
% 
% [psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);
% 
% Psy=psy(:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Calcul des facteurs
% 
% v_mat=zeros(length(y),length(x));
% V=v_mat(:);
% 
% [M,M2]=MM2(V,dx,dy,dt);
% 
% clear b c d diag diag2 diag3 f g k V v_mat
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Init theorie
% 
% [psy_th,norm_th(1)]=analy(x,y,t(1));
% 
% err(1)=trapeze_2D(abs(squeeze(abs(psy_th))-squeeze(abs(psy))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Calcul
% 
% vid = VideoWriter('dx=10e-12,dt=1e-20','MPEG-4');
% vid.FrameRate = 150;
% vid.Quality = 100;
% open(vid)
% 
% gcf=figure();
% colormap(gcf,'hot')
% for j = 1 : length(t)-1
%     tic
%     b=M2*Psy;
%     Psy = mldivide(M,b);
%     
%     psy=vec2mat(Psy,length(y))';
%     norm(j+1)=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%     
%     [psy_th,norm_th(j+1)]=analy(x,y,t(j+1));
%     
%     err(j+1)=trapeze_2D(abs(squeeze(abs(psy_th))-abs(psy)),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%     
%     surf(x,y,abs(squeeze(abs(psy_th))-abs(psy)),'edgecolor','none');
%     hold off
%     view(0,90);
%     daspect([1 1 1])
%     xlim([x(1) x(end)])
%     ylim([y(1) y(end)])
%     colorbar
%     caxis([0 max(max(abs(psy)))]);
%     
%     F=getframe(gcf);
%     writeVideo(vid,F);
%     
%     if mod(j,100)==0
%         mem3=cat(3,mem3,abs(squeeze(abs(psy_th))-abs(psy)));
%     end
%     
%     toc
% end
% 
% close(vid)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% DT = 5e-20, DX = 5e-12
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Constantes
% 
% global ly hbar m kp x0 y0 sig
% 
% hbar=1.0545718*10^-34;
% m=9.10938*10^-31;
% 
% sig=8e-11;           
% lamda=5e-11;    
% kp=2*(pi/lamda)*[1,0];
% 
% dx=5e-12;
% dy=5e-12;
% dt=5e-20;
% xf=1.2e-9;
% yf=0.6e-9;
% tf=5e-17;
% 
% x=0:dx:xf;
% y=0:dy:yf; ly=length(y);
% t=0:dt:tf;
% x0=x(end)/4;
% y0=y(end)/2;
% 
% norm=zeros(1,length(t));
% norm_th=zeros(1,length(t));
% err=zeros(1,length(t));
% mem4=[];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% CONDITIONS INITIALES
% 
% [psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);
% 
% Psy=psy(:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Calcul des facteurs
% 
% v_mat=zeros(length(y),length(x));
% V=v_mat(:);
% 
% [M,M2]=MM2(V,dx,dy,dt);
% 
% clear b c d diag diag2 diag3 f g k V v_mat
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Init theorie
% 
% [psy_th,norm_th(1)]=analy(x,y,t(1));
% 
% err(1)=trapeze_2D(abs(squeeze(abs(psy_th))-squeeze(abs(psy))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Calcul
% 
% vid = VideoWriter('dx=5e-12,dt=5e-20','MPEG-4');
% vid.FrameRate = 150;
% vid.Quality = 100;
% open(vid)
% 
% gcf=figure();
% colormap(gcf,'hot')
% for j = 1 : length(t)-1
%     tic
%     b=M2*Psy;
%     Psy = mldivide(M,b);
%     
%     psy=vec2mat(Psy,length(y))';
%     norm(j+1)=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%     
%     [psy_th,norm_th(j+1)]=analy(x,y,t(j+1));
%     
%     err(j+1)=trapeze_2D(abs(squeeze(abs(psy_th))-abs(psy)),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
%     
%     surf(x,y,abs(squeeze(abs(psy_th))-abs(psy)),'edgecolor','none');
%     hold off
%     view(0,90);
%     daspect([1 1 1])
%     xlim([x(1) x(end)])
%     ylim([y(1) y(end)])
%     colorbar
%     caxis([0 max(max(abs(psy)))]);
%     
%     F=getframe(gcf);
%     writeVideo(vid,F);
%     
%     if mod(j,100)==0
%         mem4=cat(3,mem4,abs(squeeze(abs(psy_th))-abs(psy)));
%     end
%     
%     toc
% end
% 
% close(vid)
% 
