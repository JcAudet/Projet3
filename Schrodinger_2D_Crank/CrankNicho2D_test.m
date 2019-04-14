clear all;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretisation

global ly hbar m kp sig x0 y0

hbar=1.0545718*10^-34;
m=9.10938*10^-31;

sig=8e-11;           
lamda=5e-11;    
kp=2*(pi/lamda)*[1,0];

dx=0.6e-11;
dy=0.6e-11;
dt=1e-20;
xf=2e-9;
yf=1.4e-9;
tf=1e-16;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/3;
y0=y(end)/2;

C=dt*m/(hbar*norm(kp)^2*dx*dy)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONDITIONS INITIALES
% La condition initiale est une fonciton Psy(x,y). Pour representer une 
% fonction de 2 variables dans un vecteur, il faut adopter une convention.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONDITIONS INITIALES

norm=zeros(1,length(t));
norm_th=zeros(1,length(t));
err=zeros(1,length(t));
[psy,norm(1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

Psy=[];
Psy(:,1)=psy(:);
Psy_mat=psy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul des facteurs

v_mat=barr(x,y,1000,x(end)/2,3e-11,y(end)/2,10e-11,5e-11,'Carre');
V=v_mat(:);

b = 1 + 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) + 1i*dt*V/(2*hbar);
f = 1 - 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) - 1i*dt*V/(2*hbar);

c=-1i*hbar*dt*(1/dx^2)/(4*m);
d=-1i*hbar*dt*(1/dy^2)/(4*m);

g=1i*hbar*dt*(1/dx^2)/(4*m);
k=1i*hbar*dt*(1/dy^2)/(4*m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creation de M et M2, v et v2

tic

M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
    [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
    ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c ones(1,length(V)-1)*d ones(1,length(V)-1)*d]);
M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
    [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
    ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g ones(1,length(V)-1)*k ones(1,length(V)-1)*k]);
% M2=sparse(1:length(V),1:length(V),f,length(V),length(V),2*length(V));
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul

[psy,norm_th(1)]=analy(x,y,t(1));

err(1)=trapeze_2D(abs(squeeze(abs(Psy_mat(:,:,1)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

figure('units','normalized','outerposition',[0 0 1 1])
tic
for  j = 1 : length(t)-1

    b=M2*Psy(:,j);
    Psy(:,j+1) = mldivide(M,b);
    
    Psy_mat(:,:,j+1)=vec2mat(Psy(:,j+1),length(y))';

    norm(j+1)=trapeze_2D(abs(Psy_mat(:,:,j+1)).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    [psy,norm_th(j+1)]=analy(x,y,t(j+1));

    err(j)=trapeze_2D(abs(squeeze(abs(Psy_mat(:,:,j+1)))-squeeze(abs(psy))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

%     subplot(2,2,1)
    surf(x,y,abs(Psy_mat(:,:,j)).^2,'edgecolor','none');
    hold off
    title(sprintf('Temps = %e  Norme: %.10f',t(j),norm(j)))
    view(0,90);
    daspect([1 1 1])
% 
%     subplot(2,2,2)
%     plot(t(1:j),err(1:j))
% 
%     subplot(2,2,3)
%     plot(x,abs(Psy_mat(floor(length(y)/2),:,j)).^2)
%     hold on
%     plot(x,abs(psy(floor(length(y)/2),:)).^2)
%     legend('CN','ANALY')
% 
%     hold off
% 
%     subplot(2,2,4)
%     plot(t(1:j),norm(1:j))
%     hold on
%     plot(t(1:j),norm_th(1:j))
%     legend('CN','ANALY')
% 
% 
    hold off
    pause(0.00001)

end

toc

figure()
plot(t,norm)
hold on
plot(t,norm_th)
legend('CN','ANALY')

figure()
surf(x,y,abs(Psy_mat(:,:,1)).^2,'edgecolor','none');
hold on
surf(x,y,abs(Psy_mat(:,:,end)).^2,'edgecolor','none');
hold off
title(sprintf('Temps = %e  Norme: %.10f',t(j),norm(j)))
view(0,90);
daspect([1 1 1])



