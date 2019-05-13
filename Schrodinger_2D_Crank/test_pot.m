clear all
clc
close all


sig=8e-11;           
lamda=1e-10;    
kp=2*(pi/lamda)*[1,0];

dx=5e-12;
dy=5e-12;
dt=1e-20;
xf=2e-9;
yf=1.4e-9;
tf=1.8e-17;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;
x0=x(end)/3;
y0=y(end)/2;

% barr1=barr_simple(x,y,1000,x(end)/2,3e-11,y(end)/2,10e-11,'Carre');

L_mur=1.5e-11;
D_fente=15e-11;
L_fente=2e-11;

barr1=barr(x,y,1000,x(end)/2,L_mur,y(end)/2,D_fente,L_fente,'Carre');
% barr2=barr(x,y,1000,x(end)/2,L_mur,y(end)/2,D_fente,L_fente,'Carre');

% figure()
% surf(x,y,barr1)

% figure()
% surf(x,y,barr2)
% colormap hot

figure()
hold on
surf(x,y,barr1)
% plot(x,barr2(1,:)/1000,'b')
% xlim([0.97 1.05]*1e-9)
legend('Barrière Carrée','Barrière Gaussienne')

% figure()
% hold on
% plot(y,barr1(:,floor(length(x)/2))/1000,'r')
% plot(y,barr2(:,floor(length(x)/2))/1000,'b')
% % xlim([0.97 1.05]*1e-9)
% legend('Barrière Carrée','Barrière Gaussienne')
