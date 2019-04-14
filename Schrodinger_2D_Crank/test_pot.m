clear all
clc
close all


sig=8e-11;           
lamda=1e-10;    
kp=2*(pi/lamda)*[1,0];

dx=0.94e-11;
dy=0.94e-11;
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

barr2=barr(x,y,1000,x(end)/2,3e-11,y(end)/2,10e-11,10e-11,'Carre');

% figure()
% surf(x,y,barr1)
% colormap hot

figure()
surf(x,y,barr2)
colormap hot
