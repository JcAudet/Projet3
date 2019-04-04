clear all
clc
close all


dx=10^-14;
dy=10^-14;
dt=10^-25;

Imax=150;                % Indices d'espaces max
Jmax=150;
Kmax=1000;               % Indice de temps max

x=0:dx:(Imax-1)*dx;
y=0:dy:(Jmax-1)*dy;

barr=sparse( barr(x,y,1000,x(length(x)/2),1e-13,y(length(y)/2),2e-13,1e-13,'Gauss') );

figure()
surf(x,y,barr)
colormap hot

