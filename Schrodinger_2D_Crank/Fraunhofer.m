clc;clear all;
I_0=1;

%dimensions
dx=0.6*10^-11;dy=dx;
L=217*dy;
H=-(0.5-2/3)*301*dx;
thetamax=atan(L/(2*H));
theta=-thetamax:0.01:thetamax;
% v_mat=barr(x,y,1000,x(end)/2,3e-11,y(end)/2,20e-11,15e-11,'Carre');

%a : distance qui separe les fentes
%b : epaisseur des fentes
lamda=5*10^-11;
b=20-11;a=15e-11;
alpha=(pi*a/lamda).*sin(theta);
beta=(pi*b/lamda).*sin(theta);


I=I_0.*(cos(alpha).*sin(beta)./beta).^2;


plot(H*tan(theta),I);

hold on
title('Intensite sur une tranche');
hold off