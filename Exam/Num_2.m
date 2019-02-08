clear all;
clc;

format long

dx=0.1;
ddx=0.01;

%% Derivee premiere
dy=(exp(-2*dx)-8*exp(-dx)+8*exp(dx)-exp(2*dx))/(12*dx);
dy_2=(exp(-2*ddx)-8*exp(-ddx)+8*exp(ddx)-exp(2*ddx))/(12*ddx);

ddy=(-exp(-2*dx)+2*exp(-dx)-2*exp(dx)+exp(2*dx))/(2*dx^3);
ddy_2=(-exp(-2*ddx)+2*exp(-ddx)-2*exp(ddx)+exp(2*ddx))/(2*ddx^3);

e1=1-dy
e2=1-dy_2
e3=1-ddy
e4=1-ddy_2

Rap_e1=abs(1-dy)/abs(1-dy_2)
Ordre1=log10(Rap_e1)
Rap_e2=abs(1-ddy)/abs(1-ddy_2)
Ordre2=log10(Rap_e2)
