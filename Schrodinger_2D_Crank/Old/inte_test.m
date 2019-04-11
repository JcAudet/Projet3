clear all
clc
close all



x=linspace(0,1,1000);
y=linspace(0,1,1000);

f=x'*y^2;

norm=trapeze_2D(f,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);