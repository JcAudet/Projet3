clear all;
clc;

y=[5.67 7.28 9.14 13.53 17.73 22.07 28.04];
x=[0.5 1 1.5 2 2.5 3 3.5];

x_2=x.^2;

mx_2y=mean(x_2.*y);
mx_2=mean(x_2);
my=mean(y);
mx_4=mean(x_2.^2);
mx_2_2=(mean(x_2))^2;

A=(mx_2y-mx_2*my)/(mx_4-mx_2_2);
B=my-A*mx_2;

figure
plot(x,A*x.^2+B,'r')
hold on
plot(x,y)