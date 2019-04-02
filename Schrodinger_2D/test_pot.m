clear all
clc
close all


x=linspace(0,10,1000);
y=linspace(0,10,1000);

barr=barr(x,y,1000,5,0.5,4,6,1);

figure()
hold on
surf(x,y,barr)