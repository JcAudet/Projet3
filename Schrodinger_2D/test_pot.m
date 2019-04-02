clear all
clc
close all


x=linspace(0,10,10000);
y=linspace(0,10,10000);

barr=sparse( barr(x,y,1000,5,0.2,4,6,0.2,'Carre') );

figure()
imagesc(x,y,barr)
colormap hot

