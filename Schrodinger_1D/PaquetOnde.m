clear all;
clc;
close all;
global dx dt; 
format long

x=linspace(0,10,1001);
dx=x(2)-x(1);
t=linspace(0,1000,3000);
dt=t(2)-t(1);
    
% Paquet d'onde

sig=0.3;
k=10;
x_0=5;
psy=exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);

% Normalisation

norm=trapeze(abs(psy),x(1),x(end),length(psy)-1);
psy=psy./norm;

% Déplacement temporel

a=1;
m=9.109*10^(-31);
h_barre=6.62607004*10^(-34)/(2*pi);

for i=1:length(t)-1
    
    y=sqrt(1+4*h_barre^2*t(i+1)^2/(m^2*a^4));
    
    for j=1:1:1001
        %psy=(2/(pi*a^2))^(1/4) * y^(-1/2) * exp(i*k0*(x(j)-h_barre*k0*t/(2*m)))*exp(-2/(a*y)^2*(x(j)-h_barre*k0*t(i+1)/m)^2);
        psy_norm(j)= sqrt(2/(pi*a^2)) * 1/y * exp(-2/(a*y)^2*(x(j)-x_0-h_barre*k*t(i+1)/m)^2);
    end
    
    plot(x,psy_norm)
    norm=trapeze(abs(psy),x(1),x(end),length(psy)-1)
    pause(0.001)
end








