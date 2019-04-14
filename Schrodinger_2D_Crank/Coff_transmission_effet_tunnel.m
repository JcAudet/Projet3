clear all
clc

I=1000;
V0=linspace(0,2*10^(-13),I);
T=zeros(1,I);
R=zeros(1,I);

for i=1:1:I
V=V0(i);                    %Potentiel de la barrière 
E = 8.0108831e-14;  %Energie de l'électron, 50KeV
a=10^(-12);                  %Largeur de la barrière en m
m=9.10938*10^-31;           %Masse de l'électron
hbar=6.62607*10^-34*(1/(2*pi)); 
lambda=a*sqrt(2*m*V/hbar^2);
epsilon=E/V;
    if epsilon > 1
     T1=(1+(1/(4*epsilon*(epsilon-1)))*(sin(lambda*sqrt(epsilon-1)))^2)^(-1);
     R1=(1+(4*epsilon*(epsilon-1))/(sin(lambda*sqrt(epsilon-1)))^2)^(-1);
    else 
     T1=(1+1/(4*epsilon*(1-epsilon))*(sinh(lambda*sqrt(1-epsilon)))^2)^(-1);
     R1=(T1/(4*epsilon*(1-epsilon))*(sinh(lambda*sqrt(1-epsilon)))^2);
    end
T(i)=T1;
R(i)=R1;
end

sum(T)+sum(R)

figure
hold on
plot(V0,T,'blue')
plot(V0,R,'red')
legend('Coeff de transmission','Coeff de réflexion')
xlabel('Potentiel (eV)')
ylabel('Coefficients')
title('Coefficients de réflexion et transmission pour E=50KeV et a=10^(-12) m')



