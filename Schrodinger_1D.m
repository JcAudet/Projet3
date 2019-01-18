clear all;
clc;

x=linspace(0,10,1000);
t=linspace(0,0.03,300);

sig=0.1;
k=50;

P_0_u=exp(-((x-1).^2)./(2.*sig.^2)).*exp(1i.*k.*x);

Norm=P_0_u.*conj(P_0_u);
A=sum(Norm);
P_0=P_0_u./A;
Lap=zeros(1000,length(P_0));
for i=1:1000
    Lap(2:end-1)=(P_0(3:length(P_0))+P_0(1:length(P_0)-2)-...
        2.*P_0(2:length(P_0)-1))./(x(2)-x(1));
end


figure
plot(x,P_0);

