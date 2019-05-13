function [Psy,norm]=analy_1D(x,t)

global hbar m k sig x_0


tet_x=atan(2*hbar*t/(m*sig^2));
delx=(sig/2)*sqrt(1+(4*hbar^2*t^2)/(m^2*sig^4));

Psy=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet_x/2).*...
    exp(1i*k*(x-k*hbar*t/(2*m))).*...
    exp(-(x-x_0-hbar*k*t/m).^2 / (sig^2+2*1i*hbar*t/m));

norm=trapeze(abs(Psy).^2,x(1),x(end),length(x)-1);