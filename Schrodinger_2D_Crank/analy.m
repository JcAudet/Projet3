function [Psy,norm]=analy(x,y,t)

global hbar m kp sig x0 y0


tet_x=atan(2*hbar*t/(m*sig^2));
delx=(sig/2)*sqrt(1+(4*hbar^2*t^2)/(m^2*sig^4));

Psy_x=(1/sqrt(sqrt(2*pi)*delx))*exp(-1i*tet_x/2).*...
    exp(1i*kp(1)*(x-kp(1)*hbar*t/(2*m))).*...
    exp(-(x-x0-hbar*kp(1)*t/m).^2 / (sig^2+2*1i*hbar*t/m));

tet_y=atan(2*hbar*t/(m*sig^2));
dely=(sig/2)*sqrt(1+(4*hbar^2*t^2)/(m^2*sig^4));

Psy_y=(1/sqrt(sqrt(2*pi)*dely))*exp(-1i*tet_y/2).*...
    exp(1i*kp(2)*(y-kp(2)*hbar*t/(2*m))).*...
    exp(-(y-y0-hbar*kp(2)*t/m).^2 / (sig^2+2*1i*hbar*t/m));

Psy(:,:)=Psy_y'*Psy_x;

norm=trapeze_2D(abs(Psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);