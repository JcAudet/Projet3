function [C,A] = def_crank(x)

global dt dx

D=1i;
% D=complex(0,1)*(6*10^-34)/(4*pi*9.109*10^-31);

g=1-D*dt/(dx^2);
a=1+D*dt/(dx^2);
b=-0.5*D*dt/(dx^2);

%% Construction de C et c
C=full(gallery('tridiag',length(x),-b,g,-b));

% c=zeros(1,length(x));
% c(1)=-2*b*Psy(1);
% c(end)=-2*b*Psy(end);

%Le systeme a resoudre est de la forme A*(Psy_t+deltat)=b, ou b=C*(Psy_t)+c

A=full(gallery('tridiag',length(x),b,a,b));



