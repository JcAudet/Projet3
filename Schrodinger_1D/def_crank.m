function [C,A] = def_crank(x)

global dt dx

g=1-D*dt/(dx^2);
a=1i*dt/(2*dx^2);
b=1i*dt/(dx^2)+1i*dt*V/2;

%% Construction de C et c
C=sparse( full(gallery('tridiag',-b(2:end),g,-b(1:end-1))) );

% c=zeros(1,length(x));
% c(1)=-2*b*Psy(1);
% c(end)=-2*b*Psy(end);

%Le systeme a resoudre est de la forme A*(Psy_t+deltat)=b, ou b=C*(Psy_t)+c

A=sparse( full(gallery('tridiag',length(x),b,a,b)) );



