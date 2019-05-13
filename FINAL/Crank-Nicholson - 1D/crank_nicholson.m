function P = crank_nicholson(Psy)

global A C

% a=1+1i*dt/(dx^2)+1i*dt*V/2;
% b=ones(1,length(x)-1)*(-1i*dt/(2*dx^2));
% g=1-1i*dt/(dx^2)-1i*dt*V/2;

%% Construction de C et c

% C=sparse( full( gallery('tridiag',-b,g,-b)));
% C=sparse( full( gallery('tridiag',length(x),-b,g,-b)));

% c=zeros(1,length(x));
% c(1)=-2*b*Psy(1);
% c(end)=-2*b*Psy(end);

%Le systeme a resoudre est de la forme A*(Psy_t+deltat)=b, ou b=C*(Psy_t)+c

% A=sparse( full(gallery('tridiag',b,a,b)) );
% A=sparse( full(gallery('tridiag',length(x),b,a,b)));

D = C * transpose(Psy);
P = A \ D