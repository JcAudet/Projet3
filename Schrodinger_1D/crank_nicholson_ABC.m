function P = crank_nicholson_ABC(Psy)

global A C C1 C2 C3 C4

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
  
    A(1,1)=C1;      % Conditions frontières absorbantes
    A(1,2)=C2;      % Valeurs à changer dans les matrices à résoudre
    A(end,end)=C1;
    A(end-1,end)=C2;
    
    C(1,1)=C3;
    C(1,2)=C4;
    C(end,end)=C3;
    C(end-1,end)=C4;

b= C * transpose(Psy);%+c; 

P = A \ b;