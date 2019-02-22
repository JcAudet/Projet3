% Quizz 3
clc
clear all

Ta=-10;
Ti=20;
N=101;
L=0.3;
k=1;
h=3;
C=1000;
p=2000;
q=2000;

%Discrétisation
dx=0.3/(N-1);
dt=100*dx^2*C*p/k;
i=(0:1:(N-1));
x0=0;
t0=0;

x=x0+i*dx;
t=t0+i*dt;

%Constantes frontières
C1=k;
C2=-h;
C3=h*Ta;
D1=k;
D2=h;
D3=-h*Ti;

% Martice A initiale
A=diag(-2*ones(1,N),0)+diag(ones(1,N-1),1)+diag(ones(1,N-1),-1);
A(1,1)=2*C2*dx-3*C1;
A(1,2)=4*C1;
A(1,3)=-C1;
A(N,N)=3*D1+2*D2*dx;
A(N,N-1)=-4*D1;
A(N,N-2)=D1;

% Matrice M
M=diag(ones(1,N),0);
M(1,1)=0;
M(N,N)=0;

% Température initiale
for j=1:N
    T(j,1)=Ta+(Ti-Ta)/(2*k/(h*L)+1)*(k/(h*L)+x(j)/L);
end


% S(x)
for j=1:N
    b(j,1)=q/((1+((x(j)-L/2)/dx)^2)*k);
end
alpha=C*p/k;


% Algorithme a appliquer aux points du domaine
T_n=T;
Z=0.9;
C=(M-(Z*dt/(alpha*dx^2))*A);

for j=1:100
    plot(x,T_n)
    pause (0.1)   
       
    D=(M+(1-Z)*(dt/(alpha*dx^2))*A)*T_n-dt*((Z*b)+(1-Z)*b)/alpha;
    Tn_1=C\D;
    T_n=Tn_1;  
    
end



%Teq=maxT(0)+0.99*(Tmax-maxT(0))


