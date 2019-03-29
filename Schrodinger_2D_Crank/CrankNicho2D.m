clear all;clc;

global dx dy dt

dx=0.5*10^-13;
dy=0.5*10^-13;
dt=10^-24;
Imax=50; %indices d'espaces max
Jmax=50;
Kmax=100; %indice de temps max
%centre en (x,y)=(5,5)
x_0=25*dx;y_0=25*dx;
%proprietse initiales
sig=10^-13;
%si accelerer a partir de repos a 50KeV
lamda=5.4847*10^-12;
k=2*(pi/lamda)*[1,1];


hbar=6.62607*10^-34;
m=9.10938*10^-31;


%CONDITIONS INITIALES
%LA condition initiale est une fonciton Psy(x,y). Pour representer une
%fonction de 2 variables dans un vecteur, il faut adopter une convention.
%Voir doc Overleaf ou note de cour sur matrice 2D...
Psy0=[];
for j = 1 : Jmax
    for i = 1 : Imax 
        %La boncle if sert a placer les Conditions Frontieres
        if i==1
            Psy0=[Psy0 (Psy_haut(j)+CI(i*dx,j*dy,x_0,y_0,sig,k))/2];
        elseif i==Imax
            Psy0=[Psy0 (Psy_bas(j)+CI(i*dx,j*dy,x_0,y_0,sig,k))/2];
        elseif j==1;
            Psy0=[Psy0 (Psy_gauche(i)+CI(i*dx,j*dy,x_0,y_0,sig,k))/2];
        elseif j== Jmax
            Psy0=[Psy0 (Psy_droite(i)+CI(i*dx,j*dy,x_0,y_0,sig,k))/2];
        else
    Psy0=[Psy0 CI(i*dx,j*dy,x_0,y_0,sig,k)];
        end
    end
end
%Psy0 n'est pas normalisé, il faut le faire (DOMAINE EN 2D ATTENTION)
A=dx*dy*norm(Psy0)^0.5;
Psy0=Psy0/A;
%Psy_0 en vecteur colonne
Psy0=transpose(Psy0);



b_vect=[];
f_vect=[];
%stocker b et f dans vecteur OU MATRICES??
for j = 1 : Jmax
    for i = 1 : Imax 
        var=(1+1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m)+1i*dt*V_potentiel(i*dx,j*dy)/(2*hbar));
       b_vect=[b_vect var];
         var2=(1-1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m)-1i*dt*V_potentiel(i*dx,j*dy)/(2*hbar));
       f_vect=[f_vect var2];
    end
end
c=-1i*hbar*dt*(1/dx^2)/(4*m);
a=c;
d=-1i*hbar*dt*(1/dy^2)/(4*m);
e=d;

g=1i*hbar*dt*(1/dx^2)/(4*m);
h=g;
k=1i*hbar*dt*(1/dy^2)/(4*m);
p=k;


%%Ecrire Les 2 matrices M et M2 et les deux vecteur v et v2 dans la meme boucle

%M=sparse(ones(length(b_vect),1),ones(length(b_vect),1),b_vect);
 M=diag(b_vect);
v=zeros(length(b_vect),1);

%M2=sparse(ones(length(b_vect),1),ones(length(b_vect),1),f_vect);
 M2=diag(f_vect);
v2=zeros(length(f_vect),1);

for j = 1 : Jmax
    for i = 1 : Imax 

        ind=indexeur(i,j,Imax);
        if i > 1
            M(ind,indexeur(i-1,j,Imax)) = c;
            M2(ind,indexeur(i-1,j,Imax)) = g;
        else
            v(ind)=v(ind)+c*Psy_haut(j);
            v2(ind)=v2(ind)+g*Psy_haut(j);

        end
        if i < Imax 
            M(ind,indexeur(i+1,j,Imax)) = a;
            M2(ind,indexeur(i+1,j,Imax)) = h;

        else
            v(ind)=v(ind)+a*Psy_bas(j);
            v2(ind)=v2(ind)+h*Psy_bas(j);
        end
        if j < Jmax 
            M(ind,indexeur(i,j+1,Imax)) = d;
            M2(ind,indexeur(i,j+1,Imax)) = k;
        else
            v(ind)=v(ind)+d*Psy_droite(i);
            v2(ind)=v2(ind)+k*Psy_droite(i);
        end
        if j > 1
            M(ind,indexeur(i,j-1,Imax)) = e;
            M2(ind,indexeur(i,j-1,Imax)) = p;
        else
            v(ind)=v(ind)+e*Psy_gauche(i);
            v2(ind)=v2(ind)+p*Psy_gauche(i);
        end
        
       
        end
    end

M=sparse(M);
M2=sparse(M2);
%UTILISER MATRICE CREUSE SPARSE???

%On stock tous les veteur dans MatricePsy
MatricePsy=[];
MatricePsy(:,1)=Psy0;
tic
for k = 2 : Kmax
%Re-ecrire en Ax=b
A=M;;
b=M2*MatricePsy(:,k-1)+v2-v;
x = mldivide(A,b);
%x=A\b;
MatricePsy(:,k)=x;
k
end
toc

%Re-ecrire les vecteur comme des matrices pour ensuite les ploter
cell=[];
s=0;
for k= 1:Kmax
    for j=1 : Jmax
        for i=1 : Imax
            s=s+1;
    cell(i,j,k)= MatricePsy(s,k);
    
        end
    end
    s=0;
end
% 
for k=1 : Kmax
    W=real(cell(:,:,k));
    surf(W);
    pause(0.02);
end

