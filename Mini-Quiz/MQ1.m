clear all;
clc;
close all;

%% Questions:
%% a) Qu'est ce qu'une méthode Monte-Carlo ?
%
% La methode de Monte-Carlo désigne une famille d'algorithmes qui visent à calculer
% une valeur numérique approximative en utilisant des procédés aléatoires.

%% b) Quels sont les caractéristiques principales des problèmes pour lesquels
% il est avantageux d'appliquer une méthode Monte-Carlo ?
%
% La méthode Monte-Carlo est surtout efficace pour résoudre des problèmes ayant 
% de nombreux degrés de liberté comme des calculs d'integrales de degré
% superieur à 1 (donc des aires et des volumes).

%% c)
n=9;
Mem=zeros(100,10^5);

tic
for k=1:100
    N_sph=0;
    for j=1:10^5
        r=(rand(1,n)-0.5)*2;
        if sum(r.^2)<1
            N_sph=N_sph+1; 
        end
        Mem(k,j)=(2^n)*(N_sph/j);
    end
end
toc

V_n=(pi^(n/2))/gamma(n/2+1);

Err=abs((Mem-V_n)/V_n);

Err_avg=sum(Err,1)/100;

x=1:10^5;
Fit = polyfit(log(x),log(Err_avg),1);

figure
loglog(x,Err_avg)
hold on
loglog(x,exp(Fit(2))*x.^Fit(1))
legend('Erreur relative',sprintf('Fit de l''erreur de la forme E=A*N^p \n A=%.5f , p=%.5f',Fit(2),Fit(1)))
xlabel('Nombre de points')
ylabel('Erreur relative')


N=10^5;
ErrTheorique=1/sqrt(N);
A=Fit(2);
Ptheorique=log(ErrTheorique/A)/log(N);

fprintf('La valeur de p théorique est : %f \n',Fit(1))

%% d)
Err0=10^-3;
N=exp(log(Err0/A)/(Ptheorique));

fprintf('La valeur de N moyen pour avoir un erreur de 10^-3: \n',N)