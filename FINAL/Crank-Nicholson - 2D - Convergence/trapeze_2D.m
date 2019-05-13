function I=trapeze_2D(f,a,b,c,d,n,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project          : Diffraction �lectronique
%% Function name    : Int�gration num�rique selon la m�thode des trap�zes
%
%% Author           : Jean-Christophe G.-A., Thomas R., Nouffel S.
%% Date             : 5/1/2019
%
%% Description      : Int�gration sur un domaine x=[a,b] et y[c,d]
%   -f      : Fonction � int�grer sous forme de matrice 2D
%   -a      : limite inf�rieur de x
%   -b      : limite sup�rieur de x
%   -c      : limite inf�rieur de y
%   -d      : limite sup�rieur de y
%   -n      : Nombre d'�l�ments de la discr�tisation en x
%   -m      : Nombre d'�l�ments de la discr�tisation en y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Verifications
if (n <= 0)
  error('n doit etre positif');
end

l = size(f);
if (l(2) ~= n+1 || l(1)~=m+1)
error('le nombre de points fournis ne correspond pas avec le nombre de sous-intervalles');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul
h = (b - a) / n;
k = (d - c) / m;

sum_abcd = f(1,1) + f(1,n+1) + f(m+1,1) + f(m+1,n+1);
sum_a = sum(f(2:m,1));
sum_b = sum(f(2:m,n+1));
sum_c = sum(f(1,2:n));
sum_d = sum(f(m+1,2:n));
F = f(2:m,2:n);

% R�sultat
I = ( h * k / 4) * ( sum_abcd + 2*(sum_a + sum_b + sum_c + sum_d)+ 4*sum(sum(F)) );
