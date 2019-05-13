function I=trapeze_2D(f,a,b,c,d,n,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project          : Diffraction électronique
%% Function name    : Intégration numérique selon la méthode des trapèzes
%
%% Author           : Jean-Christophe G.-A., Thomas R., Nouffel S.
%% Date             : 5/1/2019
%
%% Description      : Intégration sur un domaine x=[a,b] et y[c,d]
%   -f      : Fonction à intégrer sous forme de matrice 2D
%   -a      : limite inférieur de x
%   -b      : limite supérieur de x
%   -c      : limite inférieur de y
%   -d      : limite supérieur de y
%   -n      : Nombre d'éléments de la discrétisation en x
%   -m      : Nombre d'éléments de la discrétisation en y
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

% Résultat
I = ( h * k / 4) * ( sum_abcd + 2*(sum_a + sum_b + sum_c + sum_d)+ 4*sum(sum(F)) );
