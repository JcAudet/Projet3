    function I = simpson13(f, a, b, N)
%
%  Methode de Simpson 1/3 composee
%  Programmeur: A. Lioret (21 juin 1996)
%  Reference: Analyse numerique pour ingenieurs, A. Fortin,
%             Editions de l'Ecole Polytechnique de Montreal,
%             Section 6.4.1
%
%  Prealable:
%  Si vous connaissez la fonction a integrer, vous pouvez creer
%  un fichier .m contenant cette fonction.
%  Voir par exemple le fichier fonc.m.
%
%  Exemple d'appel:
%  I = simpson13('fonc', 0, 2, 4)
%
%  Arguments
%
%  Entree:
%  1) f: Vous devez fournir soit le nom entre apostrophes (' ') du fichier .m 
%     contenant la fonction a integrer, soit les valeurs de la fonction 
%     aux points d'integration dans un vecteur ligne, correspondant avec
%     le nombre de sous-intervalles utilises (si on a N sous-intervalles,
%     on doit avoir N+1 valeurs de la fonction a integrer).
%     Exemples: - 'fonc' (correspondant a un fichier fonc.m).
%               - [f(x0) f(x1) f(x2) ... f(xN)] 
%                 ou f(x0), f(x1), f(x2), ..., f(xN) sont les 
%                 valeurs de la fonction aux points d'integration 
%                 et N est le nombre de sous-intervalles utilises
%                 (N est un multiple positif de 2).
%  2) a: Vous devez fournir la borne inferieure d'integration.
%  3) b: Vous devez fournir la borne superieure d'integration.
%  4) N: Vous devez fournir le nombre de sous-intervalles utilises.
%     Note: N doit etre un multiple positif de 2.
%  Retour:
%  1) I est la valeur de l'approximation de l'integrale definie. 
%
    if (rem(N,2) ~= 0)|(N <= 0),
      error('N doit etre un multiple positif de 2');
    end
    h = (b - a)/N;

    if (isstr(f)),
      x = a:h:b;
      fx = feval(f, x);
      if ~(all(isfinite(fx))),
        error('la fonction n''est pas definie en certains points evalues');
      end
    else
      l = length(f);
      if (l ~= N+1),
        error('le nombre de points fournis ne correspond pas avec le nombre de sous-intervalles');
      end
      fx = f;
    end
    v1 = 2:2:N;
    v2 = 3:2:N;
    yi0 = fx(1) + fx(N+1);
    yi1 = fx(v1);
    yi2 = fx(v2);

    I = h/3*(yi0 + 4*sum(yi1) + 2*sum(yi2));

