    function I = trapeze(f, a, b, n)
%
%  Methode des trapezes composee
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
%  I = trapeze('fonc', 0, 2, 5)
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
%                 et N est le nombre de sous-intervalles utilises.
%  2) a: Vous devez fournir la borne inferieure d'integration.
%  3) b: Vous devez fournir la borne superieure d'integration.
%  4) N: Vous devez fournir le nombre de sous-intervalles utilises.
%  Retour:
%  1) I est la valeur de l'approximation de l'integrale definie. 
%
    if (n <= 0),
      error('n doit etre positif');
    end
    h = (b - a)/n;

    if (isstr(f)),
      x = a:h:b;
      fx = feval(f, x);
      if ~(all(isfinite(fx))),
        error('la fonction n''est pas definie pour certains points evalues');
      end
    else
      l = length(f);
      if (l ~= n+1),
        error('le nombre de points fournis ne correspond pas avec le nombre de sous-intervalles');
      end
      fx = f;
    end
    yi0 = fx(1) + fx(n+1);
    yi1 = fx(2:n);
 
    I = h/2*(yi0 + 2*sum(yi1));
