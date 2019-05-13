function barr=barr(x,y,V,xm,h,ym,d,a,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project          : Diffraction électronique
%% Function name    : Analyse de convergence de Crank-Nicholson 2D
%
%% Author           : Jean-Christophe G.-A., Thomas R., Nouffel S.
%% Date             : 5/1/2019
%
%% Description      : Obtention de l'ordre spatial et temporel de 
%   - x     : Vecteur de discretisation en x
%   - y     : Vecteur de discretisation en y
%   - V     : Grandeur du potentiel
%   - xm    : Position en x de la bariere 
%            (Conseiller de d'utilise la syntaxe (x(i)))
%   - h     : largeur de la barriere
%   - ym    : Position milieu entre les deux fentes
%            (Conseiller de d'utilise la syntaxe (y(length(y)/2))
%   - a     : largeur des fentes
%   - type  : Type de frontiere du potentiel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

barr = zeros(length(y),length(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Carree

if strcmp(type,'Carre')

    X = x>(xm-h/2) & x<(xm+h/2);
    Y = y>(ym+d/2-a/2) & y<(ym+d/2+a/2) | y>(ym-d/2-a/2) & y<(ym-d/2+a/2);

    barr(~Y,X)=V;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussien

elseif strcmp(type,'Gauss')

%   Parametres
    sig_x=h/10;
    sig_y=a/10;
    
%   Conditions
    X = x>(xm-h/2) & x<(xm+h/2);
    Y = y>(ym-d/2-a/2) & y<(ym-d/2+a/2) | y>(ym+d/2-a/2) & y<(ym+d/2+a/2);

%   Creation des potentiel frontiere en x et y
    G_x = exp( -(x-xm-h/2).^2 / (2*sig_x^2) ) + exp( -(x-xm+h/2).^2 / (2*sig_x^2) );
    G_y = exp( -(y-ym-d/2-a/2).^2 / (2*sig_y^2) ) + exp( -(y-ym-d/2+a/2).^2 / (2*sig_y^2) )...
        + exp( -(y-ym+d/2-a/2).^2 / (2*sig_y^2) ) + exp( -(y-ym+d/2+a/2).^2 / (2*sig_y^2) );

%   Creation des coins frontiere
    x_i=[xm-h/2 xm+h/2];
    y_j=[ym-d/2-a/2 ym-d/2+a/2 ym+d/2-a/2 ym+d/2+a/2];
    for i=1:2
        for j=1:4
            gxy = ( exp( -((y-y_j(j)).^2./(2.*sig_y.^2)) ) )'* exp(-((x-x_i(i)).^2)./(2.*sig_x.^2));
            barr = barr + V .* gxy / max(max(gxy));
        end
    end
    
    barr(~Y,:) = ones(size(barr(~Y,:),1),1) * V * G_x / max(G_x);
    
    barr(:,X) = ( ones(size(barr(:,X),2),1) * V * G_y / max(G_y) )';

    barr(~Y,X)=V;

end
