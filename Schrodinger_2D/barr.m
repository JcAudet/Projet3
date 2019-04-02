function barr=barr(x,y,V,xm,h,yf1,yf2,a,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables

% x : Vecteur de discretisation en x
% y : Vecteur de discretisation en y
% V : Grandeur du potentiel
% h : largeur de la barriere
% xm : Position en x de la bariere 
%       (Conseiller de d'utilise la syntaxe (x(i)))
% yf1, yf2 : Position des deux fentes
%       (Conseiller de d'utilise la syntaxe (y(j)))
% a : largeur des fentes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

barr = zeros(length(y),length(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Carree

if strcmp(type,'Carre')

    X = x>(xm-h/2) & x<(xm+h/2);
    Y = y>(yf1-a/2) & y<(yf1+a/2) | y>(yf2-a/2) & y<(yf2+a/2);

    barr(~Y,X)=V;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussien

elseif strcmp(type,'Gauss')

    sig=0.01;
    X = x>(xm-h/2) & x<(xm+h/2);
    Y = y>(yf1-a/2) & y<(yf1+a/2) | y>(yf2-a/2) & y<(yf2+a/2);

    G_x = exp( -(x-xm-h/2).^2 / (2*sig^2) ) + exp( -(x-xm+h/2).^2 / (2*sig^2) );
    G_y = exp( -(y-yf1-a/2).^2 / (2*sig^2) ) + exp( -(y-yf1+a/2).^2 / (2*sig^2) )...
        + exp( -(y-yf2-a/2).^2 / (2*sig^2) ) + exp( -(y-yf2+a/2).^2 / (2*sig^2) );

    x_i=[xm-h/2 xm+h/2];
    y_j=[yf1-a/2 yf1+a/2 yf2-a/2 yf2+a/2];
    for i=1:2
        for j=1:4
            gxy = ( exp( -((y-y_j(j)).^2./(2.*sig.^2)) ) )'* exp(-((x-x_i(i)).^2)./(2.*sig.^2));
            barr = barr + V .* gxy / max(max(gxy));
        end
    end
    
    barr(~Y,:) = ones(size(barr(~Y,:),1),1) * V * G_x / max(G_x);
    
    barr(:,X) = ( ones(size(barr(:,X),2),1) * V * G_y / max(G_y) )';

    barr(~Y,X)=V;

end
