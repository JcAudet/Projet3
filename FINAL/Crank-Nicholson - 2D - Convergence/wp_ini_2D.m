function [psy,nor]=wp_ini_2D(x,y,sig_x,sig_y,kp,x_0,y_0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project          : Diffraction électronique
%% Function name    : Initialisation du paquet d'onde
%
%% Author           : Jean-Christophe G.-A., Thomas R., Nouffel S.
%% Date             : 5/1/2019
%
%% Description      : Initialise un paquet d'onde selon les paramètre:
%   - x,y       : Domaine
%   - sig_x,y   : Écart-type
%   - kp        : vecteur d'onde
%   - x_0,y_0   : Position Initiale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Création de l'onde dans le domaine
psy = (exp(-((y-y_0).^2)./(sig_y.^2)).*exp(1i.*kp(2).*y))'*...
    (exp(-(((x-x_0).^2)./(sig_x.^2))).*exp(1i.*kp(1).*x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Condition frontière nulles
psy(1,:)=0;psy(end,:)=0;psy(:,1)=0;psy(:,end)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalisation
nor=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
psy=psy./sqrt(nor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalisation sortante
nor=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
