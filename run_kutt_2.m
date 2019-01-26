function P = run_kutt_2(x , psy , pot)

global dt

% % Normalisation
% norm=sum(abs(psy));
% psy=psy./norm;

k1= dt * eq_schro(x , psy , pot);
k2= dt * eq_schro(x , psy + k1/2  , pot);


P = psy + k2;