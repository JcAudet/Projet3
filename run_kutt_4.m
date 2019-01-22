function P = run_kutt_4(x , psy , pot)

global dt

% % Normalisation
% norm=sum(abs(psy));
% psy=psy./norm;

k1= dt * eq_schro(x , psy , pot);
k2= dt * eq_schro(x , psy + k1/2  , pot);
k3= dt * eq_schro(x , psy + k2/2 , pot);
k4= dt * eq_schro(x, psy + k3 , pot);

P = psy + (k1 + 2*k2 + 2*k3 + k4)/6;