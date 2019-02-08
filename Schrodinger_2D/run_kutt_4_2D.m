function P = run_kutt_4_2D(psy, pot)

global dt

% % Normalisation
% norm=sum(abs(psy));
% psy=psy./norm;

k1= dt * eq_schro_2D(psy , pot);
k2= dt * eq_schro_2D(psy + k1/2  , pot);
k3= dt * eq_schro_2D(psy + k2/2 , pot);
k4= dt * eq_schro_2D(psy + k3 , pot);

P = psy + (k1 + 2*k2 + 2*k3 + k4)/6;