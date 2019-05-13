function P = run_kutt_4(dt, dx, psy, pot)

k1= dt * eq_schro(dx, psy , pot);
k2= dt * eq_schro(dx, psy + k1/2  , pot);
k3= dt * eq_schro(dx, psy + k2/2 , pot);
k4= dt * eq_schro(dx, psy + k3 , pot);

P = psy + (k1 + 2*k2 + 2*k3 + k4)/6;