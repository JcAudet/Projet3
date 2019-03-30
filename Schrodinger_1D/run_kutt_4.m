function P = run_kutt_4(t, psy , pot)

dt=t(2)-t(1);

k1= dt * eq_schro(psy , pot);
k2= dt * eq_schro(psy + k1/2  , pot);
k3= dt * eq_schro(psy + k2/2 , pot);
k4= dt * eq_schro(psy + k3 , pot);

P = psy + (k1 + 2*k2 + 2*k3 + k4)/6;