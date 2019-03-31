    function  [y] = euler(x , psy , pot)
global dt  
    k1 = eq_schro(x, psy , pot);
    psy = psy + dt * k1;
    y= psy;
    end
    