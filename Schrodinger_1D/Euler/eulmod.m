    function   [y] = eulmod(x , psy , pot)
%
global dt
       fy0 = eq_schro(x, psy , pot);
       yint = psy + dt * fy0; 
       fyint = eq_schro(x,yint,pot);
       y = psy +  dt/2*(fy0 + fyint);
    end
    