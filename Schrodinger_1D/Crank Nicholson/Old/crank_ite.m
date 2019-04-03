function P = crank_ite(x , Psy)

global A C

b=C*transpose(Psy);% + transpose(c);
P = A\b;