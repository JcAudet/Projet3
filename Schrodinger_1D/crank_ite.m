function P = crank_ite(x , Psy)

b=C*transpose(Psy) + transpose(c);
P = A\b;