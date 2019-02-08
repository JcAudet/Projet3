function P = crank_nicholson(x , Psy)

global C A; 

% c=

b=C*transpose(Psy);% + transpose(c);
P = A\b;