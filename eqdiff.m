function f = eqdif(t,y)
f = zeros(size(y));
% Syst�me de 2 �quations diff�rentielles
%y'_1(t) = y_2(t) -y_1(t);
%y'_2(t) = t + y_1(t);
% y(1) = y_1(t)
% y(2) = y_2(t)

f(1)= i.*y./2;
