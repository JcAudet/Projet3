function Psy_0 = Psy_0(x)
deltaX=0.01;
G=exp(-((x).^2)./0.09);
produitscalaire=(deltaX*G*transpose(G))^0.5;
A=1/produitscalaire;


Psy_0=A*(exp(-((x).^2)./0.02)).*(exp(complex(0,1).*20.*x));


end