function Psy_0 = Psy_0(x)

global dx

G=exp(-((x).^2)./0.02).*(exp(1i.*20.*x));
produitscalaire=dx*sum(abs(G));
A=1/produitscalaire;

Psy_0=A*(exp(-((x).^2)./0.02)).*(exp(complex(0,1).*200.*x));


end