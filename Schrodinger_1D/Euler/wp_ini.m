function [psy,nor]=wp_ini(x,sig,k,x_0)

psy=exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);
psy(1)=0;
psy(end)=0;

psy_abs=abs(psy);
nor=trapeze(psy_abs,x(1),x(end),length(psy_abs)-1);
psy=psy./nor;

nor=trapeze(psy_abs,x(1),x(end),length(psy_abs)-1);