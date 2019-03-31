function [psy,nor]=wp_ini(x,sig,k,x_0)

psy=exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);
psy(1)=0;
psy(end)=0;

nor=trapeze(abs(psy),x(1),x(end),length(abs(psy))-1);
psy=psy./nor;

nor=trapeze(abs(psy),x(1),x(end),length(abs(psy))-1);