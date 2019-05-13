function [psy,nor]=wp_ini(x,sig,k,x_0)

% psy=(1/sqrt(sqrt(pi)*sig)).*exp(-((x-x_0).^2)./(2.*sig.^2)).*exp(1i.*k.*x);
psy=exp(-((x-x_0).^2)./(sig.^2)).*exp(1i.*k.*x);
psy(1)=0;
psy(end)=0;

nor=trapeze(abs(psy).^2,x(1),x(end),length(psy)-1);
psy=psy./sqrt(nor);
nor=trapeze(abs(psy).^2,x(1),x(end),length(psy)-1);