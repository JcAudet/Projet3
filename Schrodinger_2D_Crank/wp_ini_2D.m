function [psy,nor]=wp_ini_2D(x,y,sig_x,sig_y,kp,x_0,y_0)

psy = (exp(-(((x-x_0).^2)./(2.*sig_x.^2))).*exp(1i.*kp(1).*x))'*...
    (exp(-((y-y_0).^2)./(2.*sig_y.^2)).*exp(1i.*kp(2).*y));
psy(1,:)=0;psy(end,:)=0;psy(:,1)=0;psy(:,end)=0;

nor=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
psy=psy./sqrt(nor);
nor=trapeze_2D(abs(psy).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
