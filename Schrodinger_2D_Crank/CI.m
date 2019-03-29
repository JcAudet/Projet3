function  f = CI(x,y,x_0,y_0,sig,k)

f=exp(-((x-x_0).^2 + (y-y_0)^2)./(2.*sig.^2)).*exp(1i.*(k(1).*x + k(2).*y));

end