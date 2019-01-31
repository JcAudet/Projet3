clear all;
clc;
close all;

n=10;
Mem=zeros(100,10^5);

tic
for k=1:100
    N_sph=0;
    for j=1:10^5
        r=(rand(1,n)-0.5)*2;
        if sum(r.^2)<1
            N_sph=N_sph+1; 
        end
        Mem(k,j)=(2^n)*(N_sph/j);
    end
    clc;
    k
end
toc

V_n=(pi^(n/2))/gamma(n/2+1);

Err=abs((Mem-V_n)/V_n);

Err_avg=sum(Err,1)/100;

x=1:10^5;
Fit = polyfit(log(x),log(Err_avg),1);

figure
loglog(x,Err_avg)
% cftool(x,Err_avg)
hold on
loglog(x,exp(Fit(2))*x.^Fit(1))

% Err_0=0.002;
% N=1/(Err_0^2);
% 
% tic
% for i=1:N
%     r=(rand(1,n)-0.5)*2;
%     if sum(r.^2)<1
%        N_sph=N_sph+1; 
%     end
% end
% toc