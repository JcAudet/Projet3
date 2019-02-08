clear all;
clc;
close all;

format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes physique
k=20;
m=0.05;
g=9.8;
y_0=-0.005;
v_0=-0.025;
alp=0.16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation vecteur temps et espace
exp=-4.75:0.25:-3;
T=(2*pi)*(sqrt(m/k));
dt=(10.^exp).*T./pi;
d_max=zeros(3,length(dt));

%% Resolution interative
for j=1:length(dt)
    t=0:dt(j):0.6*T;
    Y=zeros(1,length(t));
    Y(1)=y_0;
    Y(2)=Y(1)*(1-(k*dt(j)^2)/(2*m))+v_0*dt(j)*(1-(alp*dt(j))/(2*m))-(g*dt(j)^2)/2;
    for i=2:length(t)-1
        Y(i+1)=(Y(i)*(2*m/dt(j)^2-k)+Y(i-1)*(alp/(2*dt(j))-m/dt(j)^2)-m*g)/(m/dt(j)^2+alp/(2*dt(j)));
%         Y(i+1)=(Y(i)*(2*m-k*dt(j)^2)+Y(i-1)*(alp*dt(j)/2-m)-m*g*dt(j)^2)/(m+alp*dt(j)/2);
    end
    d=abs(abs(Y)-abs(Y(1)))*100;
    d_max(1,j)=max(d);
%     Err_abs=abs(d_max(1,j)-d_true);
%     end
%     d_max(2,j)=M;
%     d_max(3,j)=Err_abs;
end

d_rel(:)=abs(d_max(1,1:end-1)-d_max(1,2:end));

for k=1:length(dt)-1
    for i=1:9
        if d_rel(k)<=0.5*10^-i
            d_max(2,k)=i;
        end
    end
end

for k=1:length(dt)-1
    val=num2str(d_max(1,k),d_max(2,k)+1);
    fprintf('Pour le pas de temps: %e :\n \t La distance max est: %s , avec une erreur de %e \n',dt(k),val,d_rel(k))
end

figure
loglog(dt(2:end),d_rel)
xlim([dt(1) dt(end)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation vecteur temps et espace
dt=(10e-3).*T./pi;
y_eq=-m*g/k;

t=0:dt:10*T;
Y=zeros(1,length(t));
Y(1)=y_0;
Y(2)=Y(1)*(1-(k*dt^2)/(2*m))+v_0*dt*(1-(alp*dt)/(2*m))-(g*dt^2)/2;
for i=2:length(t)-1
    Y(i+1)=(Y(i)*(2*m/dt^2-k)+Y(i-1)*(alp/(2*dt)-m/dt^2)-m*g)/(m/dt^2+alp/(2*dt));
end
D=Y-y_eq;
D_tresh=max(D)/10;

sta=zeros(1,length(t))+D_tresh;
figure
plot(t,D)
hold on
plot(t,sta)

