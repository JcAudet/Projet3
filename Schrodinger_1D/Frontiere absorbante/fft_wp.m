function [f,Phi,STD]=fft_wp(x,psy)

dx=x(2)-x(1);
nfft=100000;
f=(0:nfft/2-1)/(nfft*dx);

phi=fft(psy,nfft);
Phi=abs(phi(1:nfft/2));
Phi_norm=trapeze(Phi,f(1),f(end),length(Phi)-1);
Phi=Phi/Phi_norm;
STD=var(f,Phi);

% %% Verif
% [PKS,LOC]=findpeaks(Phi,'minpeakheight',0.5);
% inte_66=trapeze(Phi(2903:3465),f(LOC)-STD,f(LOC)+STD,length(Phi(2903:3465))-1)
% inte_99=trapeze(Phi(2340:4028),f(LOC)-3*STD,f(LOC)+3*STD,length(Phi(2340:4028))-1)

% psy_ver=zeros(1,length(x));
% for i=1:nfft/2
%     psy_ver=psy_ver+phi(i)*sin(2*pi*f(i)*x);
% end

% figure()
% hold on
% plot(f,Phi)
% plot(f,fit,'r');
% h=line([f(LOC)+STD f(LOC)+STD],ylim);
% g=line([f(LOC)+2*STD f(LOC)+2*STD],ylim);
% q=line([f(LOC)+3*STD f(LOC)+3*STD],ylim);
% r=line([f(LOC)-3*STD f(LOC)-3*STD],ylim);
% set(h,'color','r');
% set(g,'color','r');
% set(q,'color','r');
% 
% figure()
% plot(x,psy_ver)
