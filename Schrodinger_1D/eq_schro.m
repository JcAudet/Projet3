function diff_f = eq_schro (psy , pot)

global dx

diff_f=zeros(1,length(psy));

%% 2nd order

% diff_f(1)=0;
% % %diff_f(1) = 1i * ( ( psy(1) - 2 * psy(2) + psy(3) ) ./ ( 2*dx^2 ) - pot(1)*psy(1) );
% % 
% diff_f(2:end-1) = ( 1i / (2*dx^2) ) * ( psy(1:end-2) - 2 .* psy(2:end-1) + psy(3:end) ) - 1i .* pot(2:end-1).*psy(2:end-1);
% % 
% diff_f(end)=0;
%diff_f(end) = 1i * ( ( psy(end-2) - 2 * psy(end-1) + psy(end) ) ./ ( 2*dx^2 ) - pot(end)*psy(end));

for i=2:length(psy)-1
    diff_f(i) = 1i * ( ( psy(i+1) + psy(i-1) - 2 * psy(i) ) ./ ( 2*dx^2 ) - pot(i)*psy(i) );
end

%% 4th order

% diff_f(3:end-2) = 1i * ( ( (-1/12)*psy(1:end-4) + (4/3)*psy(2:end-3) - (5/2)*psy(3:end-2) + (4/3)*psy(4:end-1) + (-1/12)*psy(5:end) ) ./ ( 2*dx^3 ) - pot(2:end-1).*psy(2:end-1) );