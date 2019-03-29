function  [V_b,V_f] = V_potentiel(x,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zeros!!

V_b=[];
V_f=[];
for j = 1 : length(y)
    for i = 1 : length(x) 
       var_b=0;
       V_b=[V_b var_b];
       var_f=0;
       V_f=[V_f var_f];
    end
end

end