clear all;
clc;

Mant=0.1;
Bin=[];
i=1;

while Mant ~= 1
    Mant=Mant*2;
    if Mant>1
        Mant=Mant-1;
        Bin=[1 Bin];
    elseif Mant<1
        Bin=[0 Bin];
    elseif Mant==1
        Bin=[0 Bin];
    end
    i=i+1;
end

Bin