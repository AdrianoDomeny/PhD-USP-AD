function [ fi1Dif ] = fFi1DifTx( Z, y2, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
I1 = argumentos.I1mm;
R = argumentos.R3;
mi = repositorio.miAtr;
tp3 = Z(6,1);
fi1Dif = (R*y2*mi(tp3))/I1;
end