function [ dfi1DifdZT ] = fdFi1DifTxdZT( Z, y2, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
I1 = argumentos.I1mm;
R = argumentos.R3;
dmidtp = repositorio.dmiAtrdtp;
tp3 = Z(6,1);
dfi1DifdZT = [ 0, 0, 0, 0, 0, (R*y2*dmidtp(tp3))/I1];
end