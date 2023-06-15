function [ dfi2DifdZT ] = fdFi2DifTxdZT( Z, y2, yp2, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
I1 = argumentos.I1mm;
R = argumentos.R3;
K2 = argumentos.K2mm;
c3 = argumentos.bTxNM1;
mi = repositorio.miAtr;
dmidtp = repositorio.dmiAtrdtp;
d2midtp2 = repositorio.d2miAtrdtp2;
t2 = Z(2,1); t3 = Z(3,1);
tp3 = Z(6,1);
dfi2DifdZT = [ 0, (K2*R*y2*dmidtp(tp3))/I1^2, -(K2*R*y2*dmidtp(tp3))/I1^2, 0, 0, (R*yp2*dmidtp(tp3))/I1 - (R*y2*(R*y2*dmidtp(tp3) + c3)*dmidtp(tp3))/I1^2 - (R*y2*d2midtp2(tp3)*(K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3)))/I1^2 - (R*c3*y2*dmidtp(tp3))/I1^2];
end