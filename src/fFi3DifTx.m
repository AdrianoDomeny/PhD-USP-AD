function [ fi3Dif ] = fFi3DifTx( Z, y2, yp2, y2p2, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
I1 = argumentos.I1mm;
K2 = argumentos.K2mm;
c3 = argumentos.bTxNM1;
R = argumentos.R3;
mi = repositorio.miAtr;
dmidtp = repositorio.dmiAtrdtp;
d2midtp2 = repositorio.d2miAtrdtp2;
t2 = Z(2,1); t3 = Z(3,1);
tp2 = Z(5,1); tp3 = Z(6,1);
fi3Dif = (((R*y2*(R*y2*dmidtp(tp3) + c3)*dmidtp(tp3))/I1^2 + (R*y2*d2midtp2(tp3)*(K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3)))/I1^2 + (R*c3*y2*dmidtp(tp3))/I1^2)*(K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3)))/I1 - (R*yp2*(c3*mi(tp3) - K2*t2*dmidtp(tp3) + K2*t3*dmidtp(tp3) + c3*tp3*dmidtp(tp3) + 2*R*y2*mi(tp3)*dmidtp(tp3)))/I1^2 + (R*y2p2*mi(tp3))/I1 - (R*yp2*dmidtp(tp3)*(K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3)))/I1^2 + (R*y2*mi(tp3)*(c3^2 - I1*K2))/I1^3 + (K2*R*tp2*y2*dmidtp(tp3))/I1^2 - (K2*R*tp3*y2*dmidtp(tp3))/I1^2;
end