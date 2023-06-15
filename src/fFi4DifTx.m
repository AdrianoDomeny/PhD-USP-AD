function [ fi4Dif ] = fFi4DifTx( Z, y2, yp2, y2p2, y3p2, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
I1 = argumentos.I1mm;
Ib2 = argumentos.I2mm;
l2 = argumentos.l2;
dm = argumentos.dm;
K1 = argumentos.K1mm;
K2 = argumentos.K2mm;
c3 = argumentos.bTxNM1;
R = argumentos.R3;
mi = repositorio.miAtr;
dmidtp = repositorio.dmiAtrdtp;
d2midtp2 = repositorio.d2miAtrdtp2;
d3midtp3 = repositorio.d3miAtrdtp3;
t1 = Z(1,1); t2 = Z(2,1); t3 = Z(3,1);
tp2 = Z(5,1); tp3 = Z(6,1);
fi4Dif = (K2*R^2*tp3*y2^2*dmidtp(tp3)^2 + 2*K2*R*c3*tp3*y2*dmidtp(tp3) + 2*K2*R^2*tp3*y2^2*mi(tp3)*d2midtp2(tp3) + 2*K2*R*c3*tp3^2*y2*d2midtp2(tp3) - 2*K2^2*R*t2*tp3*y2*d2midtp2(tp3) + 2*K2^2*R*t3*tp3*y2*d2midtp2(tp3))/I1^3 - (K2*R^2*tp2*y2^2*dmidtp(tp3)^2 + 2*K2*R*c3*tp2*y2*dmidtp(tp3) + 2*K2*R^2*tp2*y2^2*mi(tp3)*d2midtp2(tp3) - 2*K2^2*R*t2*tp2*y2*d2midtp2(tp3) + 2*K2^2*R*t3*tp2*y2*d2midtp2(tp3) + 2*K2*R*c3*tp2*tp3*y2*d2midtp2(tp3))/I1^3 - ((K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3))*(((R*y2*dmidtp(tp3) + c3)*((R*y2*(R*y2*dmidtp(tp3) + c3)*dmidtp(tp3))/I1^2 + (R*y2*d2midtp2(tp3)*(K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3)))/I1^2 + (R*c3*y2*dmidtp(tp3))/I1^2))/I1 + ((K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3))*((R^2*y2^2*mi(tp3)*d3midtp3(tp3))/I1^2 + (3*R*c3*y2*d2midtp2(tp3))/I1^2 + (3*R^2*y2^2*dmidtp(tp3)*d2midtp2(tp3))/I1^2 - (K2*R*t2*y2*d3midtp3(tp3))/I1^2 + (K2*R*t3*y2*d3midtp3(tp3))/I1^2 + (R*c3*tp3*y2*d3midtp3(tp3))/I1^2))/I1 + (R*y2*(c3^2 - I1*K2)*dmidtp(tp3))/I1^3 - (K2*R*y2*dmidtp(tp3))/I1^2 + (K2*R*tp2*y2*d2midtp2(tp3))/I1^2 - (K2*R*tp3*y2*d2midtp2(tp3))/I1^2))/I1 + ((K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3))*((3*R^2*y2*yp2*dmidtp(tp3)^2)/I1^2 - (R*y2p2*dmidtp(tp3))/I1 + (3*R*c3*yp2*dmidtp(tp3))/I1^2 + (3*R^2*y2*yp2*mi(tp3)*d2midtp2(tp3))/I1^2 - (2*K2*R*t2*yp2*d2midtp2(tp3))/I1^2 + (2*K2*R*t3*yp2*d2midtp2(tp3))/I1^2 + (2*R*c3*tp3*yp2*d2midtp2(tp3))/I1^2))/I1 + (R*y3p2*mi(tp3))/I1 + (R*yp2*(K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3))*(2*c3*dmidtp(tp3) + 2*R*y2*dmidtp(tp3)^2 - K2*t2*d2midtp2(tp3) + K2*t3*d2midtp2(tp3) + c3*tp3*d2midtp2(tp3) + 2*R*y2*mi(tp3)*d2midtp2(tp3)))/I1^3 - (R*c3*y2p2*mi(tp3))/I1^2 + (R*yp2*mi(tp3)*((R*y2*(R*y2*dmidtp(tp3) + c3)*dmidtp(tp3))/I1^2 + (R*y2*d2midtp2(tp3)*(K2*t3 - K2*t2 + c3*tp3 + R*y2*mi(tp3)))/I1^2 + (R*c3*y2*dmidtp(tp3))/I1^2))/I1 - (3*R^2*yp2^2*mi(tp3)*dmidtp(tp3))/I1^2 + (R*yp2*mi(tp3)*(c3^2 - I1*K2))/I1^3 - (R*c3*y2*mi(tp3)*(c3^2 - 2*I1*K2))/I1^4 - (3*R^2*y2*y2p2*mi(tp3)*dmidtp(tp3))/I1^2 + (2*K2*R*t2*y2p2*dmidtp(tp3))/I1^2 - (2*K2*R*t3*y2p2*dmidtp(tp3))/I1^2 + (3*K2*R*tp2*yp2*dmidtp(tp3))/I1^2 - (3*K2*R*tp3*yp2*dmidtp(tp3))/I1^2 - (2*R*c3*tp3*y2p2*dmidtp(tp3))/I1^2 + (K2*R*y2*dmidtp(tp3)*(K1*t1 - K1*t2 - K2*t2 + K2*t3))/(I1^2*(dm*l2^2 + Ib2));
end