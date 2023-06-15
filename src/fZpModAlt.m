function [ Zp ] = fZpModAlt( t, Z, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fTx = repositorio.fTx;
AtTr = argumentos.AtTr;
BbT = argumentos.BbT0;
At = argumentos.At;
kNM1 = argumentos.kNM1;
R = argumentos.R_NM1;
K1 = argumentos.K1mm;
micAtr = repositorio.micAtr;
w0 = argumentos.w0;

micMax = micAtr(w0);
ub30 = 10^-5;
y20 = kNM1*ub30;
TatMax = micMax*y20*R;
miCc = micAtr(10^7*w0);
Tatc = miCc*y20*R;
Fc = [0; 0; 0; 0; 0; Tatc];

tX1 = Z(1,1);
tX2 = Z(2,1);
if K1*abs(tX1-tX2)<=TatMax
    Zp = AtTr*Z + BbT*fTx(t);
else
    Zp = At*Z - Fc + BbT*fTx(t);
end
t
end

