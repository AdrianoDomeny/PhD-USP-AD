function [ Zp ] = fZpModAlt23( t, Z, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Teta1 = repositorio.Teta1;
At23 = argumentos.At23;
AtTr23 = argumentos.AtTr23;
Bbteta23 = argumentos.Bbteta23;
Id3 = argumentos.Id3;

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
Tatc = 0*miCc*y20*R;
Fc = [0; 0; 0; Tatc];

tX1 = Id3(3,:)*Teta1(t);
tX2 = Z(1,1);
if K1*abs(tX1-tX2)<=TatMax
    Zp = AtTr23*Z + Bbteta23*Teta1(t);
else
    Zp = At23*Z - Fc + Bbteta23*Teta1(t);
end
t
end