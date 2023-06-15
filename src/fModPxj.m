function [ Pxj, jP, tMax ] = fModPxj( Abvt, vTilxj, AsimPd )
%UNTITLED2 Summary of this function goes here
%   Abvt: matriz de autovetores normalizados e ordenados
%   vTilxj: vetor de proje��es unit�rias
%   tTil: vetor de cos's dos �ngulos entre vTilxj e as colunas de Abvt
%   tMax: m�x cos, correspondente ao autovetor Pxj, mais alinhado a vTilxj
%   jP: posi��o do autovetor Pxj em Abvt
%   AsimPd: matriz quadrada sim�trica positiva definida, para normaliza��o
tTil = abs(vTilxj.'*AsimPd*Abvt);
[tMax,jMax] = max(tTil);
Id = eye(size(Abvt));
jP = jMax;
Pxj = Abvt*Id(:,jP);
end

