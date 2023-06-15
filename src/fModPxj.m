function [ Pxj, jP, tMax ] = fModPxj( Abvt, vTilxj, AsimPd )
%UNTITLED2 Summary of this function goes here
%   Abvt: matriz de autovetores normalizados e ordenados
%   vTilxj: vetor de projeções unitárias
%   tTil: vetor de cos's dos ângulos entre vTilxj e as colunas de Abvt
%   tMax: máx cos, correspondente ao autovetor Pxj, mais alinhado a vTilxj
%   jP: posição do autovetor Pxj em Abvt
%   AsimPd: matriz quadrada simétrica positiva definida, para normalização
tTil = abs(vTilxj.'*AsimPd*Abvt);
[tMax,jMax] = max(tTil);
Id = eye(size(Abvt));
jP = jMax;
Pxj = Abvt*Id(:,jP);
end

