function [ dFajdujTc ] = fEFdFadcjdujTc( nElem, le, uj, upj,...
    argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
indTatr = argumentos.indTatr;
uC = argumentos.uCef;
kNM1 = argumentos.kNM1;
R_NM1 = argumentos.R_NM1;
Nu = repositorio.Nu;
NtetaX = repositorio.NtetaX;
ad_N = repositorio.ad_N;
miDc = repositorio.miDc;

u_NM1 = ad_N(nElem)*Nu(1,le)*uj;
tetaXpNM1 = ad_N(nElem)*NtetaX(1,le)*upj;
%Condição de contato:
I_ip_NM1 = double(u_NM1 >= uC);
%Força normal:
du_NM1_dujT = ad_N(nElem)*Nu(1,le);
dF_N_NM1_dujT = I_ip_NM1*(kNM1*du_NM1_dujT);
%Torque de atrito:
dTxdujT = indTatr*miDc(tetaXpNM1)*dF_N_NM1_dujT*R_NM1;

dFajdujTc = ad_N(nElem)*...
    (Nu(1,le).'*dF_N_NM1_dujT + NtetaX(1,le).'*dTxdujT);
end