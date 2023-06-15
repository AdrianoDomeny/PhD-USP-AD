function [ dFajdupjT_c ] = fEFdFadcjdupjTc( nElem, le, uj, upj,...
    argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
indTatr = argumentos.indTatr;
uC = argumentos.uCef;
c_NM1 = argumentos.c_NM1;
R_NM1 = argumentos.R_NM1;
Nu = repositorio.Nu;
NtetaX = repositorio.NtetaX;
ad_N = repositorio.ad_N;
miDc = repositorio.miDc;
dmidtXpDc = repositorio.dmidtXpDc;

u_NM1 = ad_N(nElem)*Nu(1,le)*uj;
tetaXp_NM1 = ad_N(nElem)*NtetaX(1,le)*upj;
%Condição de contato:
I_ip_NM1 = double(u_NM1 >= uC);
%Força normal:
FnNM1 = repositorio.FnNM1;

dup_NM1_dupjT = ad_N(nElem)*Nu(1,le);
dF_N_NM1_dupjT = I_ip_NM1*(c_NM1*dup_NM1_dupjT);
%Torque de atrito:
dtetaXpNM1dupjT = ad_N(nElem)*NtetaX(1,le);
dmi1dupjT = @(tetaXpNM1)(dmidtXpDc(tetaXpNM1)*dtetaXpNM1dupjT);
dTxNM1dupjT = indTatr*(dmi1dupjT(tetaXp_NM1)*FnNM1(u_NM1) +...
    miDc(tetaXp_NM1)*dF_N_NM1_dupjT)*R_NM1;

dFajdupjT_c = ad_N(nElem)*...
    (Nu(1,le).'*dF_N_NM1_dupjT + NtetaX(1,le).'*dTxNM1dupjT);
end