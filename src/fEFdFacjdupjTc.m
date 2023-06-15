function [ dFajdupjT_c ] = fEFdFacjdupjTc( nElem, le, uj, upj,...
    argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
indTatr = argumentos.indTatr;
uC = argumentos.uCef;
kNM1 = argumentos.kNM1;
c_NM1 = argumentos.c_NM1;
R_NM1 = argumentos.R_NM1;
Nu = repositorio.Nu;
NtetaX = repositorio.NtetaX;
ad_N = repositorio.ad_N;
miNM1nElemC = repositorio.miNM1nElemC;
dmidtXpNM1nElemC = repositorio.dmidtXpNM1nElemC;

u_NM1 = ad_N(nElem)*Nu(1,le)*uj;
up_NM1 = ad_N(nElem)*Nu(1,le)*upj;
tetaXp_NM1 = ad_N(nElem)*NtetaX(1,le)*upj;
%Condição de contato:
I_ip_NM1 = double(u_NM1 >= uC);
%Força normal:
F_N_NM1 = I_ip_NM1*(kNM1*(u_NM1 - uC) + c_NM1*up_NM1);
dup_NM1_dupjT = ad_N(nElem)*Nu(1,le);
dF_N_NM1_dupjT = I_ip_NM1*(c_NM1*dup_NM1_dupjT);
%Torque de atrito:
dtetaXp_NM1_dupjT = ad_N(nElem)*NtetaX(1,le);
dmi1_dupjT = @(nElem, tetaXp_NM1)(dmidtXpNM1nElemC(nElem, tetaXp_NM1)*dtetaXp_NM1_dupjT);
dTx_NM1_dupjT = indTatr*(dmi1_dupjT(nElem, tetaXp_NM1)*F_N_NM1 +...
    miNM1nElemC(nElem, tetaXp_NM1)*dF_N_NM1_dupjT)*R_NM1;

dFajdupjT_c = ad_N(nElem)*...
    (Nu(1,le).'*dF_N_NM1_dupjT + NtetaX(1,le).'*dTx_NM1_dupjT);
end