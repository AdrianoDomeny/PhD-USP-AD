function [ FRubcj_c ] = fEFFRubcjc( nElem, le, uj, upj, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Nv = repositorio.Nv;
Nw = repositorio.Nw;
NtetaX = repositorio.NtetaX;
ad_N1 = repositorio.ad_N1;
ad_N1M1 = repositorio.ad_N1M1;
fFy = repositorio.fFyRubc;
fFz = repositorio.fFzRubc;
fT2 = repositorio.fT2cRub;

vN1M1 = (Nv(1,le)*ad_N1(nElem) + Nv(0,le)*ad_N1M1(nElem))*uj;
wN1M1 = (Nw(1,le)*ad_N1(nElem) + Nw(0,le)*ad_N1M1(nElem))*uj;
vpN1M1 = (Nv(1,le)*ad_N1(nElem) + Nv(0,le)*ad_N1M1(nElem))*upj;
wpN1M1 = (Nw(1,le)*ad_N1(nElem) + Nw(0,le)*ad_N1M1(nElem))*upj;
tetaXpN1M1 = (NtetaX(1,le)*ad_N1(nElem) + NtetaX(0,le)*ad_N1M1(nElem))*upj;
Fv = fFy(vN1M1,wN1M1,tetaXpN1M1,vpN1M1,wpN1M1);
Fw = fFz(vN1M1,wN1M1,tetaXpN1M1,vpN1M1,wpN1M1);
Tx = fT2(vN1M1,wN1M1,tetaXpN1M1,vpN1M1,wpN1M1);

FRubcj_c = Fv*(ad_N1(nElem)*Nv(1,le).' + ad_N1M1(nElem)*Nv(0,le).') +...
    Fw*(ad_N1(nElem)*Nw(1,le).' + ad_N1M1(nElem)*Nw(0,le).') +...
    Tx*(ad_N1(nElem)*NtetaX(1,le).' + ad_N1M1(nElem)*NtetaX(0,le).');
end