function [ dFRubcjdupjTc ] = fEFdFRubcjdupjTc( nElem, le, uj,...
    upj, argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
R2 = argumentos.R2;
fg = argumentos.fg;
Nv = repositorio.Nv;
Nw = repositorio.Nw;
NtetaX = repositorio.NtetaX;
ad_N1 = repositorio.ad_N1;
ad_N1M1 = repositorio.ad_N1M1;
fFn = repositorio.fFn;
fdvCfidteta2p = repositorio.fdvCfidteta2p;
fdvCfidyLp = repositorio.fdvCfidyLp;
fdvCfidzLp = repositorio.fdvCfidzLp;
fdmiIpcdvCfi = repositorio.dmiIpcdvCfi;
fvCfi = repositorio.fvCfi;

vN1M1 = (Nv(1,le)*ad_N1(nElem) + Nv(0,le)*ad_N1M1(nElem))*uj;
wN1M1 = (Nw(1,le)*ad_N1(nElem) + Nw(0,le)*ad_N1M1(nElem))*uj;
vpN1M1 = (Nv(1,le)*ad_N1(nElem) + Nv(0,le)*ad_N1M1(nElem))*upj;
wpN1M1 = (Nw(1,le)*ad_N1(nElem) + Nw(0,le)*ad_N1M1(nElem))*upj;
tetaXpN1M1 = (NtetaX(1,le)*ad_N1(nElem) + NtetaX(0,le)*ad_N1M1(nElem))*upj;
r = (vN1M1^2 + wN1M1^2)^(1/2);
epsZero = argumentos.epsZero;
if r <= epsZero
    r = epsZero;
end

dyLpdupjT = ad_N1(nElem)*Nv(1,le) + ad_N1M1(nElem)*Nv(0,le);
dzLpdupjT = ad_N1(nElem)*Nw(1,le) + ad_N1M1(nElem)*Nw(0,le);
dteta2pdupjT = ad_N1(nElem)*NtetaX(1,le) + ad_N1M1(nElem)*NtetaX(0,le);
dvCdteta2p = fdvCfidteta2p(vN1M1,wN1M1);
dvCdyLp = fdvCfidyLp(vN1M1,wN1M1);
dvCdzLp = fdvCfidzLp(vN1M1,wN1M1);
vCfi = fvCfi(vN1M1,wN1M1,tetaXpN1M1,vpN1M1,wpN1M1);
dmudvC = fdmiIpcdvCfi(vCfi);

dmudupjT = dmudvC*(dvCdyLp*dyLpdupjT + dvCdzLp*dzLpdupjT + dvCdteta2p*dteta2pdupjT);
Fn = fFn(r);

dFtdupjT = dmudupjT*Fn;

dFvdupjT = dFtdupjT*(wN1M1/r);
dFwdupjT = -dFtdupjT*(vN1M1/r);
dTxdupjT = -dFtdupjT*(R2+fg);

dFRubcjdupjTc = (ad_N1(nElem)*Nv(1,le).' + ad_N1M1(nElem)*Nv(0,le).')*dFvdupjT +...
    (ad_N1(nElem)*Nw(1,le).' + ad_N1M1(nElem)*Nw(0,le).')*dFwdupjT +...
    (ad_N1(nElem)*NtetaX(1,le).' + ad_N1M1(nElem)*NtetaX(0,le).')*dTxdupjT;
end