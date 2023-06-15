function [ dFRubdcjdujT_c ] = fEFdFRubdcjdujTc( nElem, le, uj,...
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
miIpdc = repositorio.miIpdc;
fdFnduL = repositorio.dFnduL;
fduLdyL = repositorio.fduLdyL;
fduLdzL = repositorio.fduLdzL;
fdyLuLpwm1dyL = repositorio.fdyLuLpwm1dyL;
fdyLuLpwm1dzL = repositorio.fdyLuLpwm1dzL;
dzLuLpwm1dyL = repositorio.fdzLuLpwm1dyL;
dzLuLpwm1dzL = repositorio.fdzLuLpwm1dzL;
fdmiIpdcdvCfi = repositorio.dmiIpdcdvCfi;
fvCfi = repositorio.fvCfi;
fFtdc = repositorio.fFtdc;
fdvCfidyL = repositorio.dvCfidyL;
fdvCfidzL = repositorio.dvCfidzL;

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

dFnduL = fdFnduL(r);
duLdyL = fduLdyL(vN1M1,wN1M1);
duLdzL = fduLdzL(vN1M1,wN1M1);
dyLuLpwm1dyL = fdyLuLpwm1dyL(vN1M1,wN1M1);
dyLuLpwm1dzL = fdyLuLpwm1dzL(vN1M1,wN1M1);
Fn = fFn(r);
vCfi = fvCfi(vN1M1,wN1M1,tetaXpN1M1,vpN1M1,wpN1M1);
mu = miIpdc(vCfi);
dyLdujT = (ad_N1(nElem)*Nv(1,le) + ad_N1M1(nElem)*Nv(0,le));
dzLdujT = (ad_N1(nElem)*Nw(1,le) + ad_N1M1(nElem)*Nw(0,le));
dFndujT = dFnduL*(duLdyL*dyLdujT + duLdzL*dzLdujT);
dyLuLpwm1dujT = dyLuLpwm1dyL*dyLdujT + dyLuLpwm1dzL*dzLdujT;
dmudvC = fdmiIpdcdvCfi(vCfi);
dvCdyL = fdvCfidyL(vN1M1,wN1M1,tetaXpN1M1,vpN1M1,wpN1M1);
dvCdzL = fdvCfidzL(vN1M1,wN1M1,tetaXpN1M1,vpN1M1,wpN1M1);
dmudujT = dmudvC*(dvCdyL*dyLdujT + dvCdzL*dzLdujT);
dFtdujT = dmudujT*Fn + mu*dFndujT;
dzLuLpwm1dujT = dzLuLpwm1dyL(vN1M1,wN1M1)*dyLdujT + dzLuLpwm1dzL(vN1M1,wN1M1)*dzLdujT;
Ft = fFtdc(r,vCfi);

dFvdujT = -dFndujT*(vN1M1/r) - Fn*dyLuLpwm1dujT +...
    dFtdujT*(wN1M1/r) + Ft*dzLuLpwm1dujT;
dFwdujT = -dFndujT*(wN1M1/r) - Fn*dzLuLpwm1dujT -...
    dFtdujT*(vN1M1/r) - Ft*dyLuLpwm1dujT;
dTxdujT = -dFtdujT*(R2+fg);

dFRubdcjdujT_c = (ad_N1(nElem)*Nv(1,le).' + ad_N1M1(nElem)*Nv(0,le).')*dFvdujT +...
    (ad_N1(nElem)*Nw(1,le).' + ad_N1M1(nElem)*Nw(0,le).')*dFwdujT +...
    (ad_N1(nElem)*NtetaX(1,le).' + ad_N1M1(nElem)*NtetaX(0,le).')*dTxdujT;
end