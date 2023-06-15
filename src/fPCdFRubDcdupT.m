function [ dFRubDcdupT ] = fPCdFRubDcdupT( uTil, upTil, argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Id8 = argumentos.Id8;
R2 = argumentos.R2;
fg = argumentos.fg;
fvL = repositorio.fvL;
fwL = repositorio.fwL;
fteta2 = repositorio.fteta2;
fFn = repositorio.fFn;
fdvCfidteta2p = repositorio.fdvCfidteta2p;
fdvCfidyLp = repositorio.fdvCfidyLp;
fdvCfidzLp = repositorio.fdvCfidzLp;
fdmiIpdcdvCfi = repositorio.dmiIpdcdvCfi;
fvCfi = repositorio.fvCfi;

vL = fvL(uTil);
wL = fwL(uTil);
vpL = fvL(upTil);
wpL = fwL(upTil);
tetap2 = fteta2(upTil);
r = (vL^2 + wL^2)^(1/2);
epsZero = argumentos.epsZero;
if r <= epsZero
    r = epsZero;
end

dyLpdupT = Id8(4,:);
dzLpdupT = Id8(5,:);
dteta2pdupT = Id8(6,:);
dvCdteta2p = fdvCfidteta2p(vL,wL);
dvCdyLp = fdvCfidyLp(vL,wL);
dvCdzLp = fdvCfidzLp(vL,wL);
vCfi = fvCfi(vL,wL,tetap2,vpL,wpL);
dmudvC = fdmiIpdcdvCfi(vCfi);

dmudupT = dmudvC*(dvCdyLp*dyLpdupT + dvCdzLp*dzLpdupT + dvCdteta2p*dteta2pdupT);
Fn = fFn(r);

dFtdupT = dmudupT*Fn;

dFvdupT = dFtdupT*(wL/r);
dFwdupT = -dFtdupT*(vL/r);
dTxdupT = -dFtdupT*(R2+fg);

dFRubDcdupT = [zeros(1,8);
               zeros(1,8); 
               zeros(1,8);
               dFvdupT;
               dFwdupT;
               dTxdupT; 
               zeros(1,8); 
               zeros(1,8)];
end