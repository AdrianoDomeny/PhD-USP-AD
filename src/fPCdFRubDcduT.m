function [ dFRubcduT ] = fPCdFRubDcduT( uTil, upTil, argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Id8 = argumentos.Id8;
R2 = argumentos.R2;
fg = argumentos.fg;
fvL = repositorio.fvL;
fwL = repositorio.fwL;
fteta2 = repositorio.fteta2;
fFn = repositorio.fFn;
miIpDc = repositorio.miIpDc;
fdFnduL = repositorio.dFnduL;
fduLdyL = repositorio.fduLdyL;
fduLdzL = repositorio.fduLdzL;
fdyLuLpwm1dyL = repositorio.fdyLuLpwm1dyL; %fdyLuLpwm1dyL = @(yL,zL)
fdyLuLpwm1dzL = repositorio.fdyLuLpwm1dzL; %fdyLuLpwm1dzL = @(yL,zL)
fdzLuLpwm1dyL = repositorio.fdzLuLpwm1dyL; %fdzLuLpwm1dyL = @(yL,zL)
fdzLuLpwm1dzL = repositorio.fdzLuLpwm1dzL; %fdzLuLpwm1dzL = @(yL,zL)
fvCfi = repositorio.fvCfi;
fFtDc = repositorio.fFtDc;
fdmiIpDcdvCfi = repositorio.dmiIpDcdvCfi;
fdvCfidyL = repositorio.dvCfidyL;
fdvCfidzL = repositorio.dvCfidzL;

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

dFnduL = fdFnduL(r);
duLdyL = fduLdyL(vL,wL);
duLdzL = fduLdzL(vL,wL);
dyLuLpwm1dyL = fdyLuLpwm1dyL(vL,wL);
dyLuLpwm1dzL = fdyLuLpwm1dzL(vL,wL);
Fn = fFn(r);
vCfi = fvCfi(vL,wL,tetap2,vpL,wpL);
mu = miIpDc(vCfi);
dyLduT = Id8(4,:);
dzLduT = Id8(5,:);
dFnduT = dFnduL*(duLdyL*dyLduT + duLdzL*dzLduT);
dyLuLpwm1duT = dyLuLpwm1dyL*dyLduT + dyLuLpwm1dzL*dzLduT;
dmudvC = fdmiIpDcdvCfi(vCfi);
dvCdyL = fdvCfidyL(vL,wL,tetap2,vpL,wpL);
dvCdzL = fdvCfidzL(vL,wL,tetap2,vpL,wpL);
dmuduT = dmudvC*(dvCdyL*dyLduT + dvCdzL*dzLduT);
dFtduT = dmuduT*Fn + mu*dFnduT;
dzLuLpwm1duT = fdzLuLpwm1dyL(vL,wL)*dyLduT +...
    fdzLuLpwm1dzL(vL,wL)*dzLduT;
Ft = fFtDc(r,vCfi);

dFvduT = -dFnduT*(vL/r) - Fn*dyLuLpwm1duT +...
    dFtduT*(wL/r) + Ft*dzLuLpwm1duT;
dFwduT = -dFnduT*(wL/r) - Fn*dzLuLpwm1duT -...
    dFtduT*(vL/r) - Ft*dyLuLpwm1duT;
dTxduT = -dFtduT*(R2+fg);

dFRubcduT = [zeros(1,8);
             zeros(1,8); 
             zeros(1,8);
             dFvduT;
             dFwduT;
             dTxduT; 
             zeros(1,8); 
             zeros(1,8)];
end