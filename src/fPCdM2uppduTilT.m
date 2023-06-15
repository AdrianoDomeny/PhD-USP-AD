function [ dM2uppduTilT ] = fPCdM2uppduTilT( uTil, uppTil, argumentos, repositorio )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here
Id8 = argumentos.Id8;
dm = argumentos.dm;
l2 = argumentos.l2;
alfaAc = argumentos.alfa_m2;
fvL = repositorio.fvL;
fwL = repositorio.fwL;
fteta2 = repositorio.fteta2;
dtetapp2duppTilT = Id8(6,:);
dM2uppduTilT = [zeros(1,8);
                zeros(1,8);
                zeros(1,8);
                fteta2(uppTil)*(-dm*l2*cos(alfaAc + fteta2(uTil))*dtetapp2duppTilT);
                fteta2(uppTil)*(-dm*l2*sin(alfaAc + fteta2(uTil))*dtetapp2duppTilT);
                fvL(uppTil)*(-dm*l2*cos(alfaAc + fteta2(uTil))*dtetapp2duppTilT) + fwL(uppTil)*(-dm*l2*sin(alfaAc + fteta2(uTil))*dtetapp2duppTilT);
                zeros(1,8);
                zeros(1,8)];
end