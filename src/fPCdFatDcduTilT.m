function [ dFatduT ] = fPCdFatDcduTilT( uTil, upTil, argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Id8 = argumentos.Id8;
R3 = argumentos.R3;
fmiDcAtr = repositorio.miDc; %miAtr = @(tetap3)
fdNu3du3 = repositorio.fdNu3du3; %fdNu3du3 = @(u3)
fteta3 = repositorio.fteta3;
fu3 = repositorio.fu3;

du3duTilT = Id8(7,:);
fdNduTilT = @(uTil)(fdNu3du3(fu3(uTil))*du3duTilT);
fdTatduTilT = @(uTil,upTil)(fmiDcAtr(fteta3(upTil))*fdNduTilT(uTil)*R3); 

dFatduT = [zeros(1,8);
           zeros(1,8);
           zeros(1,8);
           zeros(1,8);
           zeros(1,8);
           zeros(1,8);
           fdNduTilT(uTil);
           fdTatduTilT(uTil,upTil)];
end