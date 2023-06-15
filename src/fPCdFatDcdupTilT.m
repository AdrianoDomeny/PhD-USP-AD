function [ dFatdupT ] = fPCdFatDcdupTilT( uTil, upTil, argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Id8 = argumentos.Id8;
R3 = argumentos.R3;
fNu3 = repositorio.fNu3; %fNu3 = @(u3)
fdmiDcAtrdtp = repositorio.dmidtXpDc; %dmiAtrdtp = @(tetap3)
fteta3 = repositorio.fteta3;
fu3 = repositorio.fu3;

fN = @(uTil)(fNu3(fu3(uTil)));
dtetap3dupTilT = Id8(8,:);
fdTatdupTilT = @(uTil,upTil)(fdmiDcAtrdtp(fteta3(upTil))*dtetap3dupTilT*fN(uTil)*R3); 

dFatdupT = [zeros(1,8);
            zeros(1,8);
            zeros(1,8);
            zeros(1,8);
            zeros(1,8);
            zeros(1,8);
            zeros(1,8);
            fdTatdupTilT(uTil,upTil)];
end