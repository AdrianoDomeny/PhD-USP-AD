function [ fdGacdupTilT ] = fPCdGacdupTilT( uTil, upTil, argumentos, repositorio )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here
Id8 = argumentos.Id8;
dm = argumentos.dm;
l2 = argumentos.l2;
alfaAc = argumentos.alfa_m2;
fteta2 = repositorio.fteta2;
dtetap2dupTilT = Id8(6,:);
fdGacdupTilT = [zeros(1,8);
                zeros(1,8); 
                zeros(1,8); 
                -dm*l2*2*fteta2(upTil)*dtetap2dupTilT*cos(fteta2(uTil) + alfaAc); 
                -dm*l2*2*fteta2(upTil)*dtetap2dupTilT*sin(fteta2(uTil) + alfaAc); 
                zeros(1,8); 
                zeros(1,8); 
                zeros(1,8)];
end