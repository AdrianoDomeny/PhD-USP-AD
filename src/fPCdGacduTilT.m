function [ fdGacduTilT ] = fPCdGacduTilT( uTil, upTil, argumentos, repositorio )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here
Id8 = argumentos.Id8;
dm = argumentos.dm;
l2 = argumentos.l2;
alfaAc = argumentos.alfa_m2;
fteta2 = repositorio.fteta2;
dteta2duTilT = Id8(6,:);
fdGacduTilT = [zeros(1,8); 
               zeros(1,8);
               zeros(1,8);
               -dm*l2*fteta2(upTil)^2*(-sin(fteta2(uTil) + alfaAc)*dteta2duTilT); 
               -dm*l2*fteta2(upTil)^2*cos(fteta2(uTil) + alfaAc)*dteta2duTilT; 
               zeros(1,8);
               zeros(1,8);
               zeros(1,8)];
end

