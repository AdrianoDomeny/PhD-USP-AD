function [ dFnrduT ] = fPCdFnrduTilT( uTil, argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Id8 = argumentos.Id8;
fdNu3du3 = repositorio.fdNu3du3; %fdNu3du3 = @(u3)
fu3 = repositorio.fu3;

du3duTilT = Id8(7,:);
fdNduTilT = @(uTil)(fdNu3du3(fu3(uTil))*du3duTilT);

dFnrduT = [zeros(1,8);
           zeros(1,8);
           zeros(1,8);
           zeros(1,8);
           zeros(1,8);
           zeros(1,8);
           fdNduTilT(uTil);
           zeros(1,8)];
end