function [ fEF ] = fGfdcEFaprox( u, up, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
indG = argumentos.indG;
[ GTu, Fu ] = fGgCfu( u, up,...
    argumentos, repositorio );
%Matriz giroscópica:
G = GTu;
%Torque de atrito:
[ Fat ] = fGFdcAtr( u, up, argumentos, repositorio );
%Força peso:
[ W ] = fGw( argumentos, repositorio );

fEF = indG*G + Fu - W + Fat;
end