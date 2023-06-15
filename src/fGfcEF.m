function [ fEF ] = fGfcEF( u, up, upp, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ M2upp, GTu, Fu ] = fGm2uppgCfu( u, up, upp,...
    argumentos, repositorio );
%Matriz girosc�pica:
G = GTu;
%Torque de atrito:
[ Fat ] = fGFcAtr( u, up, argumentos, repositorio );
%Rubbing:
[ Frub ] = fGFRubc( u, up, argumentos, repositorio );
%For�a peso:
[ W ] = fGw( argumentos, repositorio );

fEF = G + Fu + M2upp - W + Fat + Frub;
end

