function [ dfEFdup ] = fGdfcEFdubpAtr( ub, ubp, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
indG = argumentos.indG;
uTilEq = argumentos.uTilEqEF;
[ dGTudup ] = fGdgCdup( ub + uTilEq, ubp, argumentos, repositorio );
%Matriz giroscópica:
dGdup = dGTudup;
%Torque de atrito:
[ dFatdup ] = fGdFcAtrdup( ub + uTilEq, ubp, argumentos, repositorio );

dfEFdup = indG*dGdup + dFatdup;
end