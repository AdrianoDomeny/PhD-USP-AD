function [ dfEFdu ] = fGdfcEFdubAtr( ub, ubp, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
indG = argumentos.indG;
uTilEq = argumentos.uTilEq;
[ dGTudu ] = fGdgCdu( ub + uTilEq, ubp, argumentos, repositorio );
%Matriz giroscópica:
dGdu = dGTudu;
%Torque de atrito:
[ dFatdu ] = fGdFcAtrduT( ub + uTilEq, ubp, argumentos, repositorio );

dfEFdu = indG*dGdu + dFatdu;
end