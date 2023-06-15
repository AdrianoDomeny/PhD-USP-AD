function [ fEF ] = fGfdcEFbAtr( ub, ubp, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
indG = argumentos.indG;
uTilEq = argumentos.uTilEq;
[ GTu ] = fGgC( ub + uTilEq, ubp, argumentos, repositorio );
%Matriz giroscópica:
G = GTu;
%Torque de atrito:
[ Fat ] = fGFdcAtr( ub + uTilEq, ubp, argumentos, repositorio );

fEF = indG*G + Fat;