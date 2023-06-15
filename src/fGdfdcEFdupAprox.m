function [ dfEFdup ] = fGdfdcEFdupAprox( u, up, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
indG = argumentos.indG;
[ dGTudup ] = fGdgCdup( u, up, argumentos, repositorio );
%Matriz giroscópica:
dGdup = dGTudup;
%Torque de atrito:
[ dFatdup ] = fGdFdcAtrdup( u, up, argumentos, repositorio );

dfEFdup = indG*dGdup + dFatdup;
end