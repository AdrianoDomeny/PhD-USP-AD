function [ dfEFdup ] = fGdfdcEFdup( u, up, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ dGTudup ] = fGdgCdup( u, up, argumentos, repositorio );
%Matriz giroscópica:
dGdup = dGTudup;
%Torque de atrito:
[ dFatdup ] = fGdFdcAtrdup( u, up, argumentos, repositorio );
%Rubbing:
[ dFrubdup ] = fGdFRubdcdup( u, up, argumentos, repositorio );

dfEFdup = dGdup + dFatdup + dFrubdup;
end