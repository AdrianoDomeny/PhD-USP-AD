function [ dfEFdup ] = fGdfcEFdup( u, up, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ dGTudup ] = fGdgCdup( u, up, argumentos, repositorio );
%Matriz giroscópica:
dGdup = dGTudup;
%Torque de atrito:
[ dFatdup ] = fGdFcAtrdup( u, up, argumentos, repositorio );
%Rubbing:
[ dFrubdup ] = fGdFRubcdup( u, up, argumentos, repositorio );

dfEFdup = dGdup + dFatdup + dFrubdup;
end