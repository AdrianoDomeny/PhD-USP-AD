function [ dfEFdu ] = fGdfdcEFduAprox( u, up, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
indG = argumentos.indG;
[ dGTudu, dFudu ] = fGdgCdfudu( u, up,...
    argumentos, repositorio );
%Matriz giroscópica:
dGdu = dGTudu;
%Torque de atrito:
[ dFatdu ] = fGdFdcAtrdu( u, up, argumentos, repositorio );

dfEFdu = indG*dGdu + dFudu + dFatdu;
end