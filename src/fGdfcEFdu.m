function [ dfEFdu ] = fGdfcEFdu( u, up, upp, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ dM2uppdu, dGTudu, dFudu ] = fGdm2udgCdfudu( u, up, upp,...
    argumentos, repositorio );
%Matriz giroscópica:
dGdu = dGTudu;
%Torque de atrito:
[ dFatdu ] = fGdFcAtrdu( u, up, argumentos, repositorio );
%Rubbing:
[ dFrubdu ] = fGdFRubcdu( u, up, argumentos, repositorio );

dfEFdu = dGdu + dFudu + dM2uppdu + dFatdu + dFrubdu;
end