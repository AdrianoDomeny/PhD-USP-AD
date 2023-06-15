function [ dfEFdu ] = fGdfdcEFdu( u, up, upp, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ dM2uppdu, dGTudu, dFudu ] = fGdm2udgCdfudu( u, up, upp,...
    argumentos, repositorio );
%Matriz giroscópica:
dGdu = dGTudu;
%Torque de atrito:
[ dFatdu ] = fGdFdcAtrdu( u, up, argumentos, repositorio );
%Rubbing:
[ dFrubdu ] = fGdFRubdcdu( u, up, argumentos, repositorio );

dfEFdu = dGdu + dFudu + dM2uppdu + dFatdu + dFrubdu;
end