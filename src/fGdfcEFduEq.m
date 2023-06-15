function [ dfEFdu ] = fGdfcEFduEq( u, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ dFudu ] = fGdfuduT( u, argumentos, repositorio );

dfEFdu = dFudu;
end