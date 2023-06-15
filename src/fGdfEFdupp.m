function [ dfEFdupp ] = fGdfEFdupp( u, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ M2 ] = fGm2( u, argumentos, repositorio );
dfEFdupp = M2;
end