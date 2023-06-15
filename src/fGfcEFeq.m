function [ fEF ] = fGfcEFeq( u, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ Fu ] = fGfu( u, argumentos, repositorio );
[ W ] = fGw( argumentos, repositorio );

fEF = Fu - W;
end