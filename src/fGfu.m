function [ Fu ] = fGfu( u, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
Fu = zeros(6*(N+1),1);
U_Fuj_elem = repositorio.U_Fuj_elem;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    Fu_elem = U_Fuj_elem( le, uj );
    Fu(n+1:n+12,1) = Fu(n+1:n+12,1) + Fu_elem;
end
end