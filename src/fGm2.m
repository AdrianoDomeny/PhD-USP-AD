function [ M2 ] = fGm2( u, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
M2 = zeros(6*(N+1));
T_M2j_elem = repositorio.T_M2j_elem;
T_M2j_c = repositorio.T_M2j_c;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    M2_elem = T_M2j_elem( le, uj );
    M2_c = T_M2j_c( j, le, uj );
    M2(n+1:n+12,n+1:n+12) = M2(n+1:n+12,n+1:n+12) + M2_elem + M2_c;
end
end