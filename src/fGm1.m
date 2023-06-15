function [ M1 ] = fGm1( argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
M1 = zeros(6*(N+1));
%Matriz de inércia:
T_M1j_elem = repositorio.T_M1j_elem;
T_M1j_c = repositorio.T_M1j_c;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    M1_elem = T_M1j_elem( le );
    M1_c = T_M1j_c( j, le );
    M1(n+1:n+12,n+1:n+12) = M1(n+1:n+12,n+1:n+12) + M1_elem + M1_c;
end
end