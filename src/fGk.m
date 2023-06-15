function [ K ] = fGk( argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
K = zeros(6*(N+1));
%Matriz de rigidez:
U_Kj_elem = repositorio.U_Kj_elem;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    K_elem = U_Kj_elem( le );
    K(n+1:n+12,n+1:n+12) = K(n+1:n+12,n+1:n+12) + K_elem;
end
end