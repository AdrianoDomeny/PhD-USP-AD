function [ W ] = fGw( argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
W = zeros(6*(N+1),1);
%Matriz força peso:
Wj_elem = repositorio.Wj_elem;
Wj_c = repositorio.Wj_c;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    W_elem = Wj_elem( le );
    W_c = Wj_c( j, le );
    W(n+1:n+12,1) = W(n+1:n+12,1) + W_elem + W_c;
end
end