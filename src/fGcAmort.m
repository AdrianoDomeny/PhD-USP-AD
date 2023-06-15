function [ C ] = fGcAmort( argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
C = zeros(6*(N+1),6*(N+1));
%Matriz força peso:
Cj_c = repositorio.Cj_c;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    C_c = Cj_c( j, le );
    C(n+1:n+12,n+1:n+12) = C(n+1:n+12,n+1:n+12) + C_c;
end
end