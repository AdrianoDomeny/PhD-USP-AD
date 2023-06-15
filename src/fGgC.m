function [ G ] = fGgC( u, up,...
    argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
G = zeros(6*(N+1),1);
T_Gj_elem = repositorio.T_Gj_elem;
T_Gj_c = repositorio.T_Gj_c;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    G_elem = T_Gj_elem( le, uj, upj );
    G_c = T_Gj_c( j, le, uj, upj );
    G(n+1:n+12,1) = G(n+1:n+12,1) + G_elem + G_c;
end
end