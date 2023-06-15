function [ M2upp, G, Fu ] = fGm2uppgCfu( u, up, upp,...
    argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
M2upp = zeros(6*(N+1),1);
G = zeros(6*(N+1),1);
Fu = zeros(6*(N+1),1);
T_M2juppj_elem = repositorio.T_M2juppj_elem;
T_M2juppj_c = repositorio.T_M2juppj_c;
T_Gj_elem = repositorio.T_Gj_elem;
T_Gj_c = repositorio.T_Gj_c;
U_Fuj_elem = repositorio.U_Fuj_elem;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    uppj = upp(n+1:n+12,1);
    M2upp_elem = T_M2juppj_elem( le, uj, uppj );
    M2upp_c = T_M2juppj_c( j, le, uj, uppj );
    M2upp(n+1:n+12,1) = M2upp(n+1:n+12,1) + M2upp_elem + M2upp_c;
    G_elem = T_Gj_elem( le, uj, upj );
    G_c = T_Gj_c( j, le, uj, upj );
    G(n+1:n+12,1) = G(n+1:n+12,1) + G_elem + G_c;
    Fu_elem = U_Fuj_elem( le, uj );
    Fu(n+1:n+12,1) = Fu(n+1:n+12,1) + Fu_elem;
end
end