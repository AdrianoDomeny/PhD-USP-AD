function [ dGduT, dFuduT ] = fGdgCdfudu( u, up, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
dGduT = zeros(6*(N+1));
dFuduT = zeros(6*(N+1));
T_dGjdujT_elem = repositorio.T_dGj_dujT_elem;
T_dGjdujT_c = repositorio.T_dGj_dujT_c;
U_dFujdujT_elem = repositorio.U_dFuj_dujT_elem;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    dGjdujT_elem = T_dGjdujT_elem( le, upj );
    dGjdujT_c = T_dGjdujT_c( j, le, uj, upj );
    dGduT(n+1:n+12,n+1:n+12) = dGduT(n+1:n+12,n+1:n+12) + dGjdujT_elem + dGjdujT_c;
    dFujdujT_elem = U_dFujdujT_elem( le, uj );
    dFuduT(n+1:n+12,n+1:n+12) = dFuduT(n+1:n+12,n+1:n+12) + dFujdujT_elem;
end
end