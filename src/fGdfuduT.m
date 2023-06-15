function [ dFuduT ] = fGdfuduT( u, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
dFuduT = zeros(6*(N+1));
U_dFujdujT_elem = repositorio.U_dFuj_dujT_elem;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    dFujdujT_elem = U_dFujdujT_elem( le, uj );
    dFuduT(n+1:n+12,n+1:n+12) = dFuduT(n+1:n+12,n+1:n+12) + dFujdujT_elem;
end
end