function [ dGduT ] = fGdgCdu( u, up, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
dGduT = zeros(6*(N+1));
T_dGjdujT_elem = repositorio.T_dGj_dujT_elem;
T_dGjdujT_c = repositorio.T_dGj_dujT_c;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    dGjdujT_elem = T_dGjdujT_elem( le, upj );
    dGjdujT_c = T_dGjdujT_c( j, le, uj, upj );
    dGduT(n+1:n+12,n+1:n+12) = dGduT(n+1:n+12,n+1:n+12) + dGjdujT_elem + dGjdujT_c;
end
end