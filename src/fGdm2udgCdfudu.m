function [ dM2uppduT, dGduT, dFuduT ] =...
    fGdm2udgCdfudu( u, up, upp,...
    argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
dM2uppduT = zeros(6*(N+1));
dGduT = zeros(6*(N+1));
dFuduT = zeros(6*(N+1));
T_dM2juppjdujT_elem = repositorio.T_dM2juppj_dujT_elem;
T_dM2juppjdujT_c = repositorio.T_dM2juppj_dujT_c;
T_dGjdujT_elem = repositorio.T_dGj_dujT_elem;
T_dGjdujT_c = repositorio.T_dGj_dujT_c;
U_dFujdujT_elem = repositorio.U_dFuj_dujT_elem;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    uppj = upp(n+1:n+12,1);
    dM2juppjdujT_elem = T_dM2juppjdujT_elem( le, uj, uppj );
    dM2juppjdujT_c = T_dM2juppjdujT_c( j, le, uj, uppj );
    dM2uppduT(n+1:n+12,n+1:n+12) = dM2uppduT(n+1:n+12,n+1:n+12) + dM2juppjdujT_elem + dM2juppjdujT_c;
    dGjdujT_elem = T_dGjdujT_elem( le, upj );
    dGjdujT_c = T_dGjdujT_c( j, le, uj, upj );
    dGduT(n+1:n+12,n+1:n+12) = dGduT(n+1:n+12,n+1:n+12) + dGjdujT_elem + dGjdujT_c;
    dFujdujT_elem = U_dFujdujT_elem( le, uj );
    dFuduT(n+1:n+12,n+1:n+12) = dFuduT(n+1:n+12,n+1:n+12) + dFujdujT_elem;
end
end