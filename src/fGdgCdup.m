function [ dGdupT ] = fGdgCdup( u, up, argumentos, repositorio )
N = argumentos.N;
X = argumentos.X;
dGdupT = zeros(6*(N+1));
T_dGjdupjT_elem = repositorio.T_dGj_dupjT_elem;
T_dGjdupjT_c = repositorio.T_dGj_dupjT_c;
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    dGjdupjT_elem = T_dGjdupjT_elem( le, uj, upj );
    dGjdupjT_c = T_dGjdupjT_c( j, le, uj, upj );
    dGdupT(n+1:n+12,n+1:n+12) = dGdupT(n+1:n+12,n+1:n+12) + dGjdupjT_elem + dGjdupjT_c;
end
end