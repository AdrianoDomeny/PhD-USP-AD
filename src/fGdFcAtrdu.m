function [ dFaduT ] = fGdFcAtrdu( u, up, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
dFaduT = zeros(6*(N+1));
Atr_dFajdujT_c = @(nElem, le, uj, upj)(fEFdFacjdujTc( nElem, le, uj, upj,...
    argumentos, repositorio ));
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    dFajdujT_c = Atr_dFajdujT_c( j, le, uj, upj );
    dFaduT(n+1:n+12,n+1:n+12) = dFaduT(n+1:n+12,n+1:n+12) + dFajdujT_c;
end
end