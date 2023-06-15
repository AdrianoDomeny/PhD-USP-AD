function [ dFadupT ] = fGdFcAtrdup( u, up, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
dFadupT = zeros(6*(N+1));
Atr_dFajdupjT_c = @(nElem, le, uj, upj)(fEFdFacjdupjTc( nElem, le, uj, upj,...
    argumentos, repositorio ));
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    dFajdupjT_c = Atr_dFajdupjT_c( j, le, uj, upj );
    dFadupT(n+1:n+12,n+1:n+12) = dFadupT(n+1:n+12,n+1:n+12) + dFajdupjT_c;
end
end