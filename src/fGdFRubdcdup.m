function [ dFRubdupT ] = fGdFRubdcdup( u, up,...
    argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
dFRubdupT = zeros(6*(N+1));
Rub_dFRubjdupjT_c = @(nElem, le, uj, upj)(fEFdFRubdcjdupjTc( nElem, le, uj,...
    upj, argumentos, repositorio ));
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    dFRubjdupjT_c = Rub_dFRubjdupjT_c( j, le, uj, upj );
    dFRubdupT(n+1:n+12,n+1:n+12) = dFRubdupT(n+1:n+12,n+1:n+12) + dFRubjdupjT_c;
end
end