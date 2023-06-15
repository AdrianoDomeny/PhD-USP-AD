function [ dFRubcduT ] = fGdFRubdcdu( u, up,...
    argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
dFRubcduT = zeros(6*(N+1));
Rub_dFRubjdujT_c = @(nElem, le, uj, upj)(fEFdFRubdcjdujTc( nElem, le, uj,...
    upj, argumentos, repositorio ));
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    dFRubjdujT_c = Rub_dFRubjdujT_c( j, le, uj, upj );
    dFRubcduT(n+1:n+12,n+1:n+12) = dFRubcduT(n+1:n+12,n+1:n+12) + dFRubjdujT_c;
end
end