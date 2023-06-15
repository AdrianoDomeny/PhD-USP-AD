function [ Fa ] = fGFcAtr( u, up, argumentos, repositorio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
Fa = zeros(6*(N+1),1);
Atr_Faj_c = @(nElem, le, uj, upj)(fEFFacjc( nElem, le, uj, upj,...
    argumentos, repositorio ));
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    Faj_c = Atr_Faj_c( j, le, uj, upj );
    Fa(n+1:n+12,1) = Fa(n+1:n+12,1) + Faj_c;
end
end