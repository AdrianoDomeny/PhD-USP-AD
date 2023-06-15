function [ F_Rub ] = fGFRubc( u, up, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = argumentos.N;
X = argumentos.X;
F_Rub = zeros(6*(N+1),1);
Rub_F_Rubj = @(nElem, le, uj, upj)(fEFFRubcjc( nElem, le, uj,...
    upj, repositorio ));
for j = 1:N
    le = X(j+1) - X(j);
    n = 6*(j-1);
    uj = u(n+1:n+12,1);
    upj = up(n+1:n+12,1);
    F_Rubj = Rub_F_Rubj( j, le, uj, upj );
    F_Rub(n+1:n+12,1) = F_Rub(n+1:n+12,1) + F_Rubj;
end
end