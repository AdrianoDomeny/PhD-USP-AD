function [ vTilEsp ] = fModEsp( wP, fi0, argumentos )
%UNTITLED Summary of this function goes here
%   vTilEsp: Vetor de teste com dimensão Ng X 1
%   Detecção do modo espiral
N = argumentos.N;
X = argumentos.X;
vTilTEj = @(Xj)([0; cos(wP*Xj + fi0); 0; sin(wP*Xj + fi0); 0; 0]);
vTilEsp = zeros(6*(N+1),1);
for i = 1:N+1
    n = 6*(i-1);
    vTilEsp(n+1:n+6,1) = vTilTEj(X(i));
end
end