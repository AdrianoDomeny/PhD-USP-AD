function [ Z ] = fiQsiDifU( X, Z0, repositorio )
%UNTITLED7 Summary of this function goes here
%   Valor da inversa de QsiDifTx na proximidade de Z0
QsiDifU = repositorio.QsiDifU;
dQsiDifUdZT = repositorio.dQsiDifUdZT;
Eta0 = Z0;
Psi = @(eta)(QsiDifU(eta) - X);
dPsidEtaT = @(eta)(dQsiDifUdZT(eta));
[ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
Z = EtaConv;
end