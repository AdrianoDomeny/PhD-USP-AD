function [ Z ] = fiQsiDifTx( X, y2, yp2, y2p2, y3p2, Z0, argumentos, repositorio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
QsiDifTx = repositorio.QsiDifTx;
dQsiDifTxdZT = repositorio.dQsiDifTxdZT;
Eta0 = Z0;
Psi = @(eta)(QsiDifTx(eta,y2,yp2,y2p2,y3p2) - X);
dPsidEtaT = @(eta)(dQsiDifTxdZT(eta,y2,yp2,y2p2,y3p2));
[ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
Z = EtaConv;
end

