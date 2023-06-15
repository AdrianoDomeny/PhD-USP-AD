function [ Fu ] = fPCFu( uTil, argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
L1 = argumentos.L1;
L = argumentos.L;
k1 = argumentos.k1mm;
k2 = argumentos.k2mm;
fu1 = repositorio.fu1;
fu2 = repositorio.fu2;
fvL = repositorio.fvL;
fwL = repositorio.fwL;
fu3 = repositorio.fu3;
Fu = [k1*( (fvL(uTil)^2 + fwL(uTil)^2)/(L1 - fu1(uTil) + fu2(uTil))^2 )*( fu2(uTil) - fu1(uTil) + (fvL(uTil)^2 + fwL(uTil)^2)/(L1 - fu1(uTil) + fu2(uTil)) ) - k1*(fvL(uTil)^2 + fwL(uTil)^2)/(L1 - fu1(uTil) + fu2(uTil)); 
    0; 
    k2*( (fvL(uTil)^2 + fwL(uTil)^2)/(L - L1 - fu2(uTil) + fu3(uTil))^2 )*( fu3(uTil) - fu2(uTil) ) + k2*( (fvL(uTil)^2 + fwL(uTil)^2)/(L - L1 - fu2(uTil) + fu3(uTil))^2 - 1 )*( (fvL(uTil)^2 + fwL(uTil)^2)/(L - L1 - fu2(uTil) + fu3(uTil)) ) + k1*( -(fvL(uTil)^2 + fwL(uTil)^2)/(L1 - fu1(uTil) + fu2(uTil))^2 )*( fu2(uTil) - fu1(uTil) ) + k1*( 1 - (fvL(uTil)^2 + fwL(uTil)^2)/(L1 - fu1(uTil) + fu2(uTil))^2 )*( (fvL(uTil)^2 + fwL(uTil)^2)/(L1 - fu1(uTil) + fu2(uTil)) );
    2*k2*fvL(uTil)*( fu3(uTil) - fu2(uTil) + (fvL(uTil)^2 + fwL(uTil)^2)/(L - L1 - fu2(uTil) + fu3(uTil)) )/( L - L1 - fu2(uTil) + fu3(uTil) ) + 2*k1*fvL(uTil)*( fu2(uTil) - fu1(uTil) + (fvL(uTil)^2 + fwL(uTil)^2)/(L1 - fu1(uTil) + fu2(uTil)) )/( L1 - fu1(uTil) + fu2(uTil) );
    2*k2*fwL(uTil)*( fu3(uTil) - fu2(uTil) + (fvL(uTil)^2 + fwL(uTil)^2)/(L - L1 - fu2(uTil) + fu3(uTil)) )/( L - L1 - fu2(uTil) + fu3(uTil) ) + 2*k1*fwL(uTil)*( fu2(uTil) - fu1(uTil) + (fvL(uTil)^2 + fwL(uTil)^2)/(L1 - fu1(uTil) + fu2(uTil)) )/( L1 - fu1(uTil) + fu2(uTil) ); 
    0; 
    k2*( -(fvL(uTil)^2 + fwL(uTil)^2)/(L - L1 - fu2(uTil) + fu3(uTil))^2 )*( fu3(uTil) - fu2(uTil) ) + k2*( 1 - (fvL(uTil)^2 + fwL(uTil)^2)/(L - L1 - fu2(uTil) + fu3(uTil))^2 )*( (fvL(uTil)^2 + fwL(uTil)^2)/(L - L1 - fu2(uTil) + fu3(uTil)) ); 
    0];
end