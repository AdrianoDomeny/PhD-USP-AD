function [ LmbTil ] = fLmbTil( alfaTil, argumentos )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
I1 = argumentos.Ip1_hat;
I2b = argumentos.Ip2b_hat;
l2 = argumentos.l2_hat;
dm = argumentos.dm_hat;
K1 = argumentos.K1_hat;
K2 = argumentos.K2_hat;
alfa11 = alfaTil(1,1);
alfa12 = alfaTil(2,1);
beta11 = alfaTil(3,1);
beta12 = alfaTil(4,1);
b = alfa11*K1/(beta11*(I2b + dm*l2^2)) +...
    alfa12*K2/(beta11*(I2b + dm*l2^2)) +...
    alfa12*K2/(beta12*I1);
c = (alfa11*K1*alfa12*K2)/(beta11*(I2b + dm*l2^2)*beta12*I1);
LmbTil = (1/2)*[b - sqrt(b^2 - 4*c); 
    b + sqrt(b^2 - 4*c)];
end

