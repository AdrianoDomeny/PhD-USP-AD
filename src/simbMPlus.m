clear
clc

syms Ip1 Ip2 Ip3 c1 c2 c3 K1 K2

M = diag([Ip3, Ip2, Ip1]);
%C = diag([c1, c2, c3]);
C = [c1, -c1, 0; -c1, c1+c2, -c2; 0, -c2, c2];
K = [K1, -K1, 0; -K1, K1+K2, -K2; 0, -K2, K2];

beta1 = Ip3/(Ip1+Ip2+Ip3);
beta2 = Ip2/(Ip1+Ip2+Ip3);
beta3 = Ip1/(Ip1+Ip2+Ip3);

Bteta = [1, -1, 0; 
    0, 1, -1; 
    beta1, beta2, beta3];
iBteta = [1, 1, 1; 
    0, 1, 1; 
    0, 0, 1] -...
    [beta1, beta1+beta2, 0; 
    beta1, beta1+beta2, 0; 
    beta1, beta1+beta2, 0];

Mx = simplify(iBteta.'*M*iBteta);
Cx = simplify(iBteta.'*C*iBteta);
Kx = simplify(iBteta.'*K*iBteta);

syms teta1 teta2 teta3 teta12 teta23 tetab

uTeta = [teta1; teta2; teta3];
uTetab = [teta12; teta23; tetab];
ad1 = Bteta*uTeta;
ad2 = simplify(iBteta*uTetab);

ad = 1;






