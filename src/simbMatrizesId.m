clear
clc

syms x y z w
syms f1(x,y,z,w) f2(x,y,z,w) f3(x,y,z,w) f4(x,y,z,w)
syms m1 m2 m3 k1 k2 k3 c1 c2 c3

M = [m1, m2; m2, m3];
K = [k1, k2; k2, k3];
C = [c1, c2; c2, c3];

A = [zeros(2), eye(2); -M\K, -M\C];
Z = [x; y; z; w];
Qsi = [f1(x,y,z,w), f2(x,y,z,w), f3(x,y,z,w), f4(x,y,z,w)];

dZTQsiAZref = gradient(Qsi*(A*Z),Z).';
dQsiTdZT = jacobian(Qsi,Z);
dZTQsiAZ = (A*Z).'*dQsiTdZT + Qsi*A; %ok!
erro = simplify(dZTQsiAZref - dZTQsiAZ);

ad = 1;



