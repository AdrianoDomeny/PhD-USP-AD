clear
clc

syms K1 K2 I1 Ib20 sig

K = [K1+K2, -K2; -K2, K2];
M = diag([Ib20,I1]);

I12 = (I1^-1 + Ib20^-1)^-1;

b = K1/Ib20 + K2/I12;
c = K1*K2/(Ib20*I1);
%pC = sig^2 - b*sig + c;
%{
pCRef = simplify(det(M\K - sig*eye(2)));
erro = simplify(pCRef - pC);
ad = 1;
%}

syms lmb
pC = lmb^2 - b*lmb + c;
invMK = [(K1 + K2)/Ib20, -K2/Ib20;
                 -K2/I1,    K2/I1];
%{
invMKref = M\K;
erro = simplify(invMKref - invMK);
ad = 1;
%ok!
%}
iMKmLmbId = [(K1 + K2)/Ib20 - lmb,       -K2/Ib20;
                           -K2/I1,    K2/I1 - lmb];
iMKmLmbId1 = [ 1, -(K2/Ib20)/((K1 + K2)/Ib20 - lmb);
              -1,                   1 - lmb*(I1/K2)];
%{
iMKmLmbId1Ref = diag([((K1 + K2)/Ib20 - lmb)^-1, (K2/I1)^-1])*iMKmLmbId;
erro = simplify(iMKmLmbId1Ref - iMKmLmbId1);
ad = 1;
%ok!
%}
iMKmLmbId2 = [1,                                             -(K2/Ib20)/((K1 + K2)/Ib20 - lmb);
              0, ((1 - lmb*(I1/K2))*((K1 + K2)/Ib20 - lmb) - (K2/Ib20))/((K1 + K2)/Ib20 - lmb)];
%{
iMKmLmbId2ref = [1, 0; 1, 1]*iMKmLmbId1;
erro = simplify(iMKmLmbId2ref - iMKmLmbId2);
ad = 1;
%ok!
%}
dC = (K1 + K2)/Ib20 - lmb;
iMKmLmbId3 = [1, -(K2/Ib20)/dC;
              0,         pC/dC];
%{
iMKmLmbId3ref = [1, 0; 0, K2/I1]*iMKmLmbId2;
erro = simplify(iMKmLmbId3ref - iMKmLmbId3);
ad = 1;
%ok!
%}
%{
erro = simplify(pC - ((K2/I1 - lmb)*((K1 + K2)/Ib20 - lmb) - (K2^2/(I1*Ib20))));
ad = 1;
%ok!
%}

%Raizes do polinômio característico:
alfaSig = Ib20/I1;
alfa = 1 + 2*alfaSig;
%Autovalores:
eta1 = 1 - alfa;
eta2 = 1 + alfa;
%Autovetores:
p1 = (1 + alfa^2/(1 - eta1)^2)^(-1/2);
p2 = (1 + alfa^2/(1 - eta2)^2)^(-1/2);
P1 = [alfa/(1 - eta1); 1]*p1;
P2 = [alfa/(1 - eta2); 1]*p2;
P = [P1, P2];

betaSig = 1 + 2*alfaSig - 2*sqrt(alfaSig*(1 + alfaSig));
%{
betaSigRef = (P(1,1) + P(1,2)*sqrt(alfaSig/(alfaSig + 1)))/...
    (P(2,1) + P(2,2)*sqrt(alfaSig/(alfaSig + 1)));
erro = simplify(betaSigRef - betaSigRef);
ad = 1;
%}

syms sig2
sig1 = sig2*betaSig;
deltaK2 = simplify(I12^2*(sig1 + sig2)^2 - 4*Ib20*I12*sig1*sig2);
%{
%deltaK2 == 0
ad = 1;
%ok!
%}

K2 = (1/2)*I12*(sig1 + sig2);
K1 = 2*(I1 + Ib20)*(sig1^-1 + sig2^-1)^-1;
%{
K1Ref = sig1*sig2*I1*Ib20/K2;
erro = simplify(K1Ref - K1);
ad = 1;
%ok!
%}






