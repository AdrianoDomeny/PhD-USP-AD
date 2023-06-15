clear
clc

syms Ibp3 Ibp2 Ip1 K1 K2 lmb c1 c2 c3

Id3 = eye(3);
M = diag([Ibp3, Ibp2, Ip1]);
K = [K1, -K1, 0; -K1, K1+K2, -K2; 0, -K2, K2];
C = diag([c1, c2, c3]);
B = Id3(:,1);

%%
%Equações de estado: Zp = A*Z + Bb*U
A = [zeros(3), Id3; -M\K, -M\C];
Bb = [zeros(3,1); M\B];
C1 = [zeros(1,3), Id3(3,:)];

beta0 = C1;
beta1 = beta0*A;
beta0Bb = beta0*Bb;
beta2 = beta1*A;
beta1Bb = beta1*Bb;
beta3 = beta2*A;
beta2Bb = beta2*Bb;
beta4 = beta3*A;
beta3Bb = beta3*Bb;
beta5 = beta4*A;
beta4Bb = beta4*Bb;

%%
%Modos de vibração:
I12 = (Ip1^-1 + Ibp2^-1)^-1;
I23 = (Ibp2^-1 + Ibp3^-1)^-1;
I123 = (1/(Ibp3*Ip1) + 1/(Ibp2*Ip1) + 1/(Ibp2*Ibp3))^-1;
b = K2/I12 + K1/I23;
c = K1*K2/I123;

%Polinômio característico:
pC = -lmb*(lmb^2 - b*lmb + c);
%{
pCref = simplify(det(M\K - lmb*eye(3)));
erro = simplify(pCref - pC);
ad = 1;
%ok!
%}

%Raízes de pC:
delta = (K2/I12 - K1/I23)^2 + 4*K1*K2/Ibp2^2;
%{
deltaRef = b^2 - 4*c;
erro = simplify(deltaRef - delta);
ad = 1;
%ok!
%}
%Obs.: 
%(1) delta é maior que zero
%(2) b^2 = delta + 4*K1*K2/I123 > delta
lmb1 = 0;
lmb2 = (1/2)*(b - sqrt(delta));
lmb3 = (1/2)*(b + sqrt(delta));
%{
pClmb2 = simplify(subs(pC,lmb,lmb2));
pClmb3 = simplify(subs(pC,lmb,lmb3));
ad = 1;
%ok!
%}

%Modos de vibração:
invMK = [ K1/Ibp3,       -K1/Ibp3,        0; 
         -K1/Ibp2, (K1 + K2)/Ibp2, -K2/Ibp2; 
                0,        -K2/Ip1,   K2/Ip1];
%{
invMKref = M\K;
erro = simplify(invMKref - invMK);
ad = 1;
%ok!
%}
iMKmLmbId = [ K1/Ibp3 - lmb,             -K1/Ibp3,              0; 
                   -K1/Ibp2, (K1 + K2)/Ibp2 - lmb,       -K2/Ibp2; 
                          0,              -K2/Ip1,   K2/Ip1 - lmb];

iMKmLmbId1 = [ 1,             -K1/Ibp3/(K1/Ibp3 - lmb),            0; 
              -1,     ((K1 + K2)/Ibp2 - lmb)/(K1/Ibp2),       -K2/K1; 
               0,                              -K2/Ip1, K2/Ip1 - lmb];
%{
iMKmLmbId1Ref = diag([(K1/Ibp3 - lmb)^-1, (K1/Ibp2)^-1, 1])*iMKmLmbId;
erro = simplify(iMKmLmbId1Ref - iMKmLmbId1);
ad = 1;
%ok!
%}
iMKmLmbId2 = [1,                                                             -(K1/Ibp3)/(K1/Ibp3 - lmb),            0; 
              0, (Ibp2/K1)*(K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2)/(K1/Ibp3 - lmb),       -K2/K1; 
              0,                                                                                -K2/Ip1, K2/Ip1 - lmb];
%{
iMKmLmbId2ref = [1, 0, 0; 1, 1, 0; 0, 0, 1]*iMKmLmbId1;
erro = simplify(iMKmLmbId2ref - iMKmLmbId2);
ad = 1;
%ok!
%}
iMKmLmbId3 = [1,                                                   -(K1/Ibp3)/(K1/Ibp3 - lmb),            0; 
              0, (K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2)/(K1/Ibp3 - lmb),     -K2/Ibp2; 
              0,                                                                      -K2/Ip1, K2/Ip1 - lmb];
%{
iMKmLmbId3ref = [1, 0, 0; 0, K1/Ibp2, 0; 0, 0, 1]*iMKmLmbId2;
erro = simplify(iMKmLmbId3ref - iMKmLmbId3);
ad = 1;
%ok!
%}
iMKmLmbId4 = [1, -(K1/Ibp3)/(K1/Ibp3 - lmb),                                                                                       0; 
              0,                          1, -(K2/Ibp2)*(K1/Ibp3 - lmb)/(K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2); 
              0,                    -K2/Ip1,                                                                            K2/Ip1 - lmb];
%{
iMKmLmbId4ref = [1, 0, 0; 0, (K1/Ibp3 - lmb)/(K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2), 0; 0, 0, 1]*iMKmLmbId3;
erro = simplify(iMKmLmbId4ref - iMKmLmbId4);
ad = 1;
%ok!
%}          
iMKmLmbId5 = [1, 0,                                -(K2/Ibp2)*(K1/Ibp3)/(K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2); 
              0, 1,                          -(K2/Ibp2)*(K1/Ibp3 - lmb)/(K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2); 
              0, 0, ((K2/Ip1 - lmb)*(K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2) - (K2/Ip1)*(K2/Ibp2)*(K1/Ibp3 - lmb))/(K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2)];         
%{
iMKmLmbId5ref = [1, (K1/Ibp3)/(K1/Ibp3 - lmb), 0; 0, 1, 0; 0, K2/Ip1, 1]*iMKmLmbId4;
erro = simplify(iMKmLmbId5ref - iMKmLmbId5);
ad = 1;
%ok!
%}
dC = K2*K1/(Ibp3*Ibp2) - lmb*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb^2;
iMKmLmbId6 = [1, 0,       -(K2/Ibp2)*(K1/Ibp3)/dC; 
              0, 1, -(K2/Ibp2)*(K1/Ibp3 - lmb)/dC; 
              0, 0,                         pC/dC];         
%{
erro = simplify(iMKmLmbId6 - iMKmLmbId5);
ad = 1;
%ok!
%}
%{
dClmb2 = simplify(subs(dC,lmb,lmb2));
dClmb3 = simplify(subs(dC,lmb,lmb3));
ad = 1;
%ok!
%Obs.: dClmb2 e dClmb3 são não nulos
%}
dC2 = K2*K1/(Ibp3*Ibp2) - lmb2*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb2^2;
iMKmLmbId6lmb2 = [1, 0,        -(K2/Ibp2)*(K1/Ibp3)/dC2; 
                  0, 1, -(K2/Ibp2)*(K1/Ibp3 - lmb2)/dC2; 
                  0, 0,                               0];
dC3 = K2*K1/(Ibp3*Ibp2) - lmb3*(K1/Ibp3 + (K1 + K2)/Ibp2) + lmb3^2;
iMKmLmbId6lmb3 = [1, 0,        -(K2/Ibp2)*(K1/Ibp3)/dC3; 
                  0, 1, -(K2/Ibp2)*(K1/Ibp3 - lmb3)/dC3; 
                  0, 0,                               0];
              

              
              
              
          
          
