clear
clc

%%
%Expressões e parâmetros com variáveis reais
syms I1 I2b M1 M2b dm l2 K1 K2 k1 k2 c2 c3 c3cc b2 b3 betaS k3 g
%Obs.:
%(1) I2b: momento de inércia do rotor 2, maciço e completo
%(2) M2 + m2 = M2b + dm

Mt = diag([I2b + l2^2*dm, I1]);
Kt = [K1+K2, -K2; -K2, K2];
Ct = diag([c2, c3]);
Bt = [K1; 0];

Mu = diag([M2b + dm, M1]);
Ku = [k1+k2, -k2; -k2, k2];
Cu = diag([b2, b3]);
Bu = [k1/betaS; 0];
Wu = Mu*[1; 1]*g;
uEq = Ku\Wu;
%Variável longitudinal utilizada nas equações de estado: ub = u - uEq

M = [Mt, zeros(2); zeros(2), Mu];
K = [Kt, zeros(2); zeros(2), Ku];
C = [Ct, zeros(2); zeros(2), Cu];
B4 = [Bt, zeros(2,1); zeros(2,1), Bu];

Ab = [zeros(4), eye(4); -M\K, -M\C];
B = [zeros(4,2); M\B4];
B1 = B(:,1);
B2 = B(:,2);

Id8 = eye(8);
C1 = Id8(6,:);
C2 = Id8(4,:);
%C2_u = Id8(4,:);
%C2_up = Id8(8,:);

%Termos não lineares:
syms teta2 teta3 u2b u3b teta2p teta3p u2bp u3bp
syms Tat(teta3p,u3b)
Z = [teta2; teta3; u2b; u3b; teta2p; teta3p; u2bp; u3bp];
Tc = [0; Tat(teta3p,u3b); 0; k3*u3b];
F = [zeros(4,1); M\Tc];

%Forma Normal
%(1) y1 = w3
beta0 = C1;
beta1 = beta0*Ab;
beta2 = beta1*Ab;
beta3 = beta2*Ab;

beta0B1 = beta0*B1;
beta0B2 = beta0*B2;
beta1B1 = beta1*B1;
beta1B2 = beta1*B2;
beta2B1 = beta2*B1;
beta2B2 = beta2*B2;
%{
beta3B1 = beta3*B1;
beta3B2 = beta3*B2;
%}

%(2) y2 = N(u3b) = k3*u3b; u3b = u3 - uC, sendo uC = uEq(2,1)
alfa0 = k3*C2;
%alfa0 = k3*C2_u + c3cc*C2_up;
alfa1 = alfa0*Ab;
alfa2 = alfa1*Ab;
alfa3 = alfa2*Ab;
alfa4 = alfa3*Ab;

alfa0B1 = alfa0*B1;
alfa0B2 = alfa0*B2;
alfa1B1 = alfa1*B1;
alfa1B2 = alfa1*B2;
alfa2B1 = alfa2*B1;
alfa2B2 = alfa2*B2;
alfa3B1 = alfa3*B1;
alfa3B2 = alfa3*B2;

%Saída y1 = h1(Z):
h1 = C1*Z;

Qsi0 = gradient(beta0*F,Z);
Qsi1 = gradient(beta1*F + Qsi0.'*(Ab*Z - F),Z);

%y1p:
Lf_h1 = beta1*Z - beta0*F;
Lg1_h1 = beta0B1;
Lg2_h1 = beta0B2;

%y1pp:
Lf_Lf_h1 = beta2*Z - beta1*F - Qsi0.'*(Ab*Z - F);
Lg1_Lf_h1 = beta1B1 - Qsi0.'*B1;
Lg2_Lf_h1 = beta1B2 - Qsi0.'*B2;

%y1ppp:
Lf_L2f_h1 = beta3*Z - beta2*F - Qsi1.'*(Ab*Z - F);
Lg1_L2f_h1 = beta2B1 - Qsi1.'*B1;
Lg2_L2f_h1 = beta2B2 - Qsi1.'*B2;

%Saída y2 = h2(Z):
h2 = k3*C2*Z;

Fi0 = gradient(alfa0*F,Z);
Fi1 = gradient(alfa1*F + Fi0.'*(Ab*Z - F),Z);
Fi2 = gradient(alfa2*F + Fi1.'*(Ab*Z - F),Z);

%y2p:
Lf_h2 = alfa1*Z - alfa0*F;
Lg1_h2 = alfa0B1;
Lg2_h2 = alfa0B2;

%y2pp:
Lf_Lf_h2 = alfa2*Z - alfa1*F - Fi0.'*(Ab*Z - F);
Lg1_Lf_h2 = alfa1B1 - Fi0.'*B1;
Lg2_Lf_h2 = alfa1B2 - Fi0.'*B2;

%y2ppp:
Lf_L2f_h2 = alfa3*Z - alfa2*F - Fi1.'*(Ab*Z - F);
Lg1_L2f_h2 = alfa2B1 - Fi1.'*B1;
Lg2_L2f_h2 = alfa2B2 - Fi1.'*B2;

%y2pppp:
Lf_L3f_h2 = alfa4*Z - alfa3*F - Fi2.'*(Ab*Z - F);
Lg1_L3f_h2 = alfa3B1 - Fi2.'*B1;
Lg2_L3f_h2 = alfa3B2 - Fi2.'*B2;

%%
%Dinâmica interna:
syms q1 q2 q3 q4 q5 q6 q7 q8
syms q epsilon ne
%q = epsilon;
%q1 = 1;
%q3 = 0; q4 = 0; q6 = 0; q8 = 0;
%A_mi = [q1 q-q1 q3 q4 epsilon q6 0 q8];
A_mi = [q1 q-q1 q3 q4 q5 q6 q7 q8];
%A_mi = [q1 q-q1 q3 q4 q5 q6 0 q8];
%A_mi = beta3;
%A_mi(1,5) = -epsilon;
%{
teste1 = A_mi*B;
teste2 = A_mi*Ab*B;
%ok!
%}
%sigma = [beta0; beta1; beta2; alfa0; alfa1; alfa2; alfa3; A_mi ];
sigma = [beta0; beta1; beta2; alfa0; alfa1; alfa2; alfa3; A_mi ];
det_sigma = det(sigma);
inv_sigma = simplify(sigma\Id8);
Ad = Id8(:,8);
Bd = Id8(:,1:7);
Pi_mi = simplify(A_mi*Ab*inv_sigma);
Gama_mi = Pi_mi*Ad;
Sigma_mi = Pi_mi*Bd;
FiDif = [0; 
    beta0*F; 
    beta1*F + Qsi0.'*(Ab*Z - F); 
    0; 
    alfa0*F; 
    alfa1*F + Fi0.'*(Ab*Z - F); 
    alfa2*F + Fi1.'*(Ab*Z - F); 
    0];
dFiDif_dZT = [zeros(1,8); 
    Qsi0.'; 
    Qsi1.'; 
    zeros(1,8); 
    Fi0.'; 
    Fi1.'; 
    Fi2.'; 
    zeros(1,8)];
sigma_atr = sigma - dFiDif_dZT;
det_sigma_atr = det(sigma_atr);
inv_sigma_atr = simplify(sigma_atr\Id8);
Pi_mi_atr = simplify(A_mi*Ab*inv_sigma_atr);
Gama_mi_atr = simplify(Pi_mi_atr*Ad);
Sigma_mi_atr = simplify(Pi_mi_atr*Bd);
%{
syms alfa
teste3 = A_mi*(Id8 + alfa*Ab)*B;
Pi = simplify(A_mi*(Id + alfa*Ab)*Ab*inv_sigma_atr);
Gama = simplify(Pi*Ad);
Sigma = simplify(Pi*Bd);
syms p1 p2 p3 p4 p5 p6 p7 p8
syms p
B_mi = [p1 p-p1 p3 p4 p5 p6 p7 p8];
teste4 = A_mi*B;
teste5 = (B_mi + A_mi*Ab)*B;
Pi = simplify((B_mi + A_mi*Ab)*Ab*inv_sigma);
Gama = simplify(Pi*Ad);
Sigma = simplify(Pi*Bd);
syms r1 r2 r3 r4 r5 r6 r7 r8
syms r
C_mi = [r1 r-r1 r3 r4 r5 r6 r7 r8];
teste6 = A_mi*B;
%teste6subs = subs(teste6,[q5,q7],[0,0]);
teste7 = (B_mi + A_mi*Ab)*B;
%teste7subs = simplify(subs(teste7,[q5,q7,p5,p7],[0,0,-q1,-q3]));
teste8 = (C_mi + B_mi*Ab + A_mi*Ab^2)*B;
Pi = (C_mi + B_mi*Ab + A_mi*Ab^2)*Ab*inv_sigma;
GamaRef = simplify(Pi*Ad);
Sigma = Pi*Bd;
%}
%{
ad = I1*I2b^2*p1 + I1*I2b^2*r5 + I2b^2*K2*q6 + I1*c2^2*q5 + I1*dm^2*l2^4*p1 + I1*dm^2*l2^4*r5 + K2*dm^2*l2^4*q6 - I1*I2b*K1*q5 - I1*I2b*K2*q5 - I1*I2b*c2*p5 - I1*I2b*c2*q1 + 2*I1*I2b*dm*l2^2*p1 + 2*I1*I2b*dm*l2^2*r5 - I1*K1*dm*l2^2*q5 - I1*K2*dm*l2^2*q5 + 2*I2b*K2*dm*l2^2*q6 - I1*c2*dm*l2^2*p5 - I1*c2*dm*l2^2*q1;
ad_p1 = p1*I1*(I2b + dm*l2^2)^2;
ad_r5 = r5*I1*(I2b + dm*l2^2)^2;
ad_q5 = q5*I1*(c2^2 - (K1 + K2)*(I2b + dm*l2^2));
ad_q6 = q6*K2*(I2b + dm*l2^2)^2;
ad_p5 = -p5*I1*c2*(I2b + dm*l2^2);
ad_q1 = -q1*I1*c2*(I2b + dm*l2^2);
%erro = simplify(ad - (ad_p1 + ad_r5 + ad_q5 + ad_q6 + ad_p5 + ad_q1));
%}
%Gama = -(K1/(q*(dm*l2^2 + I2b)))*...
%    (p1 + r5 - (q5/(I2b + dm*l2^2))*(K1 + K2 - c2^2/(I2b + dm*l2^2)) + (q6/I1)*K2 - (p5 + q1)*c2/(I2b + dm*l2^2));
%erro = simplify(GamaRef - Gama);
%%
%Cálculo das expressões de incerteza:
%syms teta2 teta3 teta2p teta3p
syms teta2 teta3 w2 w3
Z = [teta2; teta3; w2; w3];

syms d rho E nu L1b Dr1 L2b Dr2 
syms rhom2 Dm l2
syms mi0 mi1 mi2 beta1mi beta2mi 
syms mic0 w0 nmi lambdami
syms k3 R3
syms L L1
syms alfa11 alfa12 alfa21 alfa22
syms pf

syms dd drho dE dnu dL1b dDr1 dL2b dDr2 
syms drhom2 dDm dl2
syms dmi0 dmi1 dmi2 dbeta1mi dbeta2mi
syms dmic0 dw0 dnmi dlambdami
syms dk3 dR3
syms dL dL1
syms dalfa11 dalfa12 dalfa21 dalfa22
syms dpf

%{
G_hat = E_hat/(2*(1 + nu_hat));
J_hat = pi*d_hat^4/32;
A_hat = pi*d_hat^2/4;
K1_hat = alfa111_hat*G_hat*J_hat/L1_hat;
K2_hat = alfa112_hat*G_hat*J_hat/(L_hat - L1_hat);
k1_hat = alfa221_hat*E_hat*A_hat/L1_hat;
k2_hat = alfa222_hat*E_hat*A_hat/(L_hat - L1_hat);
M1_hat = rho_hat*(pi*Dr1_hat^2/4)*L1b_hat;
Ar2_hat = pi*Dr2_hat^2/4;
Am2_hat = pi*Dm_hat^2/4;
M2_hat = rho_hat*(Ar2_hat - Am2_hat)*L2b_hat;
m2_hat = rho_m2_hat*Am2_hat*L2b_hat;
M2b_hat = M2_hat + m2_hat;
rCG_hat = l2_hat*((m2_hat - (Dm_hat/Dr2_hat)^1*M2_hat)/(m2_hat + M2_hat));
Jp1_hat = pi*Dr1_hat^4/32;
I1_hat = rho_hat*Jp1_hat*L1b_hat;
Jp2_hat = pi*Dr2_hat^4/32;
Jpm2_hat = pi*Dm_hat^4/32 + Am2_hat*l2_hat^2;
I2_hat = rho_hat*(Jp2_hat - Jpm2_hat)*L2b_hat;
Im2_hat = rho_m2_hat*Jpm2_hat*L2b_hat;
I2b_hat = I2_hat + Im2_hat;
beta_s_hat = 2*pi/pf_hat;
%}

syms G(E,nu)
syms A(d) Ar1(Dr1) Ar2(Dr2) Am2(Dm)
syms J(d) Jp1(Dr1) Jp2(Dr2) Jpm20(Dm)
syms betaS(pf)

%G = E/(2*(1 + nu));%
%J = pi*d^4/32;%
%A = pi*d^2/4;%
K1 = alfa11*G(E,nu)*J(d)/L1;%
K2 = alfa12*G(E,nu)*J(d)/(L - L1);%
k1 = alfa21*E*A(d)/L1;%
k2 = alfa22*E*A(d)/(L - L1);%
%Ar1 = pi*Dr1^2/4;%
M1 = rho*Ar1(Dr1)*L1b;%
%Ar2 = pi*Dr2^2/4;%
%Am2 = pi*Dm^2/4;%
M2b = rho*Ar2(Dr2)*L2b; 
dm = (rhom2 - rho)*Am2(Dm)*L2b;
%Obs.:
%(1) M2 = rho*(Ar2(Dr2) - Am2(Dm))*L2b;
%(2) m2 = rhom2*Am2(Dm)*L2b;
%(3) M2b + dm = M2 + m2;

%Excentricidade final do rotor 2 (rCG):
betam = Dm/Dr2;
alfam = (rhom2 - rho)/rho;
rCG = betam^2*l2*((betam^2 + alfam)/(1 + alfam*betam^2));

%Jp1 = pi*Dr1^4/32;%
I1 = rho*Jp1(Dr1)*L1b;
%Jp2 = pi*Dr2^4/32;
%Jpm20 = pi*Dm^4/32;
%Jpm2 = Jpm20(Dm) + Am2(Dm)*l2^2;
I2b = rho*Jp2(Dr2)*L2b;
%Obs.:
%(1) I2b: momento de inércia do rotor maciço completo

Mt = diag([I2b + l2^2*dm, I1]);
Kt = [K1 + K2, -K2; -K2, K2];
Ct = diag([c2, c3]);
Bt = [K1; 0];

%{
Mu = diag([M2b + dm, M1]);
Ku = [k1 + k2, -k2; -k2, k2];
Cu = diag([b2, b3]);
Bu = [k1/betaS(pf); 0];
Wu = Mu*[1; 1]*g;
%}

M = Mt;
K = Kt;
C = Ct;
B2 = Bt;

Ab = [zeros(2), eye(2); -M\K, -M\C];
B = [zeros(2,1); M\B2];

%%
syms mi3c(mic0,w3,w0,lambdami,nmi)
syms mi3dc(mi0,w3,mi1,mi2,beta1mi,beta2mi)
syms Normalb0(w3)
Tatr3c(mic0,w3,w0,lambdami,nmi,R3) = ...
    mi3c(mic0,w3,w0,lambdami,nmi)*...
    Normalb0(w3)*...
    R3;
Tatr3dc(mi0,w3,mi1,mi2,beta1mi,beta2mi,R3) = ...
    mi3dc(mi0,w3,mi1,mi2,beta1mi,beta2mi)*...
    Normalb0(w3)*...
    R3;
Tc = [0; Tatr3c(mic0,w3,w0,lambdami,nmi,R3)];
Tdc = [0; Tatr3dc(mi0,w3,mi1,mi2,beta1mi,beta2mi,R3)];
Fc = [zeros(2,1); M\Tc];
Fdc = [zeros(2,1); M\Tdc];

Id4 = eye(4);
C1 = Id4(4,:);
beta0 = C1;
beta1 = beta0*Ab;
beta2 = beta1*Ab;
beta3 = beta2*Ab;

beta0B = beta0*B;
beta1B = beta1*B;
beta2B = beta2*B;

Qsi0c = gradient(beta0*Fc,Z);
Qsi1c = gradient(beta1*Fc + Qsi0c.'*(Ab*Z - Fc),Z);
Qsi0dc = gradient(beta0*Fdc,Z);
Qsi1dc = gradient(beta1*Fdc + Qsi0dc.'*(Ab*Z - Fdc),Z);
%Qsi2c_hat = gradient(beta2_hat*Fc_hat + Qsi1c_hat.'*(Ab_hat*Z - Fc_hat),Z);
%{
diary ('Qsi2_hat.txt');
Qsi2_hat % your symbolic variable
diary off
%}


a1c = beta3*Z - beta2*Fc - Qsi1c.'*(Ab*Z - Fc);
a1dc = beta3*Z - beta2*Fdc - Qsi1dc.'*(Ab*Z - Fdc);
b = beta2B - Qsi1c.'*B;

qsiTil = [d; rho; E; nu; 
    L1b; Dr1; L2b; Dr2; 
    rhom2; Dm; l2;
    mi0; mi1; mi2; beta1mi; beta2mi; 
    mic0; w0; lambdami; 
    R3; 
    L; L1; 
    alfa11; alfa12; alfa21; alfa22;...
    pf];
dqsiTil = [dd; drho; dE; dnu;...
    dL1b; dDr1; dL2b; dDr2;...
    drhom2; dDm; dl2;...
    dmi0; dmi1; dmi2; dbeta1mi; dbeta2mi;...
    dmic0; dw0; dlambdami;...
    dR3;...
    dL; dL1;...
    dalfa11; dalfa12; dalfa21; dalfa22;...
    dpf];

dbeta3dqsi = jacobian(beta3,qsiTil).';
dbeta2Fcdqsi = gradient(beta2*Fc,qsiTil);
dQsi1cAZdqsi = gradient(Qsi1c.'*Ab*Z,qsiTil);
dQsi1cFcdqsi = gradient(Qsi1c.'*Fc,qsiTil);
da1cdqsiRef = gradient(a1c,qsiTil);
da1cdqsi = dbeta3dqsi*Z - dbeta2Fcdqsi - dQsi1cAZdqsi + dQsi1cFcdqsi;

d2beta2FcdqsidZT = jacobian(dbeta2Fcdqsi,Z);
d2Qsi1cAZdqsidZT = jacobian(dQsi1cAZdqsi,Z);
d2Qsi1FcdqsidZT = jacobian(dQsi1cFcdqsi,Z);
d2a1cdqsidZTref = jacobian(da1cdqsiRef,Z);
d2a1cdqsidZT = dbeta3dqsi - d2beta2FcdqsidZT - d2Qsi1cAZdqsidZT + d2Qsi1FcdqsidZT;
%{
erro1 = simplify(da1cdqsiRef - da1cdqsi);
ad = 1;
%ok!
%}
%{
erro1 = simplify(d2a1cdqsidZTref - d2a1cdqsidZT);
%ok!
ad = 1;
%}
%{
diary ('da1_hat_dqsi_ref_teste.txt');
da1_hat_dqsi_ref % your symbolic variable
diary off
%}
%{
diary ('d2Qsi1cAZdqsidZT.txt');
d2Qsi1cAZdqsidZT % your symbolic variable
diary off
%}
%{
diary ('d2Qsi1FcdqsidZT.txt');
d2Qsi1FcdqsidZT % your symbolic variable
diary off
%}
dbeta2Fdcdqsi = gradient(beta2*Fdc,qsiTil);
dQsi1dcAZdqsi = gradient(Qsi1dc.'*Ab*Z,qsiTil);
dQsi1dcFdcdqsi = gradient(Qsi1dc.'*Fdc,qsiTil);
da1dcdqsiRef = gradient(a1dc,qsiTil);
da1dcdqsi = dbeta3dqsi*Z - dbeta2Fdcdqsi - dQsi1dcAZdqsi + dQsi1dcFdcdqsi;

d2beta2FdcdqsidZT = jacobian(dbeta2Fdcdqsi,Z);
d2Qsi1dcAZdqsidZT = jacobian(dQsi1dcAZdqsi,Z);
d2Qsi1FdcdqsidZT = jacobian(dQsi1dcFdcdqsi,Z);
d2a1dcdqsidZTref = jacobian(da1dcdqsiRef,Z);
d2a1dcdqsidZT = dbeta3dqsi - d2beta2FdcdqsidZT - d2Qsi1dcAZdqsidZT + d2Qsi1FdcdqsidZT;
%{
erro1 = simplify(da1dcdqsiRef - da1dcdqsi);
ad = 1;
%ok!
%}
%{
erro1 = simplify(d2a1dcdqsidZTref - d2a1dcdqsidZT);
%ok!
ad = 1;
%}
%{
diary ('d2Qsi1dcAZdqsidZT.txt');
d2Qsi1dcAZdqsidZT % your symbolic variable
diary off
%}
%{
diary ('d2Qsi1FdcdqsidZT.txt');
d2Qsi1FdcdqsidZT % your symbolic variable
diary off
%}
dbdqsi = gradient(b,qsiTil);

da1c = da1cdqsi.'*dqsiTil;
da1dc = da1dcdqsi.'*dqsiTil;
db = dbdqsi.'*dqsiTil;
F1c = abs(da1c);
F1dc = abs(da1dc);
Dlc = (1/b)*abs(db);

%Derivadas:
dQsi0cTdZT = jacobian(Qsi0c,Z);
dQsi1cTdZT = jacobian(Qsi1c,Z);
dQsi0dcTdZT = jacobian(Qsi0dc,Z);
dQsi1dcTdZT = jacobian(Qsi1dc,Z);

z1 = teta2; z2 = teta3; z3 = w2; z4 = w3;
dQsi0cdz1 = diff(Qsi0c,z1).'; %zero(1,4)
dQsi0cdz2 = diff(Qsi0c,z2).'; %zero(1,4)
dQsi0cdz3 = diff(Qsi0c,z3).'; %zero(1,4)
dQsi0cdz4 = diff(Qsi0c,z4).';
dQsi1cdz1 = diff(Qsi1c,z1).';
dQsi1cdz2 = diff(Qsi1c,z2).';
dQsi1cdz3 = diff(Qsi1c,z3).'; %zero(1,4)
dQsi1cdz4 = diff(Qsi1c,z4).';

d2dFiDifcdz1dZT = [zeros(1,4); dQsi0cdz1; dQsi1cdz1; zeros(1,4)];
d2dFiDifcdz2dZT = [zeros(1,4); dQsi0cdz2; dQsi1cdz2; zeros(1,4)];
d2dFiDifcdz3dZT = [zeros(1,4); dQsi0cdz3; dQsi1cdz3; zeros(1,4)]; %zero(4,4)
d2dFiDifcdz4dZT = [zeros(1,4); dQsi0cdz4; dQsi1cdz4; zeros(1,4)];

dQsi0dcdz1 = diff(Qsi0dc,z1).'; %zero(1,4)
dQsi0dcdz2 = diff(Qsi0dc,z2).'; %zero(1,4)
dQsi0dcdz3 = diff(Qsi0dc,z3).'; %zero(1,4)
dQsi0dcdz4 = diff(Qsi0dc,z4).';
dQsi1dcdz1 = diff(Qsi1dc,z1).';
dQsi1dcdz2 = diff(Qsi1dc,z2).';
dQsi1dcdz3 = diff(Qsi1dc,z3).'; %zero(1,4)
dQsi1dcdz4 = diff(Qsi1dc,z4).';

d2dFiDifdcdz1dZT = [zeros(1,4); dQsi0dcdz1; dQsi1dcdz1; zeros(1,4)];
d2dFiDifdcdz2dZT = [zeros(1,4); dQsi0dcdz2; dQsi1dcdz2; zeros(1,4)];
d2dFiDifdcdz3dZT = [zeros(1,4); dQsi0dcdz3; dQsi1dcdz3; zeros(1,4)]; %zero(4,4)
d2dFiDifdcdz4dZT = [zeros(1,4); dQsi0dcdz4; dQsi1dcdz4; zeros(1,4)];

dmi3cdw3 = diff(mi3c,w3);
dNormalb0dw3(w3) = diff(Normalb0,w3);
dTatr3cdw3(mic0,w3,w0,lambdami,nmi,R3) = (dmi3cdw3(mic0,w3,w0,lambdami,nmi)*Normalb0(w3) + mi3c(mic0,w3,w0,lambdami,nmi)*dNormalb0dw3(w3))*R3;
dTatr3cduTilbpT = [0, dTatr3cdw3(mic0,w3,w0,lambdami,nmi,R3)];
dTcduTilbpT = [zeros(1,2); dTatr3cduTilbpT];
duTilbpdZT = [zeros(2), eye(2)];
dTcdZT = dTcduTilbpT*duTilbpdZT;
dFcdZT = [zeros(2,4); M\dTcdZT];
%{
dFcdZTRef = jacobian(Fc,Z);
erro = simplify(dFcdZTRef - dFcdZT);
ad = 1;
%}

dFdcdZTRef = jacobian(Fdc,Z);

%%
da1cdZref = gradient(a1c,Z);
d2a1cdqsiTdZref2 = jacobian(da1cdZref,qsiTil);

dbeta3dqsiT = dbeta3dqsi.';
d2beta2FcdqsiTdZ = jacobian(dbeta2Fcdqsi,Z).';
d2Qsi1cAZdqsiTdZ = jacobian(dQsi1cAZdqsi,Z).';
d2Qsi1FcdqsiTdZ = jacobian(dQsi1cFcdqsi,Z).';
d2a1cdqsiTdZref = jacobian(da1cdqsiRef,Z).';
%erro = simplify(d2a1cdqsiTdZref2 - d2a1cdqsiTdZref);
d2a1cdqsiTdZ = dbeta3dqsiT - d2beta2FcdqsiTdZ - d2Qsi1cAZdqsiTdZ + d2Qsi1FcdqsiTdZ;
%erro = simplify(d2a1cdqsiTdZref - d2a1cdqsiTdZ); %OK! AD,16/11/2021

ddZTdbeta3dqsiTdqsi = jacobian(dbeta3dqsiT*dqsiTil,Z);
ddZTd2beta2FcdqsiTdZdqsi = jacobian(d2beta2FcdqsiTdZ*dqsiTil,Z);
%ddZTd2beta2FcdqsiTdZdqsi = simplify(ddZTd2beta2FcdqsiTdZdqsi, 'IgnoreAnalyticConstraints', true);
%{
diary ('ddZTd2b2FcdqTdZdq.txt');
ddZTd2beta2FcdqsiTdZdqsi % your symbolic variable
diary off
%}
ddZTd2Qsi1cAZdqsiTdZdqsi = jacobian(d2Qsi1cAZdqsiTdZ*dqsiTil,Z);
%ddZTd2Qsi1cAZdqsiTdZdqsi = simplify(ddZTd2Qsi1cAZdqsiTdZdqsi, 'IgnoreAnalyticConstraints', true);
%{
diary ('ddZTd2Q1cAZdqTdZdq.txt');
ddZTd2Qsi1cAZdqsiTdZdqsi % your symbolic variable
diary off
%}
ddZTd2Qsi1FcdqsiTdZdqsi = jacobian(d2Qsi1FcdqsiTdZ*dqsiTil,Z);
%ddZTd2Qsi1FcdqsiTdZdqsi = simplify(ddZTd2Qsi1FcdqsiTdZdqsi, 'IgnoreAnalyticConstraints', true);
%{
diary ('ddZTd2Q1FcdqTdZdq.txt');
ddZTd2Qsi1FcdqsiTdZdqsi % your symbolic variable
diary off
%}
ddZTd2a1cdqsiTdZdqsi = ddZTdbeta3dqsiTdqsi - ddZTd2beta2FcdqsiTdZdqsi - ddZTd2Qsi1cAZdqsiTdZdqsi + ddZTd2Qsi1FcdqsiTdZdqsi;
ddZTd2a1cdqsiTdZdqsiref = jacobian(d2a1cdqsiTdZ*dqsiTil,Z);
%{
erroC = simplify(ddZTd2a1cdqsiTdZdqsiref - ddZTd2a1cdqsiTdZdqsi, 'IgnoreAnalyticConstraints', true);
ad = 1;
%ok!
%}

d2beta2FdcdqsiTdZ = jacobian(dbeta2Fdcdqsi,Z).';
d2Qsi1dcAZdqsiTdZ = jacobian(dQsi1dcAZdqsi,Z).';
d2Qsi1FdcdqsiTdZ = jacobian(dQsi1dcFdcdqsi,Z).';
d2a1dcdqsiTdZref = jacobian(da1dcdqsiRef,Z).';
d2a1dcdqsiTdZ = dbeta3dqsiT - d2beta2FdcdqsiTdZ - d2Qsi1dcAZdqsiTdZ + d2Qsi1FdcdqsiTdZ;
%erro = simplify(d2a1dcdqsiTdZref - d2a1dcdqsiTdZ); %OK! AD,16/11/2021

ddZTd2beta2FdcdqsiTdZdqsi = jacobian(d2beta2FdcdqsiTdZ*dqsiTil,Z);
%ddZTd2beta2FdcdqsiTdZdqsi = simplify(ddZTd2beta2FdcdqsiTdZdqsi, 'IgnoreAnalyticConstraints', true);
%{
diary ('ddZTd2b2Fdc.txt');
ddZTd2beta2FdcdqsiTdZdqsi % your symbolic variable
diary off
%}
ddZTd2Qsi1dcAZdqsiTdZdqsi = jacobian(d2Qsi1dcAZdqsiTdZ*dqsiTil,Z);
%ddZTd2Qsi1dcAZdqsiTdZdqsi = simplify(ddZTd2Qsi1dcAZdqsiTdZdqsi, 'IgnoreAnalyticConstraints', true);
%{
diary ('ddZTd2Q1dcAZ.txt');
ddZTd2Qsi1dcAZdqsiTdZdqsi % your symbolic variable
diary off
%}
ddZTd2Qsi1FdcdqsiTdZdqsi = jacobian(d2Qsi1FdcdqsiTdZ*dqsiTil,Z);
%ddZTd2Qsi1FdcdqsiTdZdqsi = simplify(ddZTd2Qsi1FdcdqsiTdZdqsi, 'IgnoreAnalyticConstraints', true);
%{
diary ('ddZTd2Q1Fdc.txt');
ddZTd2Qsi1FdcdqsiTdZdqsi % your symbolic variable
diary off
%}
ddZTd2a1dcdqsiTdZdqsi = ddZTdbeta3dqsiTdqsi - ddZTd2beta2FdcdqsiTdZdqsi - ddZTd2Qsi1dcAZdqsiTdZdqsi + ddZTd2Qsi1FdcdqsiTdZdqsi;
ddZTd2a1dcdqsiTdZdqsiref = jacobian(d2a1dcdqsiTdZ*dqsiTil,Z);
%{
erroDC = simplify(ddZTd2a1dcdqsiTdZdqsiref - ddZTd2a1dcdqsiTdZdqsi, 'IgnoreAnalyticConstraints', true);
ad = 1;
%ok!
%}

%%
syms R3
syms Normalb0(w3) mi3cHat(w3) mi3dcHat(w3) dmi3cHatdw3(w3) dmi3dcHatdw3(w3) dNormalb0dw3(w3)

dTatr3cdw3 = (dmi3cHatdw3*Normalb0 + mi3cHat*dNormalb0dw3)*R3;
dTatr3cduTilbpT = [0, dTatr3cdw3];
dTc0HatduTilbpT = [zeros(1,2); dTatr3cduTilbpT];
Id4 = eye(4);
Av = Id4(:,3:4);
duTilbpdZT = Av.';
dTc0HatdZT = dTc0HatduTilbpT*duTilbpdZT;
dFcHatdZT = [zeros(2,4); Mt\dTc0HatdZT];
beta2dFcHatdZT = beta2*dFcHatdZT;
ddZTdFcTdZb2T = jacobian(beta2dFcHatdZT.',Z);

dTatr3dcdw3 = (dmi3dcHatdw3*Normalb0 + mi3dcHat*dNormalb0dw3)*R3;
dTatr3dcduTilbpT = [0, dTatr3dcdw3];
dTdc0HatduTilbpT = [zeros(1,2); dTatr3dcduTilbpT];
dTdc0HatdZT = dTdc0HatduTilbpT*duTilbpdZT;
dFdcHatdZT = [zeros(2,4); Mt\dTdc0HatdZT];
beta2dFdcHatdZT = beta2*dFdcHatdZT;
ddZTdFdcTdZb2T = jacobian(beta2dFdcHatdZT.',Z);

Qsi1cAZ = Qsi1c.'*(Ab*Z);
dQsi1cAZdZT = gradient(Qsi1cAZ,Z).';
d2Qsi1cAZdZTdZ = jacobian(dQsi1cAZdZT.',Z);

Qsi1dcAZ = Qsi1dc.'*(Ab*Z);
dQsi1dcAZdZT = gradient(Qsi1dcAZ,Z).';
d2Qsi1dcAZdZTdZ = jacobian(dQsi1dcAZdZT.',Z);

Qsi1cFc = Qsi1c.'*Fc;
dQsi1cFcdZT = gradient(Qsi1cFc,Z).';
d2Qsi1cFcdZTdZ = jacobian(dQsi1cFcdZT.',Z);

Qsi1dcFdc = Qsi1dc.'*Fdc;
dQsi1dcFdcdZT = gradient(Qsi1dcFdc,Z).';
d2Qsi1dcFdcdZTdZ = jacobian(dQsi1dcFdcdZT.',Z);

ad = 1;








