clear
clc

syms teta2 teta3 w2 w3
syms I1 I2b dm l2 c2 c3 K1 K2 
syms gamaTD betaTD R3
syms mi3c(w3) mi3dc(w3) Normalb0(w3)
Z = [teta2; teta3; w2; w3];
Tatr3c0 = mi3c(w3)*Normalb0(w3)*R3;
Tatr3dc0 = mi3dc(w3)*Normalb0(w3)*R3;

M0 = diag([I2b + dm*l2^2, I1]);
C0 = diag([c2, c3]);
K0 = [K1 + K2, -K2; -K2, K2];
Tc0 = [0; Tatr3c0];
Tdc0 = [0; Tatr3dc0];
B0 = [K1; 0];

Ab = [zeros(2), eye(2); -M0\K0, -M0\C0];
Fc = [zeros(2,1); M0\Tc0];
Fdc = [zeros(2,1); M0\Tdc0];
B = [zeros(2,1); M0\B0];

%%
%Forma normal
Id4 = eye(4);
C1 = Id4(4,:);
y = C1*Z;

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

%yp = Lf_h + Lg_h*teta1:
Lf_h = beta1*Z - beta0*Fc;
Lg_h = beta0B;

%ypp = L2f_h + Lg_Lf_h*teta1:
Lf_Lf_h = beta2*Z - beta1*Fc - Qsi0c.'*(Ab*Z - Fc);
Lg_Lf_h = beta1B - Qsi0c.'*B;

%yppp = L3f_h + Lg_L2f_h*teta1:
Lf_L2f_h = beta3*Z - beta2*Fc - Qsi1c.'*(Ab*Z - Fc);
Lg_L2f_h = beta2B - Qsi1c.'*B;

%Forma normal da dinâmica externa: yppp = a(Z) + b*teta1
a = Lf_L2f_h;
b = Lg_L2f_h;

%%
%Dinâmica interna:
syms q1 q2 q3 q4
syms q epsilon
A_mi = [q1 q-q1 q3 q4];
sigma = [beta0; beta1; beta2; A_mi ];
det_sigma = det(sigma);
inv_sigma = simplify(sigma\Id4);
Ad = Id4(:,4);
Bd = Id4(:,1:3);
Pi_mi = simplify(A_mi*Ab*inv_sigma);
Gama_mi = Pi_mi*Ad;
Sigma_mi = Pi_mi*Bd;
dFiDif_dZT = [zeros(1,4); 
    Qsi0c.'; 
    Qsi1c.'; 
    zeros(1,4)];
sigma_atr = sigma - dFiDif_dZT;
det_sigma_atr = det(sigma_atr);
inv_sigma_atr = simplify(sigma_atr\Id4);
Pi_mi_atr = simplify(A_mi*Ab*inv_sigma_atr);
Gama_mi_atr = simplify(Pi_mi_atr*Ad);
Sigma_mi_atr = simplify(Pi_mi_atr*Bd);
AmiB = A_mi*B;

%Expressão mais significativas:
AmiB = simplify(subs(AmiB,q3,((I2b + dm*l2^2)/K1)*epsilon));
GamaMiAtr = simplify(subs(Gama_mi_atr,q3,((I2b + dm*l2^2)/K1)*epsilon));
sigma = simplify(subs(sigma,[q1,q3,q4],[0,((I2b + dm*l2^2)/K1)*epsilon,0]));
inv_sigma = simplify(sigma\Id4);

%%
%Parâmetros e funções para a lei de controle
ad = 1;





