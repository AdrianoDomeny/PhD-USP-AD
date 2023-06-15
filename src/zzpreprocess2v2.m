function [ argumentos, repositorio ] = zzpreprocess2v2( argumentos, repositorio )
%UNTITLED32 Summary of this function goes here
%   Detailed explanation goes here
simulacao = argumentos.simulacao;
I1 = argumentos.Ip1;
%I1_hat = argumentos.Ip1_hat;
M1 = argumentos.M1rt;
%M1 = argumentos.M1;
%M1_hat = argumentos.M1_hat;
I2b = argumentos.Ip2b;
%I2b_hat = argumentos.Ip2b_hat;
dm = argumentos.dm;
%dm_hat = argumentos.dm_hat;
l2 = argumentos.l2;
%l2_hat = argumentos.l2_hat;
M2b = argumentos.M2b;
%M2b_hat = argumentos.M2b_hat;
%Dr2 = argumentos.Dr2;
G = argumentos.G;
%G_hat = argumentos.G_hat;
L = argumentos.L;
%L_hat = argumentos.L_hat;
L1 = argumentos.L1;
%L1_hat = argumentos.L1_hat;
J = argumentos.Jp;
%J_hat = argumentos.Jp_hat;
%E = argumentos.E;
%E_hat = argumentos.E_hat;
%A = argumentos.Area;
%A_hat = argumentos.Area_hat;
betaS = argumentos.betaS;
%betaS_hat = argumentos.betaS_hat;
k3 = argumentos.k3;
%k3_hat = argumentos.k3_hat;
R3 = argumentos.R3;
%R3_hat = argumentos.R3_hat;
g = argumentos.g;
%freqNatKbwEF = argumentos.freqNatKbwEF; %rad/s

%Função sinal suavizada:
nS = argumentos.nS;
signS = @(nu)(tanh(nS*pi*nu));
repositorio.signS = signS;
dsignSdnu = @(nu)(-pi*nS*(tanh(pi*nS*nu)^2 - 1));
repositorio.dsignSdnu = dsignSdnu;
%inc_dsignS = argumentos.inc_dsignS;
nf = argumentos.nf;
degrauS = @(nu)((1/2)*(tanh(nf*pi*nu) + ones(size(nu))));
repositorio.degrauS = degrauS;

%%
%Parâmetros para o metamodelo (utilizado na definição da Lei de Controle):
%K1_hat = argumentos.K1_hat;
%K2_hat = argumentos.K2_hat;
%c2_hat = argumentos.c2_hat;
%c3_hat = argumentos.c3_hat;
%{
%Parâmetros da dinâmica longitudinal utilizados na Lei de Controle:
k1_hat = E_hat*A_hat/L1_hat;
k2_hat = E_hat*A_hat/(L_hat - L1_hat);
%}

%%
%Parâmetros para modelos representativos da planta (parâmetros
%concentrados):
%(1) Dinâmica torcional:
K10 = G*J/L1;
argumentos.K10 = K10;
K20 = G*J/(L - L1);
argumentos.K20 = K20;
psic2 = argumentos.psic2;
psic3 = argumentos.psic3;
c2 = 2*psic2*sqrt((K10 + K20)/(I2b + dm*l2^2));
c3 = 2*psic3*sqrt(K20/I1);
%(2) Dinâmica longitudinal:
k1 = argumentos.k1;
k2 = argumentos.k2;
b2 = argumentos.b2;
b3 = argumentos.b3;
%
%Ponto de equilíbrio (modelos da planta, por parâmetros concentrados):
%Obs.: considerando u1_Eq = 0
uTilEqMods = argumentos.uTilEqMods;
u2_Eq = uTilEqMods(1,1);
u3_Eq = uTilEqMods(2,1);
%}
%Distância inicial entre o rotor 1 e o ponto de contato (u1 = u2 = u3 = 0):
uC = u3_Eq;
argumentos.uCmodS = uC;

%%
%Avanço máximo durante o contato:
%(a) Todo o eixo principal em regime de tração:
u31Ref = argumentos.u31Ref;
keq = (k1*k2)/(k1 + k2);
u11Ref = (1 + k3/keq)*u31Ref;
argumentos.u11Ref = u11Ref;
%(b) Apenas a parte do eixo de comprimento L1 em regime de traÃ§Ã£o:
u32Ref = argumentos.u32Ref;
u12Ref = (1 + k3/keq)*u32Ref;
argumentos.u12Ref = u12Ref;

%%
%{
%METAMODELO (para definição da lei de controle):

MtHat = diag([I2b_hat + dm_hat*l2_hat^2, I1_hat]);
M0hat = MtHat;
argumentos.M0hat = M0hat;
KtHat = [K1_hat + K2_hat, -K2_hat; -K2_hat, K2_hat];
K0hat = KtHat;
argumentos.K0hat = K0hat;
CtHat = diag([c2_hat, c3_hat]);
C0hat = CtHat;
argumentos.C0hat = C0hat;
BtHat = [K1_hat; 0];
B0hat = BtHat;
argumentos.B0hat = B0hat;

%Força normal:
Normalb0 = repositorio.Normalb0;
tetaP = repositorio.tetaP;
dNormalb0dw3 = repositorio.dNormalb0dw3;
dtetaPdw3 = repositorio.dtetaPdw3;

%Formas vetorizadas:
Normalb0vt = repositorio.Normalb0vt;
%}
%{
%Coeficiente de atrito:
%(a) Modelo contínuo: 
mic0_hat = argumentos.mic0_hat;
w0_hat = argumentos.w0_hat;
nmi = argumentos.nmi;
lambdami_hat = argumentos.lambdami_hat;
fiAtrHat = @(nu)(tanh(pi*nu) +...
    lambdami_hat*((2*nmi)/(2*nmi - 1))*nu*(1 + (1/(2*nmi - 1))*nu^(2*nmi))^-1);
repositorio.fiAtrHat = fiAtrHat;
mi3cHat = @(teta3p)(mic0_hat*fiAtrHat(teta3p/w0_hat));
repositorio.mi3cHat = mi3cHat;

%(b) Modelo descontÃ­nuo:
mi30dc_hat = argumentos.mi0_hat;
mi31dc_hat = argumentos.mi1_hat;
mi32dc_hat = argumentos.mi2_hat;
beta1mi3dc_hat = argumentos.beta1mi_hat;
beta2mi3dc_hat = argumentos.beta2mi_hat;
opcao_sign = argumentos.opcao_sign;
switch opcao_sign
    case 1
        mi3dcHat = @(teta3p)(sign(teta3p)*...
            (mi30dc_hat +...
            mi31dc_hat*(1 - 2*(1 + exp(beta1mi3dc_hat*abs(teta3p)))^-1) -...
            mi32dc_hat*(1 - 2*(1 + exp(beta2mi3dc_hat*abs(teta3p)))^-1)));
    case 2
        mi3dcHat = @(teta3p)(signS(teta3p)*...
            (mi30dc_hat +...
            mi31dc_hat*(1 - 2*(1 + exp(beta1mi3dc_hat*abs(teta3p)))^-1) -...
            mi32dc_hat*(1 - 2*(1 + exp(beta2mi3dc_hat*abs(teta3p)))^-1)));
    otherwise
        disp('other value')
end
repositorio.mi3dcHat = mi3dcHat;

%Torque de atrito:
Id2 = eye(2);
fteta3 = @(upTil)(Id2(2,:)*upTil);
repositorio.fteta3 = fteta3;
%(a) Modelo contÃ­nuo:
Tatr3c0Hat = @(teta3p)(mi3cHat(teta3p)*Normalb0(teta3p)*R3_hat);
Tc0Hat = @(uTilbp)([0; Tatr3c0Hat(fteta3(uTilbp))]);
repositorio.Tc0Hat = Tc0Hat;
%(b) Modelo descontÃ­nuo:
Tatr3dc0Hat = @(teta3p)(mi3dcHat(teta3p)*Normalb0(teta3p)*R3_hat);
Tdc0Hat = @(uTilbp)([0; Tatr3dc0Hat(fteta3(uTilbp))]);
repositorio.Tdc0Hat = Tdc0Hat;
%}
%{
%Derivadas (para integraÃ§Ã£o):
dfiAtrHatdnu = @(nu)(-pi*(tanh(pi*nu)^2 - 1)  +...
    (2*lambdami_hat*nmi)/(2*nmi + nu^(2*nmi) - 1) -...
    (4*lambdami_hat*nmi^2*nu^(2*nmi))/(2*nmi + nu^(2*nmi) - 1)^2);
repositorio.dfiAtrdnu = dfiAtrHatdnu;
dmi3cHatdw3 = @(teta3p)(mic0_hat*(1/w0_hat)*dfiAtrHatdnu(teta3p/w0_hat));
repositorio.dmi3cHatdw3 = dmi3cHatdw3;
switch opcao_sign
    case 1
        dmi3dcHatdw3 = @(teta3p)(...
            (2*beta1mi3dc_hat*mi31dc_hat*exp(abs(teta3p)*beta1mi3dc_hat))/(exp(abs(teta3p)*beta1mi3dc_hat) + 1)^2 +...
            (2*beta2mi3dc_hat*mi32dc_hat*exp(abs(teta3p)*beta2mi3dc_hat))/(exp(abs(teta3p)*beta2mi3dc_hat) + 1)^2);
    case 2
        dmi3dcHatdw3 = @(teta3p)(...
            (2*beta1mi3dc_hat*mi31dc_hat*exp(abs(teta3p)*beta1mi3dc_hat))/(exp(abs(teta3p)*beta1mi3dc_hat) + 1)^2 +...
            (2*beta2mi3dc_hat*mi32dc_hat*exp(abs(teta3p)*beta2mi3dc_hat))/(exp(abs(teta3p)*beta2mi3dc_hat) + 1)^2 +...
            inc_dsignS*dsignSdnu(teta3p)*(mi30dc_hat +...
            mi31dc_hat*(1 - 2*(1 + exp(beta1mi3dc_hat*abs(teta3p)))^-1) +...
            mi32dc_hat*(1 - 2*(1 + exp(beta2mi3dc_hat*abs(teta3p)))^-1)));
    otherwise
        disp('other value')
end
repositorio.dmi3dcHatdw3 = dmi3dcHatdw3;

dmi3cdmic0Hat = @(teta3p)(fiAtrHat(teta3p/w0_hat));
repositorio.dmi3cdmic0Hat = dmi3cdmic0Hat;
dmi3cdw0Hat = @(teta3p)(mic0_hat*(-teta3p/w0_hat^2)*dfiAtrHatdnu(teta3p/w0_hat));
repositorio.dmi3cdw0Hat = dmi3cdw0Hat;
dfiAtrdlambdami = @(nu)((2*nmi*nu)/(2*nmi + nu^(2*nmi) - 1));
dmi3cdlambdamiHat = @(teta3p)(mic0_hat*dfiAtrdlambdami(teta3p/w0_hat));
repositorio.dmi3cdlambdamiHat = dmi3cdlambdamiHat;
d2fiAtrHatdnu2 = @(nu)(2*pi^2*tanh(pi*nu)*(tanh(pi*nu)^2 - 1) -...
    (8*lambdami_hat*nmi^2*nu^(2*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^2 +...
    (16*lambdami_hat*nmi^3*nu^(4*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^3 -...
    (4*lambdami_hat*nmi^2*nu^(2*nmi - 1)*(2*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^2);
d2mi3cHatdw32 = @(teta3p)(mic0_hat*(1/w0_hat^2)*d2fiAtrHatdnu2(teta3p/w0_hat));
repositorio.d2mi3cHatdw32 = d2mi3cHatdw32;
d2mi3cdmic0dw3Hat = @(teta3p)((1/w0_hat)*dfiAtrHatdnu(teta3p/w0_hat));
repositorio.d2mi3cdmic0dw3Hat = d2mi3cdmic0dw3Hat;
d3mi3cdmic0dw32Hat = @(teta3p)((1/w0_hat^2)*d2fiAtrHatdnu2(teta3p/w0_hat));
repositorio.d3mi3cdmic0dw32Hat = d3mi3cdmic0dw32Hat;
d2mi3cdw0dw3Hat = @(teta3p)(mic0_hat*...
    ((-1/w0_hat^2)*dfiAtrHatdnu(teta3p/w0_hat) + (1/w0_hat)*d2fiAtrHatdnu2(teta3p/w0_hat)*(-teta3p/w0_hat^2)));
repositorio.d2mi3cdw0dw3Hat = d2mi3cdw0dw3Hat;
d3fiAtrHatdnu3 = @(nu)(-2*pi^3*(3*tanh(pi*nu)^4 - 4*tanh(pi*nu)^2 + 1) +...
    (16*lambdami_hat*nmi^3*nu^(4*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^3 +...
    (64*lambdami_hat*nmi^4*nu^(4*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^3 -...
    (96*lambdami_hat*nmi^4*nu^(6*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^4 -...
    (4*lambdami_hat*nmi^2*nu^(2*nmi - 2)*(2*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^2 -...
    (8*lambdami_hat*nmi^3*nu^(2*nmi - 2)*(2*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^2 +...
    (16*lambdami_hat*nmi^3*nu^(4*nmi - 2)*(2*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^3);
d3mi3cHatdw33 = @(teta3p)(mic0_hat*(1/w0_hat^3)*d3fiAtrHatdnu3(teta3p/w0_hat));
repositorio.d3mi3cHatdw33 = d3mi3cHatdw33;
d3mi3cdw0dw32Hat = @(teta3p)(mic0_hat*...
    ((-2/w0_hat^3)*d2fiAtrHatdnu2(teta3p/w0_hat) + (1/w0_hat^2)*d3fiAtrHatdnu3(teta3p/w0_hat)*(-teta3p/w0_hat^2)));
repositorio.d3mi3cdw0dw32Hat = d3mi3cdw0dw32Hat;
dfiAtrdlambdami = @(nu)((2*nmi*nu)/(2*nmi + nu^(2*nmi) - 1));
dmi3cdlambdami = @(teta3p)(mic0_hat*dfiAtrdlambdami(teta3p/w0_hat));
repositorio.dmi3cdlambdami = dmi3cdlambdami;
d2fiAtrdlambdamiHatdnu = @(nu)(-(2*nmi*(2*nmi - 1)*(nu^(2*nmi) - 1))/...
    (2*nmi + nu^(2*nmi) - 1)^2);
d2mi3cdlambdamidw3Hat = @(teta3p)(mic0_hat*(1/w0_hat)*d2fiAtrdlambdamiHatdnu(teta3p/w0_hat));
repositorio.d2mi3cdlambdamidw3Hat = d2mi3cdlambdamidw3Hat;
d3fiAtrdlambdamidnu2 = @(nu)(-(4*nmi^2*nu^(2*nmi - 1)*(2*nmi - 1)*(2*nmi - nu^(2*nmi) + 1))/...
    ((2*nmi + nu^(2*nmi) - 1)^3));
d3mi3cdlambdamidw32Hat = @(teta3p)(mic0_hat*(1/w0_hat^2)*d3fiAtrdlambdamidnu2(teta3p/w0_hat));
repositorio.d3mi3cdlambdamidw32Hat = d3mi3cdlambdamidw32Hat;
d4mi3cdmic0dw33Hat = @(teta3p)((1/w0_hat^3)*d3fiAtrHatdnu3(teta3p/w0_hat));
repositorio.d4mi3cdmic0dw33Hat = d4mi3cdmic0dw33Hat;
d4fiAtrHatdnu4 = @(nu)(-2*pi^3*(-4*pi*tanh(pi*nu)*(3*tanh(pi*nu)^4 - 5*tanh(pi*nu)^2 + 2)) -...
    (96*lambdami_hat*nmi^4*nu^(6*nmi - 3))/(2*nmi + nu^(2*nmi) - 1)^4 -...
    (384*lambdami_hat*nmi^5*nu^(6*nmi - 3))/(2*nmi + nu^(2*nmi) - 1)^4 +...
    (768*lambdami_hat*nmi^5*nu^(8*nmi - 3))/(2*nmi + nu^(2*nmi) - 1)^5 +...
    (16*lambdami_hat*nmi^3*nu^(4*nmi - 3)*(2*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^3 +...
    (32*lambdami_hat*nmi^4*nu^(4*nmi - 3)*(2*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^3 +...
    (16*lambdami_hat*nmi^3*nu^(4*nmi - 3)*(4*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^3 +...
    (64*lambdami_hat*nmi^4*nu^(4*nmi - 3)*(4*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^3 -...
    (96*lambdami_hat*nmi^4*nu^(6*nmi - 3)*(2*nmi - 1))/(2*nmi + nu^(2*nmi) - 1)^4 -...
    (96*lambdami_hat*nmi^4*nu^(6*nmi - 3)*(6*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^4 -...
    (4*lambdami_hat*nmi^2*nu^(2*nmi - 3)*(2*nmi - 1)*(2*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^2 -...
    (8*lambdami_hat*nmi^3*nu^(2*nmi - 3)*(2*nmi - 1)*(2*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^2 +...
    (16*lambdami_hat*nmi^3*nu^(4*nmi - 3)*(2*nmi - 1)*(4*nmi - 2))/(2*nmi + nu^(2*nmi) - 1)^3);
d4mi3cHatdw34 = @(teta3p)(mic0_hat*(1/w0_hat^4)*d4fiAtrHatdnu4(teta3p/w0_hat));
repositorio.d4mi3cHatdw34 = d4mi3cHatdw34;
d4mi3cdw0dw33Hat = @(teta3p)(mic0_hat*...
    ((-3/w0_hat^4)*d3fiAtrHatdnu3(teta3p/w0_hat) + (1/w0_hat^3)*d4fiAtrHatdnu4(teta3p/w0_hat)*(-teta3p/w0_hat^2)));
repositorio.d4mi3cdw0dw33Hat = d4mi3cdw0dw33Hat;
d4fiAtrdlambdamidnu3 = @(nu)((4*nmi^2*nu^(2*(nmi-1))*(2*nmi - 1)*...
    (2*nmi - 2*nmi*nu^(4*nmi) + 2*nu^(2*nmi) - nu^(4*nmi) +...
    16*nmi^2*nu^(2*nmi) + 4*nmi^2 - 8*nmi^3 - 1))/...
    ((2*nmi + nu^(2*nmi) - 1)^4));
d4mi3cdlambdamidw33Hat = @(teta3p)(mic0_hat*(1/w0_hat^3)*d4fiAtrdlambdamidnu3(teta3p/w0_hat));
repositorio.d4mi3cdlambdamidw33Hat = d4mi3cdlambdamidw33Hat;
d5mi3cdmic0dw34Hat = @(teta3p)((1/w0_hat^4)*d4fiAtrHatdnu4(teta3p/w0_hat));
repositorio.d5mi3cdmic0dw34Hat = d5mi3cdmic0dw34Hat;
d5fiAtrHatdnu5 = @(nu)(-2*pi^3*(-4*pi*(-pi*(17*tanh(pi*nu)^2 - 30*tanh(pi*nu)^4 + 15*tanh(pi*nu)^6 - 2))) +...
    (8*lambdami_hat*nmi^2*nu^(2*nmi-4)*(2*nmi - 1)*...
     (23*nmi - 58*nmi*nu^(2*nmi) + 36*nmi*nu^(4*nmi) + 10*nmi*nu^(6*nmi) -...
      11*nmi*nu^(8*nmi) + 12*nu^(2*nmi) - 18*nu^(4*nmi) + 12*nu^(6*nmi) -...
      3*nu^(8*nmi) + 108*nmi^2*nu^(2*nmi) - 320*nmi^3*nu^(2*nmi) -...
      60*nmi^2*nu^(4*nmi) + 1088*nmi^4*nu^(2*nmi) + 120*nmi^3*nu^(4*nmi) -...
      1632*nmi^5*nu^(2*nmi) + 20*nmi^2*nu^(6*nmi) + 528*nmi^4*nu^(4*nmi) +...
      832*nmi^6*nu^(2*nmi) + 200*nmi^3*nu^(6*nmi) - 1056*nmi^5*nu^(4*nmi) -...
      12*nmi^2*nu^(8*nmi) + 208*nmi^4*nu^(6*nmi) - 4*nmi^3*nu^(8*nmi) -...
      56*nmi^2 + 4*nmi^3 + 208*nmi^4 - 368*nmi^5 + 256*nmi^6 - 64*nmi^7 - 3))/...
      ((2*nmi + nu^(2*nmi) - 1)^6)); 
%Obs.: ExpressÃ£o obtida a partir d4fiAtrHatdnu4, mas nÃ£o verificada simb.
d5mi3cdw0dw34Hat = @(teta3p)(mic0_hat*...
    ((-3/w0_hat^5)*d4fiAtrHatdnu4(teta3p/w0_hat) +...
     (1/w0_hat^3)*(d5fiAtrHatdnu5(teta3p/w0_hat)*(1/w0_hat)*(-teta3p/w0_hat^2) +...
                   d4fiAtrHatdnu4(teta3p/w0_hat)*(-1/w0_hat^2))));
repositorio.d5mi3cdw0dw34Hat = d5mi3cdw0dw34Hat;
d5fiAtrdlambdamidnu4 = @(nu)(-(8*nmi^2*nu^(2*nmi-3)*(2*nmi - 1)*...
    (7*nmi*nu^(2*nmi) - 5*nmi + nmi*nu^(4*nmi) - 3*nmi*nu^(6*nmi) -...
     3*nu^(2*nmi) + 3*nu^(4*nmi) - nu^(6*nmi) - 14*nmi^2*nu^(2*nmi) +...
     68*nmi^3*nu^(2*nmi) + 12*nmi^2*nu^(4*nmi) - 88*nmi^4*nu^(2*nmi) +...
     44*nmi^3*nu^(4*nmi) - 2*nmi^2*nu^(6*nmi) + 4*nmi^2 + 16*nmi^3 -...
     32*nmi^4 + 16*nmi^5 + 1))/...
    ((2*nmi + nu^(2*nmi) - 1)^5));
d5mi3cdlambdamidw34Hat = @(teta3p)(mic0_hat*(1/w0_hat^4)*d5fiAtrdlambdamidnu4(teta3p/w0_hat));
repositorio.d5mi3cdlambdamidw34Hat = d5mi3cdlambdamidw34Hat;
%}

%%
%{
%EQUAÇÃO DE ESTADOS, COM PARÂMETROS ESTIMADOS: 
%Zp = AbHat*Z - FcHat(Z) + Bhat*Uc
%Obs.: Z = [uTilb; uTilbp]
AbHat = [zeros(2), eye(2); -M0hat\K0hat, -M0hat\C0hat];
argumentos.AbHat = AbHat;
Id4 = eye(4);
Au = Id4(:,1:2);
argumentos.Au = Au;
Av = Id4(:,3:4);
argumentos.Av = Av;
%Obs.: Z = [uTilb; uTilbp]
fuTil = @(Z)(Au.'*Z);
repositorio.fuTil = fuTil;
fupTil = @(Z)(Av.'*Z);
repositorio.fuTil = fuTil;
FcHat = @(Z)([zeros(2,1); M0hat\Tc0Hat(fupTil(Z))]);
repositorio.FcHat = FcHat;
FdcHat = @(Z)([zeros(2,1); M0hat\Tdc0Hat(fupTil(Z))]);
repositorio.FdcHat = FdcHat;
Bhat = [zeros(2,1); M0hat\B0hat];
argumentos.Bhat = Bhat;
%FunÃ§Ãµes Ãºteis:
fcZp = @(Z,Uc)(AbHat*Z - FcHat(Z) + Bhat*Uc);
repositorio.fcZp = fcZp;
fdcZp = @(Z,Uc)(AbHat*Z - FdcHat(Z) + Bhat*Uc);
repositorio.fdcZp = fdcZp;
%}

%%
%{
%FORMA NORMAL:
%(1) Termos auxiliares:
beta0 = argumentos.beta0;
beta1Hat = argumentos.beta1Hat;
beta2Hat = argumentos.beta2Hat;
beta3Hat = argumentos.beta3Hat;

Qsi0cHat = repositorio.Qsi0cHat;
Qsi1cHat = repositorio.Qsi1cHat;

Qsi0dcHat = @(Z)(fQsi0dcHat( Z, argumentos, repositorio ));
repositorio.Qsi0dcHat = Qsi0dcHat;
Qsi1dcHat = @(Z)(fQsi1dcHat( Z, argumentos, repositorio ));
repositorio.Qsi1dcHat = Qsi1dcHat;
%Qsi2dc_hat = @(Z)(fQsi2dc_hat( Z, argumentos, repositorio ));

%(2) DinÃ¢mica externa (forma normal): yppp = a[d]c(Z) + b*Uc
acHat = @(Z)(beta3Hat*Z - beta2Hat*FcHat(Z) - Qsi1cHat(Z)*(AbHat*Z - FcHat(Z)));
repositorio.acHat = acHat;
adcHat = @(Z)(beta3Hat*Z - beta2Hat*FdcHat(Z) - Qsi1dcHat(Z).'*(AbHat*Z - FdcHat(Z)));
repositorio.adcHat = adcHat;
bHat = beta2Hat*Bhat;
argumentos.bHat = bHat;
%}

%{
%(3) DinÃ¢mica interna: mi = Ami*Z
epsilon = 10^-5;
q = epsilon;
%q = 1;
argumentos.qAmi = q;
%q1 = 0;
q1 = 1;
q3 = ((I2b_hat + dm_hat*l2_hat^2)/K1_hat)*epsilon;
q4 = 0;
Ami = [q1, q-q1, q3, q4];
argumentos.Ami = Ami;

%(4) Difeomorfismo:
%Obs.: X = [y; yp; ypp; mi] e X = PsiDif(Z) = sigma*Z - FiDif(Z)
sigmaHat = [beta0; beta1Hat; beta2Hat; Ami];
argumentos.sigmaHat = sigmaHat;
FiDifcHat = @(Z)([0;
    beta0*FcHat(Z); 
    beta1Hat*FcHat(Z) + Qsi0cHat(Z)*(AbHat*Z - FcHat(Z));
    0]);
FiDifdcHat = @(Z)([0;
    beta0*FdcHat(Z); 
    beta1Hat*FdcHat(Z) + Qsi0dcHat(Z).'*(AbHat*Z - FdcHat(Z));
    0]);
PsiDifcHat = @(Z)(sigmaHat*Z - FiDifcHat(Z));
repositorio.PsiDifcHat = PsiDifcHat;
PsiDifdcHat = @(Z)(sigmaHat*Z - FiDifdcHat(Z));
repositorio.PsiDifdcHat = PsiDifdcHat;
dFiDifcHatdZT = @(Z)(...
    [zeros(1,4); 
    Qsi0cHat(Z); 
    Qsi1cHat(Z);
    zeros(1,4)]);
dPsiDifcHatdZT = @(Z)(sigmaHat - dFiDifcHatdZT(Z));
repositorio.dPsiDifcHatdZT = dPsiDifcHatdZT;
invPsiDifcHat = @(X,Z0)(fInvPsiDifcHat( X, Z0, argumentos, repositorio ));
repositorio.invPsiDifcHat = invPsiDifcHat;
dInvPsiDifcdXTHat = @(X,Z0)(dPsiDifcHatdZT(invPsiDifcHat(X,Z0))\Id4);
repositorio.dInvPsiDifcdXTHat = dInvPsiDifcdXTHat;
dFiDifdcHatdZT = @(Z)(...
    [zeros(1,4); 
    Qsi0dcHat(Z).'; 
    Qsi1dcHat(Z).';
    zeros(1,4)]);
dPsiDifdcHatdZT = @(Z)(sigmaHat - dFiDifdcHatdZT(Z));
repositorio.dPsiDifdcHatdZT = dPsiDifdcHatdZT;
invPsiDifdcHat = @(X,Z0)(fInvPsiDifdcHat( X, Z0, argumentos, repositorio ));
repositorio.invPsiDifdcHat = invPsiDifdcHat;
dInvPsiDifdcdXTHat = @(X,Z0)(dPsiDifdcHatdZT(invPsiDifdcHat(X,Z0))\Id4);
repositorio.dInvPsiDifdcdXTHat = dInvPsiDifdcdXTHat;
%}

%%
%LEI DE CONTROLE
%{
%(1) TrajetÃ³ria desejada: yd = w3d
t0_9wRef = argumentos.t0_9wRef;
wRef = argumentos.wRef;
alfa = argumentos.alfa_trajetoria_desejada;
fi = @(nu)(tanh(alfa*nu) - (alfa*nu)/((alfa^2*nu^2)/3 + 1));
dfidnu = @(nu)((6*alfa^3*nu^2)/(alfa^2*nu^2 + 3)^2 -...
    (3*alfa)/(alfa^2*nu^2 + 3) -...
    alfa*(tanh(alfa*nu)^2 - 1));
d2fidnu2 = @(nu)(2*alfa^2*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1) +...
    (18*alfa^3*nu)/(alfa^2*nu^2 + 3)^2 -...
    (24*alfa^5*nu^3)/(alfa^2*nu^2 + 3)^3);
d3fidnu3 = @(nu)((18*alfa^3)/(alfa^2*nu^2 + 3)^2 -...
    2*alfa^3*(tanh(alfa*nu)^2 - 1)^2 -...
    4*alfa^3*tanh(alfa*nu)^2*(tanh(alfa*nu)^2 - 1) -...
    (144*alfa^5*nu^2)/(alfa^2*nu^2 + 3)^3 +...
    (144*alfa^7*nu^4)/(alfa^2*nu^2 + 3)^4);
d4fidnu4 = @(nu)(16*alfa^4*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1)^2 -...
    (360*alfa^5*nu)/(alfa^2*nu^2 + 3)^3 +...
    8*alfa^4*tanh(alfa*nu)^3*(tanh(alfa*nu)^2 - 1) +...
    (1440*alfa^7*nu^3)/(alfa^2*nu^2 + 3)^4 -...
    (1152*alfa^9*nu^5)/(alfa^2*nu^2 + 3)^5);
yd = @(t)(wRef*fi(t/t0_9wRef));
repositorio.yd = yd;
ypd = @(t)((wRef/t0_9wRef)*dfidnu(t/t0_9wRef));
repositorio.ypd = ypd;
yppd = @(t)((wRef/t0_9wRef^2)*d2fidnu2(t/t0_9wRef));
repositorio.yppd = yppd;
ypppd = @(t)((wRef/t0_9wRef^3)*d3fidnu3(t/t0_9wRef));
repositorio.ypppd = ypppd;
yppppd = @(t)((wRef/t0_9wRef^4)*d4fidnu4(t/t0_9wRef));
repositorio.yppppd = yppppd;
%VetorizaÃ§Ã£o:
fiVt = @(nu)(tanh(alfa*nu) - (alfa*nu)./((alfa^2*nu.^2)/3 + ones(size(nu))));
dfidnuVt = @(nu)((6*alfa^3*nu.^2)./(alfa^2*nu.^2 + 3*ones(size(nu))).^2 -...
    (3*alfa*ones(size(nu)))./(alfa^2*nu.^2 + 3*ones(size(nu))) -...
    alfa*(tanh(alfa*nu).^2 - ones(size(nu))));
d2fidnu2Vt = @(nu)(2*alfa^2*tanh(alfa*nu).*(tanh(alfa*nu).^2 - ones(size(nu))) +...
    (18*alfa^3*nu)./(alfa^2*nu.^2 + 3*ones(size(nu))).^2 -...
    (24*alfa^5*nu.^3)./(alfa^2*nu.^2 + 3*ones(size(nu))).^3);
d3fidnu3Vt = @(nu)((18*alfa^3*ones(size(nu)))./(alfa^2*nu.^2 + 3*ones(size(nu))).^2 -...
    2*alfa^3*(tanh(alfa*nu).^2 - ones(size(nu))).^2 -...
    4*alfa^3*tanh(alfa*nu).^2.*(tanh(alfa*nu).^2 - ones(size(nu))) -...
    (144*alfa^5*nu.^2)./(alfa^2*nu.^2 + 3*ones(size(nu))).^3 +...
    (144*alfa^7*nu.^4)./(alfa^2*nu.^2 + 3*ones(size(nu))).^4);
ydVt = @(T)(wRef*fiVt(T/t0_9wRef));
repositorio.ydVt = ydVt;
ypdVt = @(T)((wRef/t0_9wRef)*dfidnuVt(T/t0_9wRef));
repositorio.ypdVt = ypdVt;
yppdVt = @(T)((wRef/t0_9wRef^2)*d2fidnu2Vt(T/t0_9wRef));
repositorio.yppdVt = yppdVt;
ypppdVt = @(T)((wRef/t0_9wRef^3)*d3fidnu3Vt(T/t0_9wRef));
repositorio.ypppdVt = ypppdVt;
%}

%(2) Lei de controle:
%dqsiTil = argumentos.dqsiTil;
%daChatdqsi = @(Z)(fdaChatdqsi( Z, argumentos, repositorio ));
%daC = @(Z)(dqsiTil.'*daChatdqsi(Z));
%daDcHatdqsi = @(Z)(fdaDcHatdqsi( Z, argumentos, repositorio ));
%daDc = @(Z)(dqsiTil.'*daDcHatdqsi(Z));
%dbHatdqsi = fdbHatdqsi( argumentos );
%dbHat = dqsiTil.'*dbHatdqsi;

%sigmaDmax = bHat/abs(dbHat);
%argumentos.sigmaDmax = sigmaDmax;

%{
%sigmaF = argumentos.sigmaF;
sigmaF = 1;
FcLC = @(Z)(sigmaF*abs(daC(Z)));
repositorio.FcLC = FcLC;
FdcLC = @(Z)(sigmaF*abs(daDc(Z)));
repositorio.FdcLC = FdcLC;
%sigmaD = argumentos.sigmaD;
sigmaD = 1;
Dlc = sigmaD*(1/bHat)*abs(dbHat);
argumentos.Dlc = Dlc;

%Obs.: X = [y; yp; ypp; mi]
y = @(X)(Id4(1,:)*X);
yp = @(X)(Id4(2,:)*X);
ypp = @(X)(Id4(3,:)*X);
%}

%{
%Frequência de corte:
N = argumentos.N;
rEstr = argumentos.rEstr;
sig2 = argumentos.sig2;
wbNatTil = freqNatKbwEF(3:6*(N+1)-rEstr,1);
nLmb = 0.25;
wCtr2 = sig2;
wCtr3 = wbNatTil(4,1); %Menor frequência não modelada da planta
lambda_1 = (1-nLmb)*wCtr2 + nLmb*wCtr3;
argumentos.lambda_ctr_y1 = lambda_1;
%}

%{
y2r = @(t,X)(yppd(t) - 2*lambda_1*(yp(X) - ypd(t)) -...
    lambda_1^2*(y(X) - yd(t)));
s = @(t,X)(ypp(X) - y2r(t,X));
repositorio.s = s;
y3r = @(t,X)(ypppd(t) - 2*lambda_1*(ypp(X) - yppd(t)) -...
    lambda_1^2*(yp(X) - ypd(t)));
repositorio.y3r = y3r;
eta = 1; %s
argumentos.eta = eta;
UcTilc = @(t,Z,X)(y3r(t,X) - acHat(Z));
repositorio.UcTilc = UcTilc; 
Kc = @(t,Z,X)((1 - Dlc)^-1*(FcLC(Z) + Dlc*abs(UcTilc(t,Z,X)) + eta));
repositorio.Kc = Kc;
UcTildc = @(t,Z,X)(y3r(t,X) - adcHat(Z));
Kdc = @(t,Z,X)((1 - Dlc)^-1*(FdcLC(Z) + Dlc*abs(UcTildc(t,Z,X)) + eta));
Kfi = lambda_1/(1 - Dlc);
argumentos.Kfi = Kfi;
%VetorizaÃ§Ã£o:
yVt = @(XX)(XX*Id4(:,1));
ydVt = @(T)(wRef*fiVt(T/t0_9wRef));
ypVt = @(XX)(XX*Id4(:,2));
ypdVt = @(T)((wRef/t0_9wRef)*dfidnuVt(T/t0_9wRef));
yppVt = @(XX)(XX*Id4(:,3));
y2rVt = @(T,XX)(yppdVt(T) - 2*lambda_1*(ypVt(XX) - ypdVt(T)) -...
    lambda_1^2*(yVt(XX) - ydVt(T)));
sVt = @(T,XX)(yppVt(XX) - y2rVt(T,XX));
repositorio.sVt = sVt;
%}

%{
%FunÃ§Ãµes adaptadas para a trajetoria desejada:
yTild = @(t)([yd(t); ypd(t); yppd(t)]);
Ad = Id4(:,4);
argumentos.Ad = Ad;
Bd = Id4(:,1:3);
argumentos.Bd = Bd;

fXd = @(t,mid)(Ad*mid + Bd*yTild(t));
repositorio.fXd = fXd;
fcZdHat = @(t,mid,Z0)(invPsiDifcHat(fXd(t,mid),Z0));
repositorio.fcZdHat = fcZdHat;
fdcZdHat = @(t,mid,Z0)(invPsiDifdcHat(fXd(t,mid),Z0));
repositorio.fdcZdHat = fdcZdHat;
UcTilcd = @(t,mid,Z0)(UcTilc(t,fcZdHat(t,mid,Z0),fXd(t,mid)));
FcLCd = @(t,mid,Z0)(FcLC(fcZdHat(t,mid,Z0)));
Kcd = @(t,mid,Zd0)((1 - Dlc)^-1*(FcLCd(t,mid,Zd0) + Dlc*abs(UcTilcd(t,mid,Zd0)) + eta));
repositorio.Kcd = Kcd;
UcTildcd = @(t,mid,Z0)(UcTildc(t,fdcZdHat(t,mid,Z0),fXd(t,mid)));
FdcLCd = @(t,mid,Z0)(FdcLC(fdcZdHat(t,mid,Z0)));
Kdcd = @(t,mid,Z0)((1 - Dlc)^-1*(FdcLCd(t,mid,Z0) + Dlc*abs(UcTildcd(t,mid,Z0)) + eta));
repositorio.Kdcd = Kdcd;

fZ = @(uTilb,uTilbp)(Au*uTilb + Av*uTilbp);
repositorio.fZ = fZ;
fcXhat = @(uTilb,uTilbp)(PsiDifcHat(fZ(uTilb,uTilbp)));
repositorio.fcXhat = fcXhat;
fdcXhat = @(uTilb,uTilbp)(PsiDifdcHat(fZ(uTilb,uTilbp)));
repositorio.fdcXhat = fdcXhat;

%Camada limite dinÃ¢mica:
DfiBd = @(fipBd)(double(fipBd>=0)*1/(1 + Dlc) + double(fipBd<=0)*1/(1 - Dlc));
repositorio.Dfi = DfiBd;
GamaFiBd = @(fipBd)(DfiBd(fipBd)*fipBd);
repositorio.GamaFiBd = GamaFiBd;

%ExpressÃµes de ganhos e lei de controle:
Kcb = @(t,Z,X,fiBd,mid,Zd0)(Kc(t,Z,X) - Kcd(t,mid,Zd0) + Kfi*fiBd);
repositorio.Kcb = Kcb;
Kdcb = @(t,Z,X,fiBd,mid,Zd0)(Kdc(t,Z,X) - Kdcd(t,mid,Zd0) + Kfi*fiBd);
sat = @(nu)(tanh(pi*nu));
repositorio.sat = sat;
Ucc = @(t,Z,X,fiBd,mid,Zd0)(bHat^-1*(UcTilc(t,Z,X) - Kcb(t,Z,X,fiBd,mid,Zd0)*sat(s(t,X)/fiBd)));
repositorio.Ucc = Ucc;
Ucdc = @(t,Z,X,fiBd,mid,Zd0)(bHat^-1*(UcTildc(t,Z,X) - Kdcb(t,Z,X,fiBdTil,mid,Zd0)*sat(s(t,X)/fiBd)));
repositorio.Ucdc = Ucdc;

%ExpressÃµes para integraÃ§Ã£o:
Uccu = @(t,uTilc,uTilcp,fiBd,mid,Zd0)(Ucc(t,fZ(uTilc,uTilcp),fcXhat(uTilc,uTilcp),fiBd,mid,Zd0));
repositorio.Uccu = Uccu;
Ucdcu = @(t,uTilc,uTilcp,fiBd,mid,Zd0)(Ucdc(t,fZ(uTilc,uTilcp),fdcXhat(uTilc,uTilcp),fiBd,mid,Zd0));
repositorio.Ucdcu = Ucdcu;
%}

%{
%ExpressÃµes para a dinÃ¢mica interna:
y3rMid = @(t)(ypppd(t));
UcTilcMid = @(t,mid,Zd0)(y3rMid(t) - acHat(fcZdHat(t,mid,Zd0)));
UccMid = @(t,mid,Zd0)(bHat^-1*(UcTilcMid(t,mid,Zd0)));
repositorio.UccMid = UccMid;
UcTildcMid = @(t,mid,Zd0)(y3rMid(t) - adcHat(fdcZdHat(t,mid,Zd0)));
UcdcMid = @(t,mid,Zd0)(bHat^-1*(UcTildcMid(t,mid,Zd0)));
repositorio.UcdcMid = UcdcMid;

%Derivadas (para integraÃ§Ã£o):
dGamadfipBd = @(fipBd)(DfiBd(fipBd));
repositorio.dGamadfipBd = dGamadfipBd;

%dypdXT = Id4(2,:);
%dyppdXT = Id4(3,:);
%dy3rdXT = -2*lambda_1*dyppdXT -...
%    lambda_1^2*dypdXT;
%dUcTildXT = dy3rdXT;
dKcdXdT = @(t,mid,Z0)((1 - Dlc)^-1*Dlc*dUcTildXT*sign(UcTilcd(t,mid,Z0)));
d2aChatdZTdqsi = @(Z)(fd2aChatdZTdqsi( Z, argumentos, repositorio ));
ddaCdZT = @(Z)(dqsiTil.'*d2aChatdZTdqsi(Z));
dFcLCdZT = @(Z)(ddaCdZT(Z)*sign(daC(Z)));
dQsi1cTHatdZT = @(Z)(fdQ1cHatdZT( Z, argumentos, repositorio ));
%dQsi1AZcHatdZT = @(Z)(dQsi1cHatdZT(Z).'*(AbHat*Z) + AbHat.'*Qsi1cHat(Z));
dQsi1AZcHatdZT = @(Z)((AbHat*Z).'*dQsi1cTHatdZT(Z) + Qsi1cHat(Z)*AbHat);
dTatr3cHatdw3 = @(w3)(dmi3cHatdw3(w3)*Normalb0(w3) + mi3cHat(w3)*dNormalb0dw3(w3))*R3_hat;
dTatr3cHatduTilbpT = @(w3)([0, dTatr3cHatdw3(w3)]);
fw3 = repositorio.fw3;
dTc0HatduTilbpT = @(uTilbp)([zeros(1,2); dTatr3cHatduTilbpT(fw3(uTilbp))]);
duTilbpdZT = Av.';
dTc0HatdZT = @(uTilbp)(dTc0HatduTilbpT(uTilbp)*duTilbpdZT);
dFcHatdZT = @(Z)([zeros(2,4); M0hat\dTc0HatdZT(fupTil(Z))]);
%dQsi1FcHatdZT = @(Z)(dQsi1cTHatdZT(Z).'*FcHat(Z) + dFcHatdZT(Z).'*Qsi1cHat(Z));
dQsi1FcHatdZT = @(Z)(FcHat(Z).'*dQsi1cTHatdZT(Z) + Qsi1cHat(Z)*dFcHatdZT(Z));
daChatdZT = @(Z)(beta3Hat - beta2Hat*dFcHatdZT(Z) - dQsi1AZcHatdZT(Z) + dQsi1FcHatdZT(Z));
dUcTilcdZdT = @(t,mid,Z0)(-daChatdZT(fcZdHat(t,mid,Z0)));
dFcLCdZdT = @(t,mid,Zd0)(dFcLCdZT(fcZdHat(t,mid,Zd0)));
dKcdZdT = @(t,mid,Z0)((1 - Dlc)^-1*...
    (dFcLCdZdT(t,mid,Z0) + Dlc*dUcTilcdZdT(t,mid,Z0)*sign(UcTilcd(t,mid,Z0))));
DKcDXdT = @(t,mid,Z0)(dKcdXdT(t,mid,Z0) + dKcdZdT(t,mid,Z0)*dInvPsiDifcdXTHat(fXd(t,mid),Z0));
repositorio.DKcDXdT = DKcDXdT;

dKdcdXdT = @(t,mid,Z0)((1 - Dlc)^-1*Dlc*dUcTildXT*sign(UcTildcd(t,mid,Z0)));
d2aDcHatdZTdqsi = @(Z)(f2daDcHatdZTdqsi( Z, argumentos, repositorio ));
ddaDcdZT = @(Z)(dqsiTil.'*d2aDcHatdZTdqsi(Z));
dFdcLCdZT = @(Z)(ddaDcdZT(Z)*sign(daDc(Z)));
dQsi1dcHatdZT = @(Z)(fdQsi1dcHatdZT( Z, argumentos, repositorio ));
dQsi1AZdcHatdZT = @(Z)(dQsi1dcHatdZT(Z).'*(AbHat*Z) + Qsi1dcHat(Z).'*AbHat);
dTatr3dcHatdw3 = @(w3)(dmi3dcHatdw3(w3)*Normalb0(w3) + mi3dcHat(w3)*dNormalb0dw3(w3))*R3_hat;
dTatr3dcHatduTilbpT = @(w3)([0, dTatr3dcHatdw3(w3)]);
dTdc0HatduTilbpT = @(uTilbp)([zeros(1,2); dTatr3dcHatduTilbpT(fw3(uTilbp))]);
dTdc0HatdZT = @(uTilbp)(dTdc0HatduTilbpT(uTilbp)*duTilbpdZT);
dQsi1dcHatdZT = @(Z)(fdQsi1dcHatdZT( Z, argumentos, repositorio ));
dFdcHatdZT = @(Z)([zeros(2,4); M0hat\dTdc0HatdZT(fupTil(Z))]);
dQsiFdcHatdZT = @(Z)(dQsi1dcHatdZT(Z).'*FdcHat(Z) + dFdcHatdZT(Z).'*Qsi1dcHat(Z));
daDcHatdZT = @(Z)(beta3Hat - beta2Hat*dFdcHatdZT(Z) - dQsi1AZdcHatdZT(Z) + dQsiFdcHatdZT(Z));
dUcTildcdZdT = @(t,mid,Z0)(-daDcHatdZT(fdcZdHat(t,mid,Z0)));
dFdcLCdZdT = @(t,mid,Z0)(dFdcLCdZT(fdcZdHat(t,mid,Z0)));
dKdcdZdT = @(t,mid,Z0)((1 - Dlc)^-1*...
    (dFdcLCdZdT(t,mid,Z0) + Dlc*dUcTildcdZdT(t,mid,Z0)*sign(UcTildcd(t,mid,Z0))));
DKdcDXdT = @(t,mid,Z0)(dKdcdXdT(t,mid,Z0) + dKdcdZdT(t,mid,Z0)*dInvPsiDifdcdXTHat(fXd(t,mid),Z0));
repositorio.DKdcDXdT = DKdcDXdT;

dUcTilcdZT = @(Z)(-daChatdZT(Z));
dKcdZT = @(t,Z,X)((1 - Dlc)^-1*(dFcLCdZT(Z) + Dlc*dUcTilcdZT(Z)*sign(UcTilc(t,Z,X))));
dKcbdZT = @(t,Z,X)(dKcdZT(t,Z,X));
dUccdZT = @(t,Z,X,fiBd)(bHat^-1*(dUcTilcdZT(Z) - dKcbdZT(t,Z,X)*sat(s(t,X)/fiBd)));
dKcdXT = @(t,Z,X)((1 - Dlc)^-1*(Dlc*dUcTildXT*sign(UcTilc(t,Z,X))));
dKcbdXT = @(t,Z,X,fiBd,mid,Z0)(dKcdXT(t,Z,X));
dsatdnu = @(nu)(-pi*(tanh(pi*nu)^2 - 1));
repositorio.dsatdnu = dsatdnu;
d2satdnu2 = @(nu)(2*pi^2*tanh(pi*nu)*(tanh(pi*nu)^2 - 1));
repositorio.d2satdnu2 = d2satdnu2;
%dydXT = Id4(1,:);
%dy2rdXT = -2*lambda_1*dypdXT - lambda_1^2*dydXT;
%dsdXT = dyppdXT - dy2rdXT;
dsatsFidXT = @(t,X,fiBd)(dsatdnu(s(t,X)/fiBd)*(1/fiBd)*dsdXT);
dUccdXT = @(t,Z,X,fiBd,mid,Z0)(bHat^-1*...
    (dUcTildXT -...
    (dKcbdXT(t,Z,X,fiBd,mid,Z0)*sat(s(t,X)/fiBd) + Kcb(t,Z,X,fiBd,mid,Z0)*dsatsFidXT(t,X,fiBd))));
DUccDZT = @(t,Z,X,fiBd,mid,Z0)(dUccdZT(t,Z,X,fiBd) + dUccdXT(t,Z,X,fiBd,mid,Z0)*dPsiDifcHatdZT(Z));
repositorio.DUccDZT = DUccDZT;

dKbdFiBd = Kfi;
dsatsFidFiBd = @(t,X,fiBd)(dsatdnu(s(t,X)/fiBd)*(-s(t,X)/fiBd^2));
dUccdFiBd = @(t,Z,X,fiBd,mid,Z0)(-bHat^-1*...
    (dKbdFiBd*sat(s(t,X)/fiBd) + Kcb(t,Z,X,fiBd,mid,Z0)*dsatsFidFiBd(t,X,fiBd)));
repositorio.dUccdFiBd = dUccdFiBd;
dUcdcdFiBd = @(t,Z,X,fiBd,mid,Z0)(-bHat^-1*(dKbdFiBd*sat(s(t,X)/fiBd) +...
    Kdcb(t,Z,X,fiBd,mid,Z0)*dsatsFidFiBd(t,X,fiBd)));
repositorio.dUcdcdFiBd = dUcdcdFiBd;
%}
%{
dXdmi = Ad;
fcdZdHatdmid = @(t,mid,Zd0)(dInvPsiDifcdXTHat(fXd(t,mid),Zd0)*dXdmi);
fdcdZdHatdmid = @(t,mid,Zd0)(dInvPsiDifdcdXTHat(fXd(t,mid),Zd0)*dXdmi);
dFcLCddmid = @(t,mid,Z0)(dFcLCdZT(fcZdHat(t,mid,Z0))*dInvPsiDifcdXTHat(fXd(t,mid),Z0)*dXdmi);
dUcTildXdT = dy3rdXT;
DUcTilcDXdT = @(t,mid,Z0)(dUcTildXdT + dUcTilcdZdT(t,mid,Z0)*dInvPsiDifcdXTHat(fXd(t,mid),Z0));
dUcTilcddmid = @(t,mid,Z0)(DUcTilcDXdT(t,mid,Z0)*dXdmi);
dKcddmid = @(t,mid,Z0)((1 - Dlc)^-1*(dFcLCddmid(t,mid,Z0) + Dlc*dUcTilcddmid(t,mid,Z0)*sign(UcTilcd(t,mid,Z0))));
repositorio.dKcddmid = dKcddmid;
dKcbdmid = @(t,mid,Z0)(-dKcddmid(t,mid,Z0));
dUccdmid = @(t,X,fiBd,mid,Z0)(-bHat^-1*(dKcbdmid(t,mid,Z0)*sat(s(t,X)/fiBd)));
repositorio.dUccdmid = dUccdmid;
dFdcLCddmid = @(t,mid,Z0)(dFdcLCdZT(fdcZdHat(t,mid,Z0))*dInvPsiDifdcdXTHat(fXd(t,mid),Z0)*dXdmi);
DUcTildcDXdT = @(t,mid,Z0)(dUcTildXdT + dUcTildcdZdT(t,mid,Z0)*dInvPsiDifdcdXTHat(fXd(t,mid),Z0));
dUcTildcddmid = @(t,mid,Z0)(DUcTildcDXdT(t,mid,Z0)*dXdmi);
dKdcddmid = @(t,mid,Z0)((1 - Dlc)^-1*(dFdcLCddmid(t,mid,Z0) + Dlc*dUcTildcddmid(t,mid,Z0)*sign(UcTildcd(t,mid,Z0))));
dKdcbdmid = @(t,mid,Z0)(-dKdcddmid(t,mid,Z0));
dUcdcdmid = @(t,X,fiBd,mid,Z0)(-bHat^-1*(dKdcbdmid(t,mid,Z0)*sat(s(t,X)/fiBd)));
repositorio.dUcdcdmid = dUcdcdmid;

dUcdcdZT = @(t,Z,X,fiBd,mid,Z0)(bHat^-1*(dUcTildcdZT(t,Z,X) - dKdcbdZT(t,Z,X,fiBd,mid,Z0)*sat(s(t,X)/fiBd)));
dUcdcdXT = @(t,Z,X,fiBd,mid,Z0)(bHat^-1*...
    (dUcTildcdXT(t,Z,X) -...
    (dKdcbdXT(t,Z,X,fiBd,mid,Z0)*sat(s(t,X)/fiBd) + Kdcb(t,Z,X,fiBd,mid,Z0)*dsatsFidXT(t,X,fiBd))));
DUcdcDZT = @(t,Z,X,fiBd,mid,Z0)(dUcdcdZT(t,Z,X,fiBd,mid,Z0) + dUcdcdXT(t,Z,X,fiBd,mid,Z0)*dPsiDifdcHatdZT(Z));
repositorio.DUcdcDZT = DUcdcDZT;

DUccuDZT = @(t,uTilc,uTilcp,fiBd,mid,Z0)(DUccDZT(t,fZ(uTilc,uTilcp),fcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.DUccuDZT = DUccuDZT;
dZducT = Au;
DUccuDucT = @(t,uTilc,uTilcp,fiBd,mid,Zd0)(DUccuDZT(t,uTilc,uTilcp,fiBd,mid,Zd0)*dZducT);
repositorio.DUccuDucT = DUccuDucT;
dZdupcT = Av;
DUccuDupcT = @(t,uTilc,uTilcp,fiBd,mid,Zd0)(DUccuDZT(t,uTilc,uTilcp,fiBd,mid,Zd0)*dZdupcT);
repositorio.DUccuDupcT = DUccuDupcT;

DUcdcuDZT = @(t,uTilc,uTilcp,fiBd,mid,Z0)(DUcdcDZT(t,fZ(uTilc,uTilcp),fdcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.DUcdcuDZT = DUcdcuDZT;
DUcdcuDucT = @(t,uTilc,uTilcp,fiBd,mid,Zd0)(DUcdcuDZT(t,uTilc,uTilcp,fiBd,mid,Zd0)*dZducT);
repositorio.DUcdcuDucT = DUcdcuDucT;
DUcdcuDupcT = @(t,uTilc,uTilcp,fiBd,mid,Zd0)(DUcdcuDZT(t,uTilc,uTilcp,fiBd,mid,Zd0)*dZdupcT);
repositorio.DUcdcuDupcT = DUcdcuDupcT;
dUccudFiBd = @(t,uTilc,uTilcp,fiBd,mid,Z0)(dUccdFiBd(t,fZ(uTilc,uTilcp),fcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.dUccudFiBd = dUccudFiBd;

dUcdcudFiBd = @(t,uTilc,uTilcp,fiBd,mid,Z0)(dUcdcdFiBd(t,fZ(uTilc,uTilcp),fdcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.dUcdcudFiBd = dUcdcudFiBd;
dUccudmid = @(t,uTilc,uTilcp,fiBd,mid,Z0)(dUccdmid(t,fcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.dUccudmid = dUccudmid;
dUcdcudmid = @(t,uTilc,uTilcp,fiBd,mid,Z0)(dUcdcdmid(t,fdcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.dUcdcudmid = dUcdcudmid;

%ExpressÃµes para a dinÃ¢mica interna:
dUcTilcMiddmid = @(t,mid,Z0)(-daChatdZT(fcZdHat(t,mid,Z0))*dInvPsiDifcdXTHat(fXd(t,mid),Z0)*dXdmi);
dUccMiddmid = @(t,mid,Z0)(bHat^-1*(dUcTilcMiddmid(t,mid,Z0)));
repositorio.dUccMiddmid = dUccMiddmid;
dUcTildcMiddmid = @(t,mid,Z0)(-daDcHatdZT(fdcZdHat(t,mid,Z0))*dInvPsiDifdcdXTHat(fXd(t,mid),Z0)*dXdmi);
dUcdcMiddmid = @(t,mid,Z0)(bHat^-1*(dUcTildcMiddmid(t,mid,Z0)));
repositorio.dUcdcMiddmid = dUcdcMiddmid;

%Outras expressÃµes Ãºteis:
fcMidpHat = @(t,mid,Zd0)(Ami*AbHat*fcZdHat(t,mid,Zd0) + Ami*Bhat*UccMid(t,mid,Zd0));
repositorio.fcMidpHat = fcMidpHat;
fcdMidpdmid = @(t,mid,Zd0)(Ami*AbHat*fcdZdHatdmid(t,mid,Zd0) + Ami*Bhat*dUccMiddmid(t,mid,Zd0));
repositorio.fcdMidpdmid = fcdMidpdmid;
fdcMidpHat = @(t,mid,Zd0)(Ami*AbHat*fdcZdHat(t,mid,Zd0) + Ami*Bhat*UcdcMid(t,mid,Zd0));
repositorio.fdcMidpHat = fdcMidpHat;
fdcdMidpdmid = @(t,mid,Zd0)(Ami*AbHat*fdcdZdHatdmid(t,mid,Zd0) + Ami*Bhat*dUcdcMiddmid(t,mid,Zd0));
repositorio.fdcdMidpdmid = fdcdMidpdmid;
yTilpd = @(t)([ypd(t); yppd(t); ypppd(t)]);
fcXpdHat = @(t,mid,Z0)(Ad*fcMidpHat(t,mid,Z0) + Bd*yTilpd(t));
repositorio.fcXpdHat = fcXpdHat;
fdcXpdHat = @(t,mid,Z0)(Ad*fdcMidpHat(t,mid,Z0) + Bd*yTilpd(t));
repositorio.fdcXpdHat = fdcXpdHat;
fcZpdHat = @(t,mid,Z0)(dInvPsiDifcdXTHat(fXd(t,mid),Z0)*fcXpdHat(t,mid,Z0));
repositorio.fcZpdHat = fcZpdHat;
fdcZpdHat = @(t,mid,Z0)(dInvPsiDifdcdXTHat(fXd(t,mid),Z0)*fdcXpdHat(t,mid,Z0));
repositorio.fdcZpdHat = fdcZpdHat;
fFippBd = @(t,fipBd,mid,Zd0)(...
    (DfiBd(fipBd)^-1)*...
    (-Kfi*fipBd + DKcDXdT(t,mid,Zd0)*fcXpdHat(t,mid,Zd0)));
repositorio.fFippBd = fFippBd;
%}
%{
%Derivada parcial de Ucc e Ucdc em relaÃ§Ã£o a t:
dy3rdt = @(t)(yppppd(t) + 2*lambda_1*ypppd(t) + lambda_1^2*yppd(t));
dUcTildt = @(t)(dy3rdt(t));
dKcdt = @(t,Z,X)((1 - Dlc)^-1*(Dlc*dUcTildt(t)*sign(UcTilc(t,Z,X))));
dFcLCddt = @(t,mid,Z0)(dFcLCdZT(fcZdHat(t,mid,Z0))*fcZpdHat(t,mid,Z0));
dy3rddt = @(t)(yppppd(t));
dUcTilcddt = @(t,mid,Z0)(dy3rddt(t) - daChatdZT(fcZdHat(t,mid,Z0))*fcZpdHat(t,mid,Z0));
dKcddt = @(t,mid,Z0)((1 - Dlc)^-1*...
    (dFcLCddt(t,mid,Z0) + Dlc*dUcTilcddt(t,mid,Z0)*sign(UcTilcd(t,mid,Z0))));
dKcbdt = @(t,Z,X,mid,Z0)(dKcdt(t,Z,X) - dKcddt(t,mid,Z0));
dy2rdt = @(t)(ypppd(t) + 2*lambda_1*yppd(t) + lambda_1^2*ypd(t));
dsdt = @(t)(-dy2rdt(t));
dsatdt = @(t,X,fiBd)(dsatdnu(s(t,X)/fiBd)*(dsdt(t)/fiBd));
dUccdt = @(t,Z,X,fiBd,mid,Z0)(bHat^-1*...
    (dUcTildt(t) -...
    (dKcbdt(t,Z,X,mid,Z0)*sat(s(t,X)/fiBd) + Kcb(t,Z,X,fiBd,mid,Z0)*dsatdt(t,X,fiBd))));
repositorio.dUccdt = dUccdt;
dUccudt = @(t,uTilc,fiBd,mid,uTilcp,Z0)(dUccdt(t,fZ(uTilc,uTilcp),fcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.dUccudt = dUccudt;

dKdcdt = @(t,Z,X)((1 - Dlc)^-1*(Dlc*dUcTildt(t)*sign(UcTildc(t,Z,X))));
dFdcLCddt = @(t,mid,Z0)(dFdcLCdZT(fdcZdHat(t,mid,Z0))*fdcZpdHat(t,mid,Z0));
dUcTildcddt = @(t,mid,Z0)(dy3rddt(t) - daDcHatdZT(fdcZdHat(t,mid,Z0))*fdcZpdHat(t,mid,Z0));
dKdcddt = @(t,mid,Z0)((1 - Dlc)^-1*(dFdcLCddt(t,mid,Z0) + Dlc*dUcTildcddt(t,mid,Z0)*sign(UcTildcd(t,mid,Z0))));
dKdcbdt = @(t,Z,X,mid,Z0)(dKdcdt(t,Z,X) - dKdcddt(t,mid,Z0));
dUcdcdt = @(t,Z,X,fiBd,mid,Z0)(bHat^-1*(dUcTildt(t) -...
    (dKdcbdt(t,Z,X,mid,Z0)*sat(s(t,X)/fiBd) + Kdcb(t,Z,X,fiBd,mid,Z0)*dsatdnu(s(t,X)/fiBd)*(dsdt(t)/fiBd))));
repositorio.dUcdcdt = dUcdcdt;
dUcdcudt = @(t,uTilb,fiBd,mid,uTilbp,Z0)(dUcdcdt(t,fZ(uTilb,uTilbp),fdcXhat(uTilb,uTilbp),fiBd,mid,Z0));
repositorio.dUcdcudt = dUcdcudt;

%Derivadas temporais de Uc (lei de controle):
teste = 0;
argumentos.testeUpcc = teste;
%Obs.: Estudo das parcelas de fUpcc:
%(1) teste = 0: teste desativado
%(2) teste = 1: teste ativado
DUccDZTZp = @(t,uTilb,fiBd,mid,uTilbp,Z0)(...
    DUccDZT(t,fZ(uTilb,uTilbp),fcXhat(uTilb,uTilbp),fiBd,mid,Z0)*...
    fcZp(fZ(uTilb,uTilbp),Uccu(t,uTilb,uTilbp,fiBd,mid,Z0)));
repositorio.DUccDZTZp = DUccDZTZp;
dUccdFibdFipBd = @(t,uTilb,fiBd,mid,uTilbp,fipBd,Z0)(...
    dUccdFiBd(t,fZ(uTilb,uTilbp),fcXhat(uTilb,uTilbp),fiBd,mid,Z0)*fipBd);
repositorio.dUccdFibdFipBd = dUccdFibdFipBd;
dUccdmidmipd = @(t,uTilb,fiBd,mid,uTilbp,mipd,Z0)(...
    dUccdmid(t,fcXhat(uTilb,uTilbp),fiBd,mid,Z0)*mipd);
repositorio.dUccdmidmipd = dUccdmidmipd;
%}

%Derivada temporal de Uc sem uTilcpp:
fUpcc = @(t,uTilc,fiBd,mid,uTilcp,fipBd,mipd,Z0)(...
    dUccudt(t,uTilc,fiBd,mid,uTilcp,Z0) + ...
    DUccDZTZp(t,uTilc,fiBd,mid,uTilcp,Z0) +...
    dUccdFibdFipBd(t,uTilc,fiBd,mid,uTilcp,fipBd,Z0) +...
    dUccdmidmipd(t,uTilc,fiBd,mid,uTilcp,mipd,Z0));
repositorio.fUpcc = fUpcc;
fUpcdc = @(t,uTilb,fiBd,mid,uTilbp,fipBd,mipd,Z0)(...
    dUcdcudt(t,uTilb,fiBd,mid,uTilbp,Z0) + ...
    DUcdcDZT(t,fZ(uTilb,uTilbp),fdcXhat(uTilb,uTilbp),fiBd,mid,Z0)*...
    fdcZp(fZ(uTilb,uTilbp),Ucdcu(t,uTilb,uTilbp,fiBd,mid,Z0)) +...
    dUcdcdFiBd(t,fZ(uTilb,uTilbp),fcXhat(uTilb,uTilbp),fiBd,mid,Z0)*fipBd +...
    dUcdcdmid(t,fcXhat(uTilb,uTilbp),fiBd,mid,Z0)*mipd);
repositorio.fUpcdc = fUpcdc;
%{
%Derivada temporal de Uc com uTilcpp:
fZp = @(uTilcp,uTilcpp)(fZ(uTilcp,uTilcpp));
DUccDZTZpInt = @(t,uTilc,uTilcp,uTilcpp,fiBd,mid,Zd0)(...
    DUccDZT(t,fZ(uTilc,uTilcp),fcXhat(uTilc,uTilcp),fiBd,mid,Zd0)*...
    fZp(uTilcp,uTilcpp));
fUpccInt = @(t,uTilc,uTilcp,uTilcpp,fiBd,fipBd,mid,mipd,Zd0)(...
    dUccudt(t,uTilc,fiBd,mid,uTilcp,Zd0) + ...
    DUccDZTZpInt(t,uTilc,uTilcp,uTilcpp,fiBd,mid,Zd0) +...
    dUccdFibdFipBd(t,uTilc,fiBd,mid,uTilcp,fipBd,Zd0) +...
    dUccdmidmipd(t,uTilc,fiBd,mid,uTilcp,mipd,Zd0));
repositorio.fUpccInt = fUpccInt;
%}
%{
%Derivadas relativas a dUccdt:
d2UccdZTdt = @(t,Z,X,fiBd)(-(1/bHat)*dKcbdZT(t,Z,X)*dsatdt(t,X,fiBd));
d2UccducTdt = @(t,Z,X,fiBd)(d2UccdZTdt(t,Z,X,fiBd)*dZducT);
d2UccuducTdt = @(t,uc,upc,fiBd)(d2UccducTdt(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd));
repositorio.d2UccuducTdt = d2UccuducTdt;
d2UccdupcTdt = @(t,Z,X,fiBd)(d2UccdZTdt(t,Z,X,fiBd)*dZdupcT);
d2UccudupcTdt = @(t,uc,upc,fiBd)(d2UccdupcTdt(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd));
repositorio.d2UccudupcTdt = d2UccudupcTdt;
dsatdfiBd = @(t,X,fiBd)(dsatdnu(s(t,X)/fiBd)*(-s(t,X)/fiBd^2));
ddfiBddKcbdtsat = @(t,Z,X,fiBd,mid,Zd0)(dKcbdt(t,Z,X,mid,Zd0)*dsatdfiBd(t,X,fiBd));
d2satdfiBddt = @(t,X,fiBd)(d2satdnu2(s(t,X)/fiBd)*(-s(t,X)/fiBd^2)*(dsdt(t)/fiBd) +...
    dsatdnu(s(t,X)/fiBd)*(-dsdt(t)/fiBd^2));
dKcbdfiBd = Kfi;
ddfiBdKcbdsatdt = @(t,Z,X,fiBd,mid,Zd0)(dKcbdfiBd*dsatdt(t,X,fiBd) +...
    Kcb(t,Z,X,fiBd,mid,Zd0)*d2satdfiBddt(t,X,fiBd));
d2UccdfiBddt = @(t,Z,X,fiBd,mid,Zd0)(-(1/bHat)*...
    (ddfiBddKcbdtsat(t,Z,X,fiBd,mid,Zd0) + ddfiBdKcbdsatdt(t,Z,X,fiBd,mid,Zd0)));
d2UccudfiBddt = @(t,uc,upc,fiBd,mid,Zd0)(d2UccdfiBddt(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,Zd0));
repositorio.d2UccudfiBddt = d2UccudfiBddt;

ddZTd2acdqTdZdq = @(Z)(fddZTd2acdqTdZdq( Z, argumentos, repositorio ));
d2FcLCdZTdZ = @(Z)(ddZTd2acdqTdZdq(Z)*sign(daC(Z)));
fcdZddmid = @(t,mid,Zd0)(dInvPsiDifcdXTHat(fXd(t,mid),Zd0)*dXdmi);
d2FcLCdmiddZd = @(t,mid,Zd0)(d2FcLCdZTdZ(fcZdHat(t,mid,Zd0))*fcdZddmid(t,mid,Zd0));
d2FcLCdmiddZdT = @(t,mid,Zd0)((d2FcLCdmiddZd(t,mid,Zd0)).');
fXpd = @(t,mipd)(Ad*mipd + Bd*yTilpd(t));
fcZpd = @(t,mid,mipd,Zd0)(dInvPsiDifcdXTHat(fXd(t,mid),Zd0)*fXpd(t,mipd));
d2FiDifcdz1dZT = @(Z)(fd2dFiDifcdz1dZT( Z, argumentos, repositorio ));
d2FiDifcdz2dZT = @(Z)(fd2dFiDifcdz2dZT( Z, argumentos, repositorio ));
d2FiDifcdz3dZT = @(Z)(zeros(4,4));
d2FiDifcdz4dZT = @(Z)(fd2dFiDifcdz4dZT( Z, argumentos, repositorio ));
fcdz1ddmid = @(t,mid,Zd0)(Id4(1,:)*fcdZddmid(t,mid,Zd0));
fcdz2ddmid = @(t,mid,Zd0)(Id4(2,:)*fcdZddmid(t,mid,Zd0));
fcdz3ddmid = @(t,mid,Zd0)(Id4(3,:)*fcdZddmid(t,mid,Zd0));
fcdz4ddmid = @(t,mid,Zd0)(Id4(4,:)*fcdZddmid(t,mid,Zd0));
d2FiDifcdmiddZdT = @(t,mid,Zd0)(d2FiDifcdz1dZT(fcZdHat(t,mid,Zd0))*fcdz1ddmid(t,mid,Zd0) +...
    d2FiDifcdz2dZT(fcZdHat(t,mid,Zd0))*fcdz2ddmid(t,mid,Zd0) +...
    d2FiDifcdz3dZT(fcZdHat(t,mid,Zd0))*fcdz3ddmid(t,mid,Zd0) +...
    d2FiDifcdz4dZT(fcZdHat(t,mid,Zd0))*fcdz4ddmid(t,mid,Zd0));
d2PsiDifcdmiddZdT = @(t,mid,Zd0)(-d2FiDifcdmiddZdT(t,mid,Zd0));
dInvdPsiDifcdZdTdmid = @(t,mid,Zd0)(...
    -dInvPsiDifcdXTHat(fXd(t,mid),Zd0)*d2PsiDifcdmiddZdT(t,mid,Zd0)*dInvPsiDifcdXTHat(fXd(t,mid),Zd0));
d2InvPsiDifcdmiddXdT = @(t,mid,Zd0)(dInvdPsiDifcdZdTdmid(t,mid,Zd0));
fcdZpddmid = @(t,mid,mipd,Zd0)(d2InvPsiDifcdmiddXdT(t,mid,Zd0)*fXpd(t,mipd));
d2FcLCddmiddt = @(t,mid,mipd,Zd0)(d2FcLCdmiddZdT(t,mid,Zd0)*fcZpd(t,mid,mipd,Zd0) + dFcLCdZdT(t,mid,Zd0)*fcdZpddmid(t,mid,mipd,Zd0));

ddZTdFcTdZb2T = @(Z)(fddZTdFcTdZb2T( Z, argumentos, repositorio ));
d2Q1cAZdZTdZ = @(Z)(fd2Q1cAZdZTdZ( Z, argumentos, repositorio ));
d2Q1FcdZTdZ = @(Z)(fd2Q1FcHatdZTdZ( Z, argumentos, repositorio ));
d2aChatdZTdZ = @(Z)(-ddZTdFcTdZb2T(Z) - d2Q1cAZdZTdZ(Z) - d2Q1FcdZTdZ(Z));
d2UcTilcddmiddt = @(t,mid,mipd,Zd0)(-...
    (daChatdZT(fcZdHat(t,mid,Zd0))*fcdZpddmid(t,mid,mipd,Zd0) +...
    fcZpd(t,mid,mipd,Zd0).'*d2aChatdZTdZ(fcZdHat(t,mid,Zd0))*fcdZddmid(t,mid,Zd0)));
d2Kcddmiddt = @(t,mid,mipd,Zd0)((1 - Dlc)^-1*...
    (d2FcLCddmiddt(t,mid,mipd,Zd0) + Dlc*d2UcTilcddmiddt(t,mid,mipd,Zd0)*sign(UcTilcd(t,mid,Zd0))));
d2Kcbdmiddt = @(t,Z,X,mid,mipd,Zd0)(-d2Kcddmiddt(t,mid,mipd,Zd0));
d2Uccdmiddt = @(t,Z,X,fiBd,mid,mipd,Zd0)(-(1/bHat)*...
    (d2Kcbdmiddt(t,Z,X,mid,mipd,Zd0)*sat(s(t,X)/fiBd) + dKcbdmid(t,mid,Zd0)*dsatdt(t,X,fiBd)));
d2Uccudmiddt = @(t,uc,upc,fiBd,mid,mipd,Zd0)(d2Uccdmiddt(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,mipd,Zd0));
repositorio.d2Uccudmiddt = d2Uccudmiddt;

ddZTdFdcTdZb2T = @(Z)(fddZTdFdcTdZb2T( Z, argumentos, repositorio ));
d2Q1dcAZdZTdZ = @(Z)(fd2Q1dcAZdZTdZ( Z, argumentos, repositorio ));
d2Q1FdcdZTdZ = @(Z)(fd2Q1FdcHatdZTdZ( Z, argumentos, repositorio ));
d2aDcHatdZTdZ = @(Z)(-ddZTdFdcTdZb2T(Z) - d2Q1dcAZdZTdZ(Z) - d2Q1FdcdZTdZ(Z));
d2FiDifdcdz1dZT = @(Z)(fd2dFiDifdcdz1dZT( Z, argumentos, repositorio ));
d2FiDifdcdz2dZT = @(Z)(fd2dFiDifdcdz2dZT( Z, argumentos, repositorio ));
d2FiDifdcdz3dZT = @(Z)(zeros(4,4));
d2FiDifdcdz4dZT = @(Z)(fd2dFiDifdcdz4dZT( Z, argumentos, repositorio ));
fdcdZddmid = @(t,mid,Zd0)(dInvPsiDifdcdXTHat(fXd(t,mid),Zd0)*dXdmi);
fdcdz1ddmid = @(t,mid,Zd0)(Id4(1,:)*fdcdZddmid(t,mid,Zd0));
fdcdz2ddmid = @(t,mid,Zd0)(Id4(2,:)*fdcdZddmid(t,mid,Zd0));
fdcdz3ddmid = @(t,mid,Zd0)(Id4(3,:)*fdcdZddmid(t,mid,Zd0));
fdcdz4ddmid = @(t,mid,Zd0)(Id4(4,:)*fdcdZddmid(t,mid,Zd0));
d2FiDifdcdmiddZdT = @(t,mid,Zd0)(d2FiDifdcdz1dZT(fdcZdHat(t,mid,Zd0))*fdcdz1ddmid(t,mid,Zd0) +...
    d2FiDifdcdz2dZT(fdcZdHat(t,mid,Zd0))*fdcdz2ddmid(t,mid,Zd0) +...
    d2FiDifdcdz3dZT(fdcZdHat(t,mid,Zd0))*fdcdz3ddmid(t,mid,Zd0) +...
    d2FiDifdcdz4dZT(fdcZdHat(t,mid,Zd0))*fdcdz4ddmid(t,mid,Zd0));
d2PsiDifdcdmiddZdT = @(t,mid,Zd0)(-d2FiDifdcdmiddZdT(t,mid,Zd0));
dInvdPsiDifdcdZdTdmid = @(t,mid,Zd0)(...
    -dInvPsiDifdcdXTHat(fXd(t,mid),Zd0)*d2PsiDifdcdmiddZdT(t,mid,Zd0)*dInvPsiDifdcdXTHat(fXd(t,mid),Zd0));
d2InvPsiDifdcdmiddXdT = @(t,mid,Zd0)(dInvdPsiDifdcdZdTdmid(t,mid,Zd0));
%fXpd = @(t,mipd)(Ad*mipd + Bd*yTilpd(t));
fdcZpd = @(t,mid,mipd,Zd0)(dInvPsiDifdcdXTHat(fXd(t,mid),Zd0)*fXpd(t,mipd));
fdcdZpddmid = @(t,mid,mipd,Zd0)(d2InvPsiDifdcdmiddXdT(t,mid,Zd0)*fXpd(t,mipd));

ddZTd2adcdqTdZdq = @(Z)(fddZTd2adcdqTdZdq( Z, argumentos, repositorio ));
d2FdcLCdZTdZ = @(Z)(ddZTd2adcdqTdZdq(Z)*sign(daC(Z)));
d2FdcLCdmiddZd = @(t,mid,Zd0)(d2FdcLCdZTdZ(fdcZdHat(t,mid,Zd0))*fdcdZddmid(t,mid,Zd0));
d2FdcLCdmiddZdT = @(t,mid,Zd0)((d2FdcLCdmiddZd(t,mid,Zd0)).');
d2FdcLCddmiddt = @(t,mid,Z0)(d2FdcLCdmiddZdT(t,mid,Zd0)*fdcZpd(t,mid,mipd,Zd0) + dFdcLCdZdT(t,mid,Zd0)*fdcdZpddmid(t,mid,mipd,Zd0));
d2UcTildcddmiddt = @(t,mid,mipd,Zd0)(-...
    (daDcHatdZT(fdcZdHat(t,mid,Zd0))*fdcdZpddmid(t,mid,mipd,Zd0) +...
    fdcZpd(t,mid,mipd,Zd0).'*d2aDcHatdZTdZ(fdcZdHat(t,mid,Zd0))*fdcdZddmid(t,mid,Zd0)));
d2Kdcddmiddt = @(t,mid,mipd,Zd0)((1 - Dlc)^-1*...
    (d2FdcLCddmiddt(t,mid,Zd0) + Dlc*d2UcTildcddmiddt(t,mid,mipd,Zd0)*sign(UcTildcd(t,mid,Zd0))));
d2Kdcbdmiddt = @(t,Z,X,mid,mipd,Zd0)(-d2Kdcddmiddt(t,mid,mipd,Zd0));
d2Ucdcdmiddt = @(t,Z,X,fiBd,mid,mipd,Zd0)(-(1/bHat)*...
    (d2Kdcbdmiddt(t,Z,X,mid,mipd,Zd0)*sat(s(t,X)/fiBd) + dKdcbdmid(t,mid,Zd0)*dsatdt(t,X,fiBd)));
d2Ucdcudmiddt = @(t,uc,upc,fiBd,mid,mipd,Zd0)(d2Ucdcdmiddt(t,fZ(uc,upc),PsiDifdcHat(fZ(uc,upc)),fiBd,mid,mipd,Zd0));
repositorio.d2Ucdcudmiddt = d2Ucdcudmiddt;

%Derivadas relativas a dUccdfiBdfipBd:
ddZTdUccdfiBdfipBd = @(t,Z,X,fiBd,fipBd)(-(1/bHat)*dKcbdZT(t,Z,X)*dsatdfiBd(t,X,fiBd)*fipBd);
dsatdXT = @(t,X,fiBd)(dsatdnu(s(t,X)/fiBd)*(dsdXT/fiBd));
d2satdXTdfiBd = @(t,X,fiBd)(d2satdnu2(s(t,X)/fiBd)*(dsdXT/fiBd)*(-s(t,X)/fiBd^2) +...
    dsatdnu(s(t,X)/fiBd)*(-dsdXT/fiBd^2));
ddXTdUccdfiBdfipBd = @(t,Z,X,fiBd,mid,fipBd,Zd0)(-(1/bHat)*...
    (dKcbdfiBd*dsatdXT(t,X,fiBd) +...
    dKcbdXT(t,Z,X,fiBd,mid,Zd0)*dsatdfiBd(t,X,fiBd) +...
    Kcb(t,Z,X,fiBd,mid,Zd0)*d2satdXTdfiBd(t,X,fiBd))*fipBd);
DDZTdUccdfiBdfipBd = @(t,Z,X,fiBd,mid,fipBd,Zd0)(...
    ddZTdUccdfiBdfipBd(t,Z,X,fiBd,fipBd) +...
    ddXTdUccdfiBdfipBd(t,Z,X,fiBd,mid,fipBd,Zd0)*dPsiDifcHatdZT(Z));
dducTdUccdfiBdfipBd = @(t,Z,X,fiBd,mid,fipBd,Zd0)(...
    DDZTdUccdfiBdfipBd(t,Z,X,fiBd,mid,fipBd,Zd0)*dZducT);
dducTdUccudfiBdfipBd = @(t,uc,upc,fiBd,mid,fipBd,Zd0)(...
    dducTdUccdfiBdfipBd(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,fipBd,Zd0));
repositorio.dducTdUccudfiBdfipBd = dducTdUccudfiBdfipBd;
ddupcTdUccdfiBdfipBd = @(t,Z,X,fiBd,mid,fipBd,Zd0)(...
    DDZTdUccdfiBdfipBd(t,Z,X,fiBd,mid,fipBd,Zd0)*dZdupcT);
ddupcTdUccudfiBdfipBd = @(t,uc,upc,fiBd,mid,fipBd,Zd0)(...
    ddupcTdUccdfiBdfipBd(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,fipBd,Zd0));
repositorio.ddupcTdUccudfiBdfipBd = ddupcTdUccudfiBdfipBd;
d2satdfiBd2 = @(t,X,fiBd)(...
    d2satdnu2(s(t,X)/fiBd)*(-s(t,X)/fiBd^2)^2 +...
    dsatdnu(s(t,X)/fiBd)*(2*s(t,X)/fiBd^3));
ddfiBddUccdfiBdfipBd = @(t,Z,X,fiBd,mid,fipBd,Zd0)(-(1/bHat)*...
    (2*dKcbdfiBd*dsatdfiBd(t,X,fiBd) + Kcb(t,Z,X,fiBd,mid,Zd0)*d2satdfiBd2(t,X,fiBd))*...
    fipBd);
ddfiBddUccudfiBdfipBd = @(t,uc,upc,fiBd,mid,fipBd,Zd0)(...
    ddfiBddUccdfiBdfipBd(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,fipBd,Zd0));
repositorio.ddfiBddUccudfiBdfipBd = ddfiBddUccudfiBdfipBd;
dUccdfiBd = @(t,Z,X,fiBd,mid,Zd0)(-(1/bHat)*...
    (dKcbdfiBd*sat(s(t,X)/fiBd) + Kcb(t,Z,X,fiBd,mid,Zd0)*dsatdfiBd(t,X,fiBd)));
ddfipBddUccdfiBdfipBd = @(t,Z,X,fiBd,mid,Zd0)(dUccdfiBd(t,Z,X,fiBd,mid,Zd0));
ddfipBddUccudfiBdfipBd = @(t,uc,upc,fiBd,mid,Zd0)(...
    ddfipBddUccdfiBdfipBd(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,Zd0));
repositorio.ddfipBddUccudfiBdfipBd = ddfipBddUccudfiBdfipBd;
ddmiddUccdfiBdfipBd = @(t,X,fiBd,mid,fipBd,Zd0)(-(1/bHat)*...
    (dKcbdmid(t,mid,Zd0)*dsatdfiBd(t,X,fiBd))*fipBd);
ddmiddUccudfiBdfipBd = @(t,uc,upc,fiBd,mid,fipBd,Zd0)(...
    ddmiddUccdfiBdfipBd(t,PsiDifcHat(fZ(uc,upc)),fiBd,mid,fipBd,Zd0));
repositorio.ddmiddUccudfiBdfipBd = ddmiddUccudfiBdfipBd;

%Derivadas relativas a dUccdmidmipd:
ddXTdUccdmidmipd = @(t,X,fiBd,mid,mipd,Zd0)(-(1/bHat)*...
    (dKcbdmid(t,mid,Zd0)*dsatdXT(t,X,fiBd))*...
    mipd);
DDZTdUccdmidmipd = @(t,Z,X,fiBd,mid,mipd,Zd0)(...
    ddXTdUccdmidmipd(t,X,fiBd,mid,mipd,Zd0)*dPsiDifcHatdZT(Z));
dducTdUccdmidmipd = @(t,Z,X,fiBd,mid,mipd,Zd0)(...
    DDZTdUccdmidmipd(t,Z,X,fiBd,mid,mipd,Zd0)*dZducT);
dducTdUccudmidmipd = @(t,uc,upc,fiBd,mid,mipd,Zd0)(...
    dducTdUccdmidmipd(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,mipd,Zd0));
repositorio.dducTdUccudmidmipd = dducTdUccudmidmipd;
ddupcTdUccdmidmipd = @(t,Z,X,fiBd,mid,mipd,Zd0)(...
    DDZTdUccdmidmipd(t,Z,X,fiBd,mid,mipd,Zd0)*dZdupcT);
ddupcTdUccudmidmipd = @(t,uc,upc,fiBd,mid,mipd,Zd0)(...
    ddupcTdUccdmidmipd(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,mipd,Zd0));
repositorio.ddupcTdUccudmidmipd = ddupcTdUccudmidmipd;
ddfiBddUccdmidmipd = @(t,X,fiBd,mid,mipd,Zd0)(-(1/bHat)*...
    (dKcbdmid(t,mid,Zd0)*dsatdfiBd(t,X,fiBd))*...
    mipd);
ddfiBddUccudmidmipd = @(t,uc,upc,fiBd,mid,mipd,Zd0)(...
    ddfiBddUccdmidmipd(t,PsiDifcHat(fZ(uc,upc)),fiBd,mid,mipd,Zd0));
repositorio.ddfiBddUccudmidmipd = ddfiBddUccudmidmipd;
ddmiddFcLCdZd = @(t,mid,Zd0)(d2FcLCdZTdZ(fcZdHat(t,mid,Zd0))*fcdZddmid(t,mid,Zd0));
ddmiddFcLCdZdT = @(t,mid,Zd0)(ddmiddFcLCdZd(t,mid,Zd0).');
fcd2Zddmid2 = @(t,mid,Zd0)(d2InvPsiDifcdmiddXdT(t,mid,Zd0)*dXdmi);
d2FcLCddmid2 = @(t,mid,Zd0)(ddmiddFcLCdZdT(t,mid,Zd0)*fcdZddmid(t,mid,Zd0) +...
    dFcLCdZdT(t,mid,Zd0)*fcd2Zddmid2(t,mid,Zd0));
d2aChatdmiddZdT = @(t,mid,Zd0)((d2aChatdZTdZ(fcZdHat(t,mid,Zd0))*fcdZddmid(t,mid,Zd0)).');
d2UcTilcddmid2 = @(t,mid,Zd0)(-...
    (d2aChatdmiddZdT(t,mid,Zd0)*dInvPsiDifcdXTHat(fXd(t,mid),Zd0) +...
    daChatdZT(fcZdHat(t,mid,Zd0))*d2InvPsiDifcdmiddXdT(t,mid,Zd0))*...
    dXdmi);
d2Kcddmid2 = @(t,mid,Zd0)((1-Dlc)^-1*...
    (d2FcLCddmid2(t,mid,Zd0) + Dlc*d2UcTilcddmid2(t,mid,Zd0)*sign(UcTilcd(t,mid,Zd0))));
d2Kcbdmid2 = @(t,mid,Zd0)(-d2Kcddmid2(t,mid,Zd0));
ddmiddUccdmidmipd = @(t,X,fiBd,mid,mipd,Zd0)(-(1/bHat)*...
    (d2Kcbdmid2(t,mid,Zd0)*sat(s(t,X)/fiBd))*mipd);
ddmiddUccudmidmipd = @(t,uc,upc,fiBd,mid,mipd,Zd0)(...
    ddmiddUccdmidmipd(t,PsiDifcHat(fZ(uc,upc)),fiBd,mid,mipd,Zd0));
repositorio.ddmiddUccudmidmipd = ddmiddUccudmidmipd;
ddmipddUccdmidmipd = @(t,X,fiBd,mid,Zd0)(dUccdmid(t,X,fiBd,mid,Zd0));
ddmipddUccudmidmipd = @(t,uc,upc,fiBd,mid,Zd0)(...
    ddmipddUccdmidmipd(t,PsiDifcHat(fZ(uc,upc)),fiBd,mid,Zd0));
repositorio.ddmipddUccudmidmipd = ddmipddUccudmidmipd;

%Derivadas relativas a dducTDUccDZTZp e ddupcTDUccDZTZp:
d2UcTilcdZTdZ = @(Z)(-d2aChatdZTdZ(Z));
d2KcdZTdZ = @(t,Z,X)((1-Dlc)^-1*...
    (d2FcLCdZTdZ(Z) + Dlc*d2UcTilcdZTdZ(Z)*sign(UcTilc(t,Z,X))));
d2KcbdZTdZ = @(t,Z,X)(d2KcdZTdZ(t,Z,X));
d2UccdZTdZ = @(t,Z,X,fiBd)((1/bHat)*...
    (d2UcTilcdZTdZ(Z) - d2KcbdZTdZ(t,Z,X)*sat(s(t,X)/fiBd)));
dFiDifcHatTdZ = @(Z)(dFiDifcHatdZT(Z).');
dsatdX = @(t,X,fiBd)(dsatdXT(t,X,fiBd).');
d2UccdZTdX = @(t,Z,X,fiBd)(-(1/bHat)*dsatdX(t,X,fiBd)*dKcbdZT(t,Z,X));
dUccdX = @(t,Z,X,fiBd,mid,Zd0)(dUccdXT(t,Z,X,fiBd,mid,Zd0).');
dUccdx1 = @(t,Z,X,fiBd,mid,Zd0)(Id4(1,:)*dUccdX(t,Z,X,fiBd,mid,Zd0));
dUccdx2 = @(t,Z,X,fiBd,mid,Zd0)(Id4(2,:)*dUccdX(t,Z,X,fiBd,mid,Zd0));
dUccdx3 = @(t,Z,X,fiBd,mid,Zd0)(Id4(3,:)*dUccdX(t,Z,X,fiBd,mid,Zd0));
dUccdx4 = @(t,Z,X,fiBd,mid,Zd0)(Id4(4,:)*dUccdX(t,Z,X,fiBd,mid,Zd0));
d2FiDifcTdZTdZdUccdX = @(t,Z,X,fiBd,mid,Zd0)(...
    d2FiDifcdz1dZT(Z)*dUccdx1(t,Z,X,fiBd,mid,Zd0) +...
    d2FiDifcdz2dZT(Z)*dUccdx2(t,Z,X,fiBd,mid,Zd0) +...
    d2FiDifcdz3dZT(Z)*dUccdx3(t,Z,X,fiBd,mid,Zd0) +...
    d2FiDifcdz4dZT(Z)*dUccdx4(t,Z,X,fiBd,mid,Zd0));
ddZTdPsiDifcTdZdUccdX = @(t,Z,X,fiBd,mid,Zd0)(-...
    (dFiDifcHatTdZ(Z)*d2UccdZTdX(t,Z,X,fiBd) +...
    d2FiDifcTdZTdZdUccdX(t,Z,X,fiBd,mid,Zd0)));
dKcbdZ = @(t,Z,X)(dKcbdZT(t,Z,X).');
d2UccdXTdZ = @(t,Z,X,fiBd)(-(1/bHat)*dKcbdZ(t,Z,X)*dsatdXT(t,X,fiBd));
dPsiDifcHatTdZ = @(Z)(dPsiDifcHatdZT(Z).');
dKcbdX = @(t,Z,X,fiBd,mid,Zd0)(dKcbdXT(t,Z,X,fiBd,mid,Zd0).');
%dsdX = dsdXT.';
d2satdXTdX = @(t,X,fiBd)(d2satdnu2(s(t,X)/fiBd)*(1/fiBd^2)*(dsdX*dsdXT));
d2UccdXTdX = @(t,Z,X,fiBd,mid,Zd0)(-(1/bHat)*...
    (dsatdX(t,X,fiBd)*dKcbdXT(t,Z,X,fiBd,mid,Zd0) +...
    dKcbdX(t,Z,X,fiBd,mid,Zd0)*dsatdXT(t,X,fiBd) +...
    Kcb(t,Z,X,fiBd,mid,Zd0)*d2satdXTdX(t,X,fiBd)));
ddXTdPsiDifcTdZdUccdX = @(t,Z,X,fiBd,mid,Zd0)(...
    dPsiDifcHatTdZ(Z)*d2UccdXTdX(t,Z,X,fiBd,mid,Zd0));
D2UccDZTDZ = @(t,Z,X,fiBd,mid,Zd0)(...
    d2UccdZTdZ(t,Z,X,fiBd) + ddZTdPsiDifcTdZdUccdX(t,Z,X,fiBd,mid,Zd0) +...
    (d2UccdXTdZ(t,Z,X,fiBd) + ddXTdPsiDifcTdZdUccdX(t,Z,X,fiBd,mid,Zd0))*...
    dPsiDifcHatdZT(Z));
dducTDUccDZTZp = @(t,Z,X,fiBd,mid,Zp,Zd0)(...
    Zp.'*D2UccDZTDZ(t,Z,X,fiBd,mid,Zd0)*dZducT);
dducTDUccuDZTZp = @(t,uc,upc,uppc,fiBd,mid,Zd0)(...
    dducTDUccDZTZp(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,fZ(upc,uppc),Zd0));
repositorio.dducTDUccuDZTZp = dducTDUccuDZTZp;
dZpdupcT = Au;
ddupcTDUccDZTZp = @(t,Z,X,fiBd,mid,Zp,Zd0)(...
    DUccDZT(t,Z,X,fiBd,mid,Zd0)*dZpdupcT +...
    Zp.'*D2UccDZTDZ(t,Z,X,fiBd,mid,Zd0)*dZdupcT);
ddupcTDUccuDZTZp = @(t,uc,upc,uppc,fiBd,mid,Zd0)(...
    ddupcTDUccDZTZp(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,fZ(upc,uppc),Zd0));
repositorio.ddupcTDUccuDZTZp = ddupcTDUccuDZTZp;
dZpduppcT = Av;
dduppcTDUccDZTZp = @(t,Z,X,fiBd,mid,Zd0)(...
    DUccDZT(t,Z,X,fiBd,mid,Zd0)*dZpduppcT);
dduppcTDUccuDZTZp = @(t,uc,upc,fiBd,mid,Zd0)(...
    dduppcTDUccDZTZp(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,Zd0));
repositorio.dduppcTDUccuDZTZp = dduppcTDUccuDZTZp;

%Derivadas relativas a ddfiBdDUccDZTZp:
d2UccdfiBddZT = @(t,Z,X,fiBd)(-(1/bHat)*dKcbdZT(t,Z,X)*dsatdfiBd(t,X,fiBd));
ddfiBddKcbdXTsat = @(t,Z,X,fiBd,mid,Zd0)(...
    dKcbdXT(t,Z,X,fiBd,mid,Zd0)*dsatdfiBd(t,X,fiBd));
d2satdfiBddXT = @(t,X,fiBd)(...
    (d2satdnu2(s(t,X)/fiBd)*(-s(t,X)/fiBd^2)*(1/fiBd) +...
    dsatdnu(s(t,X)/fiBd)*(-1/fiBd^2))*...
    dsdXT);
ddfiBdKcbdsatdXT = @(t,Z,X,fiBd,mid,Zd0)(dKcbdfiBd*dsatdXT(t,X,fiBd) +...
    Kcb(t,Z,X,fiBd,mid,Zd0)*d2satdfiBddXT(t,X,fiBd));
d2UccdfiBddXT = @(t,Z,X,fiBd,mid,Zd0)(-(1/bHat)*...
    (ddfiBddKcbdXTsat(t,Z,X,fiBd,mid,Zd0) +...
    ddfiBdKcbdsatdXT(t,Z,X,fiBd,mid,Zd0)));
ddfiBdDUccDZTZp = @(t,Z,X,fiBd,mid,Zp,Zd0)((d2UccdfiBddZT(t,Z,X,fiBd) +...
    d2UccdfiBddXT(t,Z,X,fiBd,mid,Zd0)*dPsiDifcHatdZT(Z))*Zp);
ddfiBdDUccuDZTZp = @(t,uc,upc,uppc,fiBd,mid,Zd0)(...
    ddfiBdDUccDZTZp(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,fZ(upc,uppc),Zd0));
repositorio.ddfiBdDUccuDZTZp = ddfiBdDUccuDZTZp;

%Derivadas relativas a ddmidDUccDZTZp:
ddmidKcbdsatdXT = @(t,X,fiBd,mid,Zd0)(dKcbdmid(t,mid,Zd0)*dsatdXT(t,X,fiBd));
d2UccdmiddXT = @(t,X,fiBd,mid,Zd0)(-(1/bHat)*...
    (ddmidKcbdsatdXT(t,X,fiBd,mid,Zd0)));
ddmidDUccDZTZp = @(t,Z,X,fiBd,mid,Zp,Zd0)(...
    d2UccdmiddXT(t,X,fiBd,mid,Zd0)*dPsiDifcHatdZT(Z)*Zp);
ddmidDUccuDZTZp = @(t,uc,upc,uppc,fiBd,mid,Zd0)(...
    ddmidDUccDZTZp(t,fZ(uc,upc),PsiDifcHat(fZ(uc,upc)),fiBd,mid,fZ(upc,uppc),Zd0));
repositorio.ddmidDUccuDZTZp = ddmidDUccuDZTZp;

%Derivadas de Upcc:
dUpccuducT = @(t,uc,upc,uppc,fiBd,mid,fipBd,mipd,Zd0)(...
    d2UccuducTdt(t,uc,upc,fiBd) +...
    dducTDUccuDZTZp(t,uc,upc,uppc,fiBd,mid,Zd0) +...
    dducTdUccudfiBdfipBd(t,uc,upc,fiBd,mid,fipBd,Zd0) +...
    dducTdUccudmidmipd(t,uc,upc,fiBd,mid,mipd,Zd0));
repositorio.dUpccuducT = dUpccuducT;
dUpccudupcT = @(t,uc,upc,uppc,fiBd,mid,fipBd,mipd,Zd0)(...
    d2UccudupcTdt(t,uc,upc,fiBd) +...
    ddupcTDUccuDZTZp(t,uc,upc,uppc,fiBd,mid,Zd0) +...
    ddupcTdUccudfiBdfipBd(t,uc,upc,fiBd,mid,fipBd,Zd0) +...
    ddupcTdUccudmidmipd(t,uc,upc,fiBd,mid,mipd,Zd0));
repositorio.dUpccudupcT = dUpccudupcT;
dUpccuduppcT = @(t,uc,upc,fiBd,mid,Zd0)(...
    dduppcTDUccuDZTZp(t,uc,upc,fiBd,mid,Zd0));
repositorio.dUpccuduppcT = dUpccuduppcT;
dUpccudfiBd = @(t,uc,upc,uppc,fiBd,mid,fipBd,mipd,Zd0)(...
    d2UccudfiBddt(t,uc,upc,fiBd,mid,Zd0) +...
    ddfiBdDUccuDZTZp(t,uc,upc,uppc,fiBd,mid,Zd0) +...
    ddfiBddUccudfiBdfipBd(t,uc,upc,fiBd,mid,fipBd,Zd0) +...
    ddfiBddUccudmidmipd(t,uc,upc,fiBd,mid,mipd,Zd0));
repositorio.dUpccudfiBd = dUpccudfiBd;
dUpccudfipBd = @(t,uc,upc,fiBd,mid,Zd0)(...
    ddfipBddUccudfiBdfipBd(t,uc,upc,fiBd,mid,Zd0));
repositorio.dUpccudfipBd = dUpccudfipBd;
dUpccudmid = @(t,uc,upc,uppc,fiBd,mid,fipBd,mipd,Zd0)(...
    d2Uccudmiddt(t,uc,upc,fiBd,mid,mipd,Zd0) +...
    ddmidDUccuDZTZp(t,uc,upc,uppc,fiBd,mid,Zd0) +...
    ddmiddUccudfiBdfipBd(t,uc,upc,fiBd,mid,fipBd,Zd0) +...
    ddmiddUccudmidmipd(t,uc,upc,fiBd,mid,mipd,Zd0));
repositorio.dUpccudmid = dUpccudmid;
dUpccudmipd = @(t,uc,upc,fiBd,mid,Zd0)(...
    ddmipddUccudmidmipd(t,uc,upc,fiBd,mid,Zd0));
repositorio.dUpccudmipd = dUpccudmipd;
%}

%%
%MODELOS DA PLANTA

%Modelo sem acoplamento das dinâmicas lateral e longitudinal:
Mt = diag([I2b + dm*l2^2, I1]);
Kt = [K10 + K20, -K20; -K20, K20];
Ct = diag([c2, c3]);
Bt = [K10; 0];

Mu = diag([M2b + dm, M1]);
Ku = [k1 + k2, -k2; -k2, k2];
Cu = diag([b2, b3]);
Bu = [k1/betaS; 0];

%Força normal:
Normal = @(u3)(k3*(u3 - uC)*double(u3>=uC));
%Normal = @(u3)(k3*(u3 - uC)*degrauS(u3 - uC));
Normalb = @(u3b)(Normal(u3b + u3_Eq));
repositorio.Normalb = Normalb;
dNormaldu3 = @(u3)(k3*double(u3>=uC));
dNormalbdu3b = @(u3b)(dNormaldu3(u3b + u3_Eq));
repositorio.dNormalbdu3b = dNormalbdu3b;
%Formas vetorizadas:
NormalVt = @(u3)(k3*(u3 - uC*ones(size(u3))).*double(u3>=uC*ones(size(u3))));
NormalbVt = @(u3b)(NormalVt(u3b + u3_Eq*ones(size(u3b))));
repositorio.NormalbVt = NormalbVt;

%Coeficiente de atrito:
%(a) Modelo contínuo:
mic0 = argumentos.mic0;
w0 = argumentos.w0;
lambdami = argumentos.lambdami;
nmi = argumentos.nmi;
fiAtr = @(nu)(tanh(pi*nu) +...
    lambdami*((2*nmi)/(2*nmi - 1))*nu*(1 + (1/(2*nmi - 1))*nu^(2*nmi))^-1);
repositorio.fiAtr = fiAtr;
mi3c = @(teta3p)(mic0*fiAtr(teta3p/w0));
repositorio.mi3c = mi3c;
dfiAtrdnu = @(nu)(-pi*(tanh(pi*nu)^2 - 1)  +...
    (2*lambdami*nmi)/(2*nmi + nu^(2*nmi) - 1) -...
    (4*lambdami*nmi^2*nu^(2*nmi))/(2*nmi + nu^(2*nmi) - 1)^2);
repositorio.dfiAtrdnu = dfiAtrdnu;
dmi3cdw3 = @(teta3p)(mic0*(1/w0)*dfiAtrdnu(teta3p/w0));
repositorio.dmi3cdw3 = dmi3cdw3;
%Formas vetorizadas:
fiAtrVt = @(nu)(tanh(pi*nu) +...
    lambdami*((2*nmi)/(2*nmi - 1))*nu.*(ones(size(nu)) + (1/(2*nmi - 1))*nu.^(2*nmi)).^-1);
mi3cVt = @(teta3p)(mic0*fiAtrVt(teta3p/w0));
repositorio.mi3cVt = mi3cVt;

%(b) Modelo descontínuo:
mi30dc = argumentos.mi0;
mi31dc = argumentos.mi1;
mi32dc = argumentos.mi2;
beta1mi3dc = argumentos.beta1mi;
beta2mi3dc = argumentos.beta2mi;
opcao_sign = argumentos.opcao_sign;
switch opcao_sign
    case 1
        mi3dc = @(teta3p)(sign(teta3p)*...
            (mi30dc +...
            mi31dc*(1 - 2*(1 + exp(beta1mi3dc*abs(teta3p)))^-1) -...
            mi32dc*(1 - 2*(1 + exp(beta2mi3dc*abs(teta3p)))^-1)));
        %Forma vetorizada:
        mi3dcVt = @(teta3p)(sign(teta3p).*...
            (mi30dc*ones(size(teta3p)) +...
            mi31dc*(ones(size(teta3p)) - 2*(ones(size(teta3p)) + exp(beta1mi3dc*abs(teta3p))).^-1) -...
            mi32dc*(ones(size(teta3p)) - 2*(ones(size(teta3p)) + exp(beta2mi3dc*abs(teta3p))).^-1)));
    case 2
        mi3dc = @(teta3p)(signS(teta3p)*...
            (mi30dc +...
            mi31dc*(1 - 2*(1 + exp(beta1mi3dc*abs(teta3p)))^-1) -...
            mi32dc*(1 - 2*(1 + exp(beta2mi3dc*abs(teta3p)))^-1)));
        %Forma vetorizada:
        mi3dcVt = @(teta3p)(signS(teta3p).*...
            (mi30dc*ones(size(teta3p)) +...
            mi31dc*(ones(size(teta3p)) - 2*(ones(size(teta3p)) + exp(beta1mi3dc*abs(teta3p))).^-1) -...
            mi32dc*(ones(size(teta3p)) - 2*(ones(size(teta3p)) + exp(beta2mi3dc*abs(teta3p))).^-1)));
    otherwise
        disp('other value')
end
repositorio.mi3dc = mi3dc;
repositorio.mi3dcVt = mi3dcVt;

%Torque de atrito:
%(a) Modelo contínuo:
Tatr3c = @(teta3p,u3b)(mi3c(teta3p)*Normalb(u3b)*R3);
repositorio.Tatr3c = Tatr3c;
%(b) Modelo descontínuo:
Tatr3dc = @(teta3p,u3b)(mi3dc(teta3p)*Normalb(u3b)*R3);
repositorio.Tatr3dc = Tatr3dc;
%Formas vetorizadas:
Tatr3cVt = @(teta3p,u3b)(mi3cVt(teta3p).*NormalbVt(u3b)*R3);
repositorio.Tatr3cVt = Tatr3cVt;
Tatr3dcVt = @(teta3p,u3b)(mi3dcVt(teta3p).*NormalbVt(u3b)*R3);
repositorio.Tatr3dcVt = Tatr3dcVt;

%%
%SISTEMA SIMPLIFICADO, COM ACOPLAMENTO DAS DINÂMICAS LONGITUDINAL E LATERAL
%Obs.: uTilAc = [teta2; yL; zL; u2; u3; teta3];

alfam2 = argumentos.alfam2;
Id6 = eye(6);
fteta2 = repositorio.fteta2;
fyL = repositorio.fyL;
fzL = repositorio.fzL;
fu2 = repositorio.fu2;
fu3 = repositorio.fu3;
fteta3 = repositorio.fteta3;
%Vetor posição de equilíbrio:
uTilAc_Eq = argumentos.uTilAc_Eq;

M1ac = [      I2b + l2^2*dm, -dm*l2*sin(alfam2), dm*l2*cos(alfam2),        0,  0,  0;
         -dm*l2*sin(alfam2),           M2b + dm,                 0,        0,  0,  0;
          dm*l2*cos(alfam2),                  0,          M2b + dm,        0,  0,  0;
                          0,                  0,                 0, M2b + dm,  0,  0;
                          0,                  0,                 0,        0, M1,  0;
                          0,                  0,                 0,        0,  0, I1];
argumentos.M1ac = M1ac;

fM2ac = @(uTilAc)([                                                 0, -dm*l2*(sin(fteta2(uTilAc) + alfam2) - sin(alfam2)), dm*l2*(cos(fteta2(uTilAc) + alfam2) - cos(alfam2)), 0, 0, 0;
                  -dm*l2*(sin(fteta2(uTilAc) + alfam2) - sin(alfam2)),                                                   0,                                                  0, 0, 0, 0;
                   dm*l2*(cos(fteta2(uTilAc) + alfam2) - cos(alfam2)),                                                   0,                                                  0, 0, 0, 0;
                                                                    0,                                                   0,                                                  0, 0, 0, 0;
                                                                    0,                                                   0,                                                  0, 0, 0, 0;
                                                                    0,                                                   0,                                                  0, 0, 0, 0]);
fM2bAc = @(uTilbAc)(fM2ac(uTilbAc + uTilAc_Eq));
repositorio.fM2bAc = fM2bAc;

dM2uppduT = @(uTilAc,uTilppAc)([ -dm*l2*fyL(uTilppAc)*cos(alfam2 + fteta2(uTilAc)) - dm*l2*fzL(uTilppAc)*sin(alfam2 + fteta2(uTilAc)), 0, 0, 0, 0, 0;
    -dm*l2*fteta2(uTilppAc)*cos(alfam2 + fteta2(uTilAc)), 0, 0, 0, 0, 0;
    -dm*l2*fteta2(uTilppAc)*sin(alfam2 + fteta2(uTilAc)), 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0]);
dM2bubppdubT = @(uTilbAc,uTilbppAc)(dM2uppduT(uTilbAc + uTilAc_Eq,uTilbppAc));
repositorio.dM2bubppdubT = dM2bubppdubT;

Kac = [K10 + K20, 0, 0,       0,   0, -K20;
               0, 0, 0,       0,   0,    0;
               0, 0, 0,       0,   0,    0;
               0, 0, 0, k1 + k2, -k2,    0;
               0, 0, 0,     -k2,  k2,    0;
            -K20, 0, 0,       0,   0,  K20];
argumentos.Kac = Kac;

Wac = [0; 0; 0; (M2b + dm)*g; M1*g; 0];

L2 = L - L1;
FuAc = @(uTilAc)([0;
    2*fyL(uTilAc)*fu2(uTilAc)*((k1/L1) - (k2/L2)) + 2*fyL(uTilAc)*fu3(uTilAc)*(k2/L2) + 2*fyL(uTilAc)*(fyL(uTilAc)^2 + fzL(uTilAc)^2)*((k1/L1^2) + (k2/L2^2)); 
    2*fzL(uTilAc)*fu2(uTilAc)*((k1/L1) - (k2/L2)) + 2*fzL(uTilAc)*fu3(uTilAc)*(k2/L2) + 2*fzL(uTilAc)*(fyL(uTilAc)^2 + fzL(uTilAc)^2)*((k1/L1^2) + (k2/L2^2));
    ((k1/L1) - (k2/L2))*(fyL(uTilAc)^2 + fzL(uTilAc)^2);
    (k2/L2)*(fyL(uTilAc)^2 + fzL(uTilAc)^2);
    0]);
FbuAc = @(uTilbAc)(Kac*uTilAc_Eq + FuAc(uTilbAc + uTilAc_Eq) - Wac);
repositorio.FbuAc = FbuAc;
dFuduTAc = @(uTilAc)([0, 0, 0, 0, 0, 0;
    0, 2*fu2(uTilAc)*((k1/L1) - (k2/L2)) + 2*fu3(uTilAc)*(k2/L2) + (6*fyL(uTilAc)^2 + 2*fzL(uTilAc)^2)*((k1/L1^2) + (k2/L2^2)), 4*fyL(uTilAc)*fzL(uTilAc)*((k1/L1^2) + (k2/L2^2)), 2*fyL(uTilAc)*((k1/L1) - (k2/L2)), 2*fyL(uTilAc)*(k2/L2), 0;
    0, 4*fyL(uTilAc)*fzL(uTilAc)*((k1/L1^2) + (k2/L2^2)), 2*fu2(uTilAc)*((k1/L1) - (k2/L2)) + 2*fu3(uTilAc)*(k2/L2) + (2*fyL(uTilAc)^2 + 6*fzL(uTilAc)^2)*((k1/L1^2) + (k2/L2^2)), 2*fzL(uTilAc)*((k1/L1) - (k2/L2)), 2*fzL(uTilAc)*(k2/L2), 0;
    0, 2*fyL(uTilAc)*((k1/L1) - (k2/L2)), 2*fzL(uTilAc)*((k1/L1) - (k2/L2)), 0, 0, 0;
    0, 2*fyL(uTilAc)*(k2/L2), 2*fzL(uTilAc)*(k2/L2), 0, 0, 0;
    0, 0, 0, 0, 0, 0]);
dFbudubTAc = @(uTilbAc)(dFuduTAc(uTilbAc + uTilAc_Eq));
repositorio.dFbudubTAc = dFbudubTAc;

%Matriz de amortecimento:
uTilbAc0 = zeros(6,1);
KbAc = Kac + dFbudubTAc(uTilbAc0);
%Frequências naturais (rad/s):
freqNatKbw2 = sort(real(eig(M1ac\KbAc)));
freqNatKbw2(1,1) = 0;
freqNatKbw = sqrt(freqNatKbw2);
%Frequências naturais (Hz):
freqNatWTab = [sort(real(eig(M1ac\KbAc))),sort(real(eig(M1ac\Kac)))];
freqNatWTab = (2*pi)^-1*sqrt(freqNatWTab(2:6,:));
argumentos.freqNatWTab = freqNatWTab;
%Matriz de amortecimento:
%Obs: 
%(1) C = alfa*M1 + beta*K
%(2) alfa e beta foram escolhidos de modo que todos os modos fossem subamort:
%   (a) alfa = 2*(wMin*wMax)/(wMin + wMax)
%   (b) beta = 2/(wMin + wMax)
%   (c) wMin e wMax são tais que f(wMin) = f(wMax) = 1
%   (d) f(w) = (1/2)*(alfa/w + beta*w)
%   (e) qsiMin: menor coeficiente de amortecimento nominal
%(3) Parâmetros para determinação de alfa e beta: qsiMin, wRmin, wRmax
%   (a) wRmin: menor frequência natural do sistema
%   (b) wRmax: maior frequência natural do sistema
%   (c) wMin < wRmin <= w <= wRmax < wMax
%(4) Sequência de cálculo:
%   (a) determinação de qsiMin, tal que 0 <= qsiMin <= qsiMinMax
%   (b) deltaQsi = qsi*(1 + sqrt(1 - qsi^2))^-1
%   (c) deltaQsi0 = wRmin/wRmax
%   (d) deltaQsiMax = sqrt(wRmin/wRmax)
%   (e) wQsi: wRmin <= wQsi <= wRmax, se 0 <= qsiMin < qsiMin0
%       wQsi: wRmax*deltaQsi <= wQsi <= wRmin/deltaQsi, se qsiMin0 <= qsiMin <= qsiMinMax
wRmin = freqNatKbw(2,1); %rad/s
wRmax = max(freqNatKbw); %rad/s
deltaQsi0 = wRmin/wRmax;
deltaQsiMax = sqrt(wRmin/wRmax);
qsiMax = 0.01; %Valor máximo dos coeficientes de amortecimento
argumentos.qsiMax = qsiMax;
qsiMin0 = 2*qsiMax*(deltaQsi0 + deltaQsi0^-1)^-1;
qsiMinMax = 2*qsiMax*(deltaQsiMax + deltaQsiMax^-1)^-1;
sigQsi = 0.001; %Condição: 0 < sigQsi <= 1
%{
antes_de_QsiMin0 = 1;
qsiMin = sigQsi*qsiMin0*double(antes_de_QsiMin0) +...
    ((1 - sigQsi)*qsiMin0 + sigQsi*qsiMinMax)*double(~antes_de_QsiMin0);
%Obs.: Escolha de qsiMin, especificando previamente se menor ou maior que
%qsiMin0
%}
%
qsiMin = sigQsi*qsiMinMax;
%Obs.: Escolha de qsiMin, especificando apenas um valor entre 0 e qsiMinMax
%}
argumentos.qsiMin = qsiMin;
deltaQsi = qsiMin*(qsiMax + sqrt(qsiMax^2 - qsiMin^2))^-1;
sigWqsi = 0.1;
wQsiMin = (wRmin*(1-sigWqsi) + wRmax*sigWqsi)*...
    double(qsiMin<qsiMin0) +...
    (deltaQsi*wRmax*(1-sigWqsi) + (wRmin/deltaQsi)*sigWqsi)*...
    double(qsiMin>=qsiMin0);
argumentos.wQsiMin = wQsiMin; %Freq nominal associada a qsiMin
wMin = deltaQsi*wQsiMin;
wMax = wQsiMin/deltaQsi;
alfaAc = 2*qsiMax*(wMin*wMax)/(wMin + wMax);
betaAc = 2*qsiMax/(wMin + wMax);
deltaMin = log10(wRmin) - log10(wMin);
argumentos.deltaMin = deltaMin;
deltaMax = log10(wMax) - log10(wRmax);
argumentos.deltaMax = deltaMax;
Cestr = alfaAc*M1ac + betaAc*KbAc;
wNatTil = 2*pi*freqNatWTab(:,1);
psiAmort = (1/2)*(alfaAc./wNatTil + betaAc*wNatTil);
argumentos.psiAmort = psiAmort;
Cconc = diag([c2, 0, 0, b2, b3, c3]);
adCac = argumentos.adCac;
Cac = adCac*(Cestr + Cconc);
argumentos.Cac = Cac;

GuupAc = @(uTilAc,uTilpAc)([0;
    -dm*fteta2(uTilpAc)^2*l2*cos(fteta2(uTilAc) + alfam2); 
    -dm*fteta2(uTilpAc)^2*l2*sin(fteta2(uTilAc) + alfam2);
    0;
    0;
    0]);
GubupbAc = @(uTilbAc,uTilbpAc)(GuupAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.GubupbAc = GubupbAc;
dGuupduTAc = @(uTilAc,uTilpAc)([0, 0, 0, 0, 0, 0;
    dm*l2*fteta2(uTilpAc)^2*sin(alfam2 + fteta2(uTilAc)), 0, 0, 0, 0, 0;
    -dm*l2*fteta2(uTilpAc)^2*cos(alfam2 + fteta2(uTilAc)), 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0]);
dGubupbdubTAc = @(uTilbAc,uTilbpAc)(dGuupduTAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.dGubupbdubTAc = dGubupbdubTAc;
dGuupdupTAc = @(uTilAc,uTilpAc)([0, 0, 0, 0, 0, 0;
    -2*dm*l2*fteta2(uTilpAc)*cos(alfam2 + fteta2(uTilAc)), 0, 0, 0, 0, 0;
    -2*dm*l2*fteta2(uTilpAc)*sin(alfam2 + fteta2(uTilAc)), 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0]);
dGubupbdubpTAc = @(uTilbAc,uTilbpAc)(dGuupdupTAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.dGubupbdubpTAc = dGubupbdubpTAc;

TcAc = @(uTilbAc,uTilbpAc)([0;
    0;
    0;
    0;
    Normalb(fu3(uTilbAc));
    Tatr3c(fteta3(uTilbpAc),fu3(uTilbAc))]);
repositorio.TcAc = TcAc;
TdcAc = @(uTilbAc,uTilbpAc)([0;
    0;
    0;
    0;
    Normalb(fu3(uTilbAc));
    Tatr3dc(fteta3(uTilbpAc),fu3(uTilbAc))]);
repositorio.TdcAc = TdcAc;

du3bdubTAc = Id6(5,:);
dNormalbdubTAc = @(uTilbAc)(dNormalbdu3b(fu3(uTilbAc))*du3bdubTAc);

dTatr3cdubTAc = @(uTilbAc,uTilbpAc)(mi3c(fteta3(uTilbpAc))*dNormalbdubTAc(uTilbAc)*R3);
dTcdubTAc = @(uTilbAc,uTilbpAc)([zeros(1,6);
    zeros(1,6);
    zeros(1,6);
    zeros(1,6);
    dNormalbdubTAc(uTilbAc);
    dTatr3cdubTAc(uTilbAc,uTilbpAc)]);
repositorio.dTcdubTAc = dTcdubTAc;
dTatr3dcdubTAc = @(uTilbAc,uTilbpAc)(mi3dc(fteta3(uTilbpAc))*dNormalbdubTAc(uTilbAc)*R3);
dTdcdubTAc = @(uTilbAc,uTilbpAc)([zeros(1,6);
    zeros(1,6)
    zeros(1,6);
    zeros(1,6);
    dNormalbdubTAc(uTilbAc);
    dTatr3dcdubTAc(uTilbAc,uTilbpAc)]);
repositorio.dTdcdubTAc = dTdcdubTAc;

dteta3pduTilbpT = Id6(6,:);
dmi3cduTilbpTAc = @(uTilbpAc)(dmi3cdw3(fteta3(uTilbpAc))*dteta3pduTilbpT);
dTatr3cdubpTAc = @(uTilbAc,uTilbpAc)(dmi3cduTilbpTAc(uTilbpAc)*Normalb(fu3(uTilbAc))*R3);
dTcdubpTAc = @(uTilbAc,uTilbpAc)([zeros(1,6);
    zeros(1,6);
    zeros(1,6);
    zeros(1,6);
    zeros(1,6);
    dTatr3cdubpTAc(uTilbAc,uTilbpAc)]);
repositorio.dTcdubpTAc = dTcdubpTAc;
dmi3dcduTilbpTAc = @(uTilbpAc)(dmi3dcdw3(fteta3(uTilbpAc))*dteta3pduTilbpT);
dTatr3dcdubpTAc = @(uTilbAc,uTilbpAc)(dmi3dcduTilbpTAc(uTilbpAc)*Normalb(fu3(uTilbAc))*R3);
dTdcdubpTAc = @(uTilbAc,uTilbpAc)([zeros(1,6);
    zeros(1,6);
    zeros(1,6);
    zeros(1,6);
    zeros(1,6);
    dTatr3dcdubpTAc(uTilbAc,uTilbpAc)]);
repositorio.dTdcdubpTAc = dTdcdubpTAc;

%{
%Rubbing:
Kip = argumentos.Kip;
R2 = Dr2/2;
fg = argumentos.fg;
delta0 = argumentos.delta0;
n = (1/pi)*(1 - fg/delta0)^-1;
fuL = @(yL,zL)(sqrt(yL^2 + zL^2)*double(sqrt(yL^2 + zL^2)>eps) +...
    eps*double(sqrt(yL^2 + zL^2)<=eps));
repositorio.fuL = fuL;
switch opcao_sign
    case 1
        fFn = @(uL)(Kip*(uL - fg)*double(uL>=fg));
    case 2
        fFn = @(uL)(Kip*(uL - fg)*(tanh(n*pi*(uL/delta0 - 1)) + 1));
    otherwise
        disp('other value')
end
repositorio.fFn = fFn;
fFy = @(yL,zL)(fFn(fuL(yL,zL))*(yL/fuL(yL,zL)));
repositorio.fFyRub = fFy;
fFz = @(yL,zL)(fFn(fuL(yL,zL))*(zL/fuL(yL,zL)));
repositorio.fFzRub = fFz;
Alc = @(yL,zL)(...
    [1 + (R2/fuL(yL,zL))*(zL/fuL(yL,zL))^2, -yL*zL*R2/fuL(yL,zL)^3; 
    -zL*yL*R2/fuL(yL,zL)^3, 1 + (R2/fuL(yL,zL))*(yL/fuL(yL,zL))^2]);
AvCfi = @(yL,zL)([R2, zeros(1,2); zeros(2,1), Alc(yL,zL)]);
%{
fuCp = @(yL,zL,uL,yLp,zLp)(Alc(yL,zL,uL)*[yLp;zLp]);
fuTilHatfiL = @(fiL)([-sin(fiL); cos(fiL)]);
ffiL = @(yL,zL)(atan2(zL,yL));
fvCfi = @(yL,zL,teta2p,yLp,zLp)(fuTilHatfiL(ffiL(yL,zL)).'*...
    fuCp(yL,zL,fuL(yL,zL),yLp,zLp) +...
    R2*teta2p);
%}

%Versor unitário na direção tangencial:
fuTilHatfiL = @(yL,zL)([-zL/fuL(yL,zL); yL/fuL(yL,zL)]);

fvCfi = @(yL,zL,teta2p,yLp,zLp)([1, fuTilHatfiL(yL,zL).']*AvCfi(yL,zL)*[teta2p; yLp; zLp]);
repositorio.fvCfi = fvCfi;
miIp0c = argumentos.miIp0c;
v0Ipc = argumentos.v0Ipc;
nmiIpc = argumentos.nmiIpc;
lambdamiIpc = argumentos.lambdamiIpc;
fiIpc = @(nu)(tanh(pi*nu) +...
    lambdamiIpc*((2*nmiIpc)/(2*nmiIpc - 1))*...
    (1 + (1/(2*nmiIpc - 1))*nu^(2*nmiIpc))^-1);
miIpc = @(vCfi)(miIp0c*fiIpc(vCfi/v0Ipc));
repositorio.miIpc = miIpc;
switch opcao_sign
    case 1
        miIpdc = @(vCfi)(sign(vCfi)*...
            (miIp0dc +...
            miIp1dc*(1 - 2*(1 + exp(beta1miIpdc*abs(vCfi)))^-1) +...
            miIp2dc*(1 - 2*(1 + exp(beta2miIpdc*abs(vCfi)))^-1)));
    case 2
        miIpdc = @(vCfi)(signS(vCfi)*...
            (miIp0dc +...
            miIp1dc*(1 - 2*(1 + exp(beta1miIpdc*abs(vCfi)))^-1) +...
            miIp2dc*(1 - 2*(1 + exp(beta2miIpdc*abs(vCfi)))^-1)));
    otherwise
        disp('other value')
end
fT2c = @(yL,zL,teta2p,yLp,zLp)(miIpc(fvCfi(yL,zL,teta2p,yLp,zLp))*...
    fFn(fuL(yL,zL))*(R2 + fg));
repositorio.fT2cRub = fT2c;
fT2dc = @(yL,zL,teta2p,yLp,zLp)(miIpdc(fvCfi(yL,zL,teta2p,yLp,zLp))*...
    fFn(fuL(yL,zL))*(R2 + fg));
repositorio.fT2dcRub = fT2dc;
TrubcAc = @(uTilAc,uTilpAc)(...
    [fT2c(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc));
    fFy(fyL(uTilAc),fzL(uTilAc));
    fFz(fyL(uTilAc),fzL(uTilAc));
    0;
    0;
    0]);
TrubbcAc = @(uTilbAc,uTilbpAc)(TrubcAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.TrubbcAc = TrubbcAc;
TrubdcAc = @(uTilAc,uTilpAc)(...
    [fT2dc(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc));
    fFy(fyL(uTilAc),fzL(uTilAc));
    fFz(fyL(uTilAc),fzL(uTilAc));
    0;
    0;
    0]);
TrubbdcAc = @(uTilbAc,uTilbpAc)(TrubdcAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.TrubbdcAc = TrubbdcAc;

%(a) Rubbing (derivadas):
switch opcao_sign
    case 1
        fdFnduL = @(uL)(Kip*double(uL>=fg));
    case 2
        fdFnduL = @(uL)(Kip*(tanh(n*pi*(uL/delta0 - 1)) + 1) +...
            inc_dsignS*Kip*(fg - uL)*(pi*n/delta0)*(tanh(pi*n*(uL/delta0 - 1))^2 - 1));
    otherwise
        disp('other value')
end
repositorio.dFnduL = fdFnduL;
duLdyL = @(yL,zL)(yL/fuL(yL,zL));
repositorio.duLdyL = duLdyL;
duLdzL = @(yL,zL)(zL/fuL(yL,zL));
repositorio.duLdzL = duLdzL;
dyLduT = Id6(2,:);
dzLduT = Id6(3,:);
duLduT = @(uTilAc)(duLdyL(fyL(uTilAc),fzL(uTilAc))*dyLduT +...
    duLdzL(fyL(uTilAc),fzL(uTilAc))*dzLduT);
dFnduT = @(uTilAc)(fdFnduL(fuL(fyL(uTilAc),fzL(uTilAc)))*duLduT(uTilAc));
dfiIpcdvCfi = @(nu)(- pi*(tanh(pi*nu)^2 - 1) -...
    (4*lambdamiIpc*nmiIpc^2*nu^(2*nmiIpc - 1))/(2*nmiIpc + nu^(2*nmiIpc) - 1)^2);
dmiIpcdvCfi = @(vCfi)(miIp0c*(1/v0Ipc)*dfiIpcdvCfi(vCfi/v0Ipc));
repositorio.dmiIpcdvCfi = dmiIpcdvCfi;

duTilHatfiLdyL = @(yL,zL)([(yL*zL)/fuL(yL,zL)^3; zL^2/fuL(yL,zL)^3]);
dAlcdyL = @(yL,zL)([-(3*R2*yL*zL^2)/fuL(yL,zL)^5, -(R2*zL*(fuL(yL,zL)^2 - 3*yL^2))/fuL(yL,zL)^5; 
    -(R2*zL*(fuL(yL,zL)^2 - 3*yL^2))/fuL(yL,zL)^5, (R2*yL*(2*fuL(yL,zL)^2 - 3*yL^2))/fuL(yL,zL)^5]);
dAvCfidyL = @(yL,zL)([R2, zeros(1,2); zeros(2,1), dAlcdyL(yL,zL)]);
dvCfidyL = @(yL,zL,teta2p,yLp,zLp)(([0, duTilHatfiLdyL(yL,zL).']*AvCfi(yL,zL) +...
    [0, fuTilHatfiL(yL,zL).']*dAvCfidyL(yL,zL))*...
    [teta2p; yLp; zLp]);
repositorio.dvCfidyL = dvCfidyL;
duTilHatfiLdzL = @(yL,zL)([-yL^2/fuL(yL,zL)^3; -(yL*zL)/fuL(yL,zL)^3]);
dAlcdzL = @(yL,zL)([(R2*zL*(2*fuL(yL,zL)^2 - 3*zL^2))/fuL(yL,zL)^5, -(R2*yL*(fuL(yL,zL)^2 - 3*zL^2))/fuL(yL,zL)^5;
    -(R2*yL*(fuL(yL,zL)^2 - 3*zL^2))/fuL(yL,zL)^5, -(3*R2*yL^2*zL)/fuL(yL,zL)^5]);
dAvCfidzL = @(yL,zL)([R2, zeros(1,2); zeros(2,1), dAlcdzL(yL,zL)]);
dvCfidzL = @(yL,zL,teta2p,yLp,zLp)(([0, duTilHatfiLdzL(yL,zL).']*AvCfi(yL,zL) +...
    [0, fuTilHatfiL(yL,zL).']*dAvCfidzL(yL,zL))*...
    [teta2p; yLp; zLp]);
repositorio.dvCfidzL = dvCfidzL;
dvCfiduT = @(uTilAc,uTilpAc)(dvCfidyL(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc))*...
    dyLduT +...
    dvCfidzL(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc))*...
    dzLduT);

fdT2cdyL = @(yL,zL,teta2p,yLp,zLp)(...
    (dmiIpcdvCfi(fvCfi(yL,zL,teta2p,yLp,zLp))*dvCfidyL(yL,zL,teta2p,yLp,zLp)*...
    fFn(fuL(yL,zL)) +...
    miIpc(fvCfi(yL,zL,teta2p,yLp,zLp))*...
    fdFnduL(fuL(yL,zL))*duLdyL(yL,zL))*...
    (R2 + fg));
repositorio.fdT2cRubdyL = fdT2cdyL;
fdT2cdzL = @(yL,zL,teta2p,yLp,zLp)(...
    (dmiIpcdvCfi(fvCfi(yL,zL,teta2p,yLp,zLp))*dvCfidzL(yL,zL,teta2p,yLp,zLp)*...
    fFn(fuL(yL,zL)) +...
    miIpc(fvCfi(yL,zL,teta2p,yLp,zLp))*...
    fdFnduL(fuL(yL,zL))*duLdzL(yL,zL))*...
    (R2 + fg));
repositorio.fdT2cRubdzL = fdT2cdzL;
fdvCfidteta2p = @(yL,zL)([1, fuTilHatfiL(yL,zL).']*AvCfi(yL,zL)*[1; 0; 0]);
fdT2cdteta2p = @(yL,zL,teta2p,yLp,zLp)(...
    dmiIpcdvCfi(fvCfi(yL,zL,teta2p,yLp,zLp))*fdvCfidteta2p(yL,zL)*fFn(fuL(yL,zL))*...
    (R2 + fg));
repositorio.fdT2cRubdteta2p = fdT2cdteta2p;
fdvCfidyLp = @(yL,zL)([1, fuTilHatfiL(yL,zL).']*AvCfi(yL,zL)*[0; 1; 0]);
fdT2cdyLp = @(yL,zL,teta2p,yLp,zLp)(...
    dmiIpcdvCfi(fvCfi(yL,zL,teta2p,yLp,zLp))*fdvCfidyLp(yL,zL)*fFn(fuL(yL,zL))*...
    (R2 + fg));
repositorio.fdT2cRubdyLp = fdT2cdyLp;
fdvCfidzLp = @(yL,zL)([1, fuTilHatfiL(yL,zL).']*AvCfi(yL,zL)*[0; 0; 1]);
fdT2cdzLp = @(yL,zL,teta2p,yLp,zLp)(...
    dmiIpcdvCfi(fvCfi(yL,zL,teta2p,yLp,zLp))*fdvCfidzLp(yL,zL)*fFn(fuL(yL,zL))*...
    (R2 + fg));
repositorio.fdT2cRubdzLp = fdT2cdzLp;

dmiIpcduT = @(uTilAc,uTilpAc)(...
    dmiIpcdvCfi(fvCfi(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc)))*...
    dvCfiduT(uTilAc,uTilpAc));
dT2cduTAc = @(uTilAc,uTilpAc)(...
    (miIpc(fvCfi(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc)))*...
    dFnduT(uTilAc) +...
    dmiIpcduT(uTilAc,uTilpAc)*fFn(fuL(fyL(uTilAc),fzL(uTilAc))))*...
    (R2 + fg));

dyLuLm1dyL = @(yL,zL)(fuL(yL,zL)^-1*(1 - (yL/fuL(yL,zL))^2));
dyLuLm1dzL = @(yL,zL)(-fuL(yL,zL)^-1*(yL*zL/fuL(yL,zL)^2));
dyLuLm1duT = @(uTilAc)(dyLuLm1dyL(fyL(uTilAc),fzL(uTilAc))*dyLduT + dyLuLm1dzL(fyL(uTilAc),fzL(uTilAc))*dzLduT);
dFyduTAc = @(uTilAc)(dyLuLm1duT(uTilAc)*fFn(fuL(fyL(uTilAc),fzL(uTilAc))) +...
    (fyL(uTilAc)/fuL(fyL(uTilAc),fzL(uTilAc)))*dFnduT(uTilAc));
dzLuLm1dyL = @(yL,zL)(-fuL(yL,zL)^-1*(zL*yL/fuL(yL,zL)^2));
dzLuLm1dzL = @(yL,zL)(fuL(yL,zL)^-1*(1 - (zL/fuL(yL,zL))^2));
dzLuLm1duT = @(uTilAc)(dzLuLm1dyL(fyL(uTilAc),fzL(uTilAc))*dyLduT + dzLuLm1dzL(fyL(uTilAc),fzL(uTilAc))*dzLduT);
dFzduTAc = @(uTilAc)(dzLuLm1duT(uTilAc)*fFn(fuL(fyL(uTilAc),fzL(uTilAc))) +...
    (fzL(uTilAc)/fuL(fyL(uTilAc),fzL(uTilAc)))*dFnduT(uTilAc));

dTrubcduTAc = @(uTilAc,uTilpAc)([dT2cduTAc(uTilAc,uTilpAc); 
    dFyduTAc(uTilAc); 
    dFzduTAc(uTilAc); 
    zeros(1,6); 
    zeros(1,6); 
    zeros(1,6)]);
dTrubbcdubTAc = @(uTilbAc,uTilbpAc)(dTrubcduTAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.dTrubbcdubTAc = dTrubbcdubTAc;

dvCfidupT = @(uTilAc)([[1, fuTilHatfiL(fyL(uTilAc),fzL(uTilAc)).']*AvCfi(fyL(uTilAc),fzL(uTilAc)), zeros(1,3)]);
dmiIpcdupT = @(uTilAc,uTilpAc)(...
    dmiIpcdvCfi(fvCfi(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc)))*...
    dvCfidupT(uTilAc));
dT2cdupTAc = @(uTilAc,uTilpAc)(dmiIpcdupT(uTilAc,uTilpAc)*fFn(fuL(fyL(uTilAc),fzL(uTilAc)))*...
    (R2 + fg));
dTrubcdupTAc = @(uTilAc,uTilpAc)([dT2cdupTAc(uTilAc,uTilpAc); 
    zeros(1,6); 
    zeros(1,6); 
    zeros(1,6); 
    zeros(1,6); 
    zeros(1,6)]);
dTrubbcdubpTAc = @(uTilbAc,uTilbpAc)(dTrubcdupTAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.dTrubbcdubpTAc = dTrubbcdubpTAc;

%Modelo de atrito descontÃ­nuo na origem:
switch opcao_sign
    case 1
        dmiIpdcdvCfi = @(vCfi)(miIp1dc*...
            (2*beta1miIpdc*exp(beta1miIpdc*abs(vCfi))/(1 + exp(beta1miIpdc*abs(vCfi)))^2) +...
            miIp2dc*...
            (2*beta2miIpdc*exp(beta2miIpdc*abs(vCfi))/(1 + exp(beta2miIpdc*abs(vCfi)))^2));
    case 2
        dmiIpdcdvCfi = @(vCfi)(miIp1dc*...
            (2*beta1miIpdc*exp(beta1miIpdc*abs(vCfi))/(1 + exp(beta1miIpdc*abs(vCfi)))^2) +...
            miIp2dc*...
            (2*beta2miIpdc*exp(beta2miIpdc*abs(vCfi))/(1 + exp(beta2miIpdc*abs(vCfi)))^2) +...
            inc_dsignS*dsignSdvCfi(vCfi)*...
            (miIp0dc +...
            miIp1dc*(1 - 2*(1 + exp(beta1miIpdc*abs(vCfi)))^-1) +...
            miIp2dc*(1 - 2*(1 + exp(beta2miIpdc*abs(vCfi)))^-1)));
    otherwise
        disp('other value')
end
dmiIpdcduT = @(uTilAc,uTilpAc)(...
    dmiIpdcdvCfi(fvCfi(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc)))*...
    dvCfiduT(uTilAc,uTilpAc));
dT2dcduTAc = @(uTilAc,uTilpAc)(...
    (miIpdc(fvCfi(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc)))*...
    dFnduT(uTilAc) +...
    dmiIpdcduT(uTilAc,uTilpAc)*fFn(fuL(fyL(uTilAc),fzL(uTilAc))))*...
    (R2 + fg));
dTrubdcduTAc = @(uTilAc,uTilpAc)([dT2dcduTAc(uTilAc,uTilpAc); 
    dFyduTAc(uTilAc); 
    dFzduTAc(uTilAc); 
    zeros(1,6); 
    zeros(1,6); 
    zeros(1,6)]);
dTrubbdcdubTAc = @(uTilbAc,uTilbpAc)(dTrubdcduTAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.dTrubbdcdubTAc = dTrubbdcdubTAc;

dmiIpdcdupT = @(uTilAc,uTilpAc)(...
    dmiIpdcdvCfi(fvCfi(fyL(uTilAc),fzL(uTilAc),fteta2(uTilpAc),fyL(uTilpAc),fzL(uTilpAc)))*...
    dvCfidupT(uTilAc,uTilpAc));
dT2dcdupTAc = @(uTilAc,uTilpAc)(dmiIpdcdupT(uTilAc,uTilpAc)*...
    fFn(fuL(fyL(uTilAc),fzL(uTilAc)))*...
    (R2 + fg));
dTrubdcdupTAc = @(uTilAc,uTilpAc)([dT2dcdupTAc(uTilAc,uTilpAc); 
    zeros(1,6); 
    zeros(1,6); 
    zeros(1,6); 
    zeros(1,6); 
    zeros(1,6)]);
dTrubbdcdubpTAc = @(uTilbAc,uTilbpAc)(dTrubdcdupTAc(uTilbAc + uTilAc_Eq,uTilbpAc));
repositorio.dTrubbdcdubpTAc = dTrubbdcdubpTAc;
%}

%Matriz de ganho:
Bac = @(uTilAc)([K10, 0;
    0, (2*fyL(uTilAc)/L1)*k1/betaS;
    0, (2*fzL(uTilAc)/L1)*k1/betaS;
    0, k1/betaS;
    0, 0;
    0, 0]);
BbAc = @(uTilbAc)(Bac(uTilbAc + uTilAc_Eq));
repositorio.BbAc = BbAc;

dBb1duTAc = zeros(6,6);
dBb2duTAc = [ 0, 0, 0, 0, 0, 0;
    0, (2*k1)/(L1*betaS), 0, 0, 0, 0;
    0, 0, (2*k1)/(L1*betaS), 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0];
%Obs.: dBbUcduTAc = dBb1duTAc*u1b + dBb2duTAc*tetaP + Bb*dUcduTAc
argumentos.dBb1duTAc = dBb1duTAc;
argumentos.dBb2duTAc = dBb2duTAc;

%Termos nã0 lineares concentrados numa única função:
TrubbcAc = repositorio.TrubbcAc;
fTilbc = @(uTilbAc,uTilbpAc,uTilbppAc)(GubupbAc(uTilbAc,uTilbpAc) +...
    FbuAc(uTilbAc) + TcAc(uTilbAc,uTilbpAc) + fM2bAc(uTilbAc)*uTilbppAc +...
    TrubbcAc(uTilbAc,uTilbpAc));
argumentos.fTilbc = fTilbc;
TrubbdcAc = repositorio.TrubbdcAc;
fTilbdc = @(uTilbAc,uTilbpAc,uTilbppAc)(GubupbAc(uTilbAc,uTilbpAc) +...
    FbuAc(uTilbAc) + TdcAc(uTilbAc,uTilbpAc) + fM2bAc(uTilbAc)*uTilbppAc +...
    TrubbdcAc(uTilbAc,uTilbpAc));
argumentos.fTilbdc = fTilbdc;

dTrubbcdubTAc = repositorio.dTrubbcdubTAc;
dfTilbcdubTAc = @(uTilbAc,uTilbpAc,uTilbppAc)(dGubupbdubTAc(uTilbAc,uTilbpAc) +...
    dFbudubTAc(uTilbAc) + dTcdubTAc(uTilbAc,uTilbpAc) +...
    dM2bubppdubT(uTilbAc,uTilbppAc) + dTrubbcdubTAc(uTilbAc,uTilbpAc));
argumentos.dfTilbcdubTAc = dfTilbcdubTAc;
dTrubbdcdubTAc = repositorio.dTrubbdcdubTAc;
dfTilbdcdubTAc = @(uTilbAc,uTilbpAc,uTilbppAc)(dGubupbdubTAc(uTilbAc,uTilbpAc) +...
    dFbudubTAc(uTilbAc) + dTdcdubTAc(uTilbAc,uTilbpAc) +...
    dM2bubppdubT(uTilbAc,uTilbppAc) + dTrubbdcdubTAc(uTilbAc,uTilbpAc));
argumentos.dfTilbdcdubTAc = dfTilbdcdubTAc;
dTrubbcdubpTAc = repositorio.dTrubbcdubpTAc;
dfTilbcdubpTAc = @(uTilbAc,uTilbpAc)(dGubupbdubpTAc(uTilbAc,uTilbpAc) +...
    dTcdubpTAc(uTilbAc,uTilbpAc) + dTrubbcdubpTAc(uTilbAc,uTilbpAc));
argumentos.dfTilbcdubpTAc = dfTilbcdubpTAc;
dTrubbdcdubpTAc = repositorio.dTrubbdcdubpTAc;
dfTilbdcdubpTAc = @(uTilbAc,uTilbpAc)(dGubupbdubpTAc(uTilbAc,uTilbpAc) +...
    dTdcdubpTAc(uTilbAc,uTilbpAc) + dTrubbdcdubpTAc(uTilbAc,uTilbpAc));
argumentos.dfTilbdcdubpTAc = dfTilbdcdubpTAc;
dfTilbdubppTAc = @(uTilbAc)(fM2bAc(uTilbAc));
argumentos.dfTilbdubppTAc = dfTilbdubppTAc;

%%
%SIMULAÇÃO
Id2 = eye(2);
Id4 = eye(4);
switch simulacao
    case 1
        M = Mt;
        argumentos.M = M;
        K = Kt;
        argumentos.K = K;
        C = Ct;
        argumentos.C = C;
        B2 = Bt;
        argumentos.B2 = B2;
        Normalb0 = repositorio.Normalb0;
        Tatr3c0 = @(teta3p)(mi3c(teta3p)*Normalb0(teta3p)*R3);
        fteta3 = @(upTil)(Id2(2,:)*upTil);
        Tc0 = @(uTilbp)([0; Tatr3c0(fteta3(uTilbp))]);
        repositorio.Tc0 = Tc0;
        %Entradas:
        Tal = argumentos.Tal; %s
        wRef = argumentos.wRef;
        fteta1pUc = @(t)(wRef*tanh(pi*t/Tal));
        repositorio.fteta1pUc = fteta1pUc;
        fUc = @(t)(wRef*(pi/Tal)^-1*log(cosh(pi*t/Tal)));
        repositorio.fUc = fUc;
        fUpc = @(t)(wRef*tanh(pi*t/Tal));
        repositorio.fUpc = fUpc;
        %Derivadas e outras expressões (para integração):
        dNormalb0dw3 = repositorio.dNormalb0dw3;
        dw3duTilpT = Id2(2,:);
        dTatr3c0dw3 = @(teta3p)((dmi3cdw3(teta3p)*Normalb0(teta3p) + mi3c(teta3p)*dNormalb0dw3(teta3p))*R3);
        repositorio.dTatr3c0dw3 = dTatr3c0dw3;
        dTatr3c0duTilpT = @(uTilbp)(dTatr3c0dw3(fteta3(uTilbp))*dw3duTilpT);
        dTc0duTilpT = @(uTilbp)([zeros(1,2); dTatr3c0duTilpT(uTilbp)]);
        repositorio.dTc0duTilpT = dTc0duTilpT;
        %Formas vetorizadas:
        Normalb0vt = repositorio.Normalb0vt;
        Tatr3c0vt = @(teta3p)(mi3cVt(teta3p).*Normalb0vt(teta3p)*R3);
        repositorio.Tatr3c0vt = Tatr3c0vt;
    case 2
        M = [Mt, zeros(2); zeros(2), Mu];
        argumentos.M = M;
        K = [Kt, zeros(2); zeros(2), Ku];
        argumentos.K = K;
        C = [Ct, zeros(2); zeros(2), Cu];
        argumentos.C = C;
        B4 = [Bt, zeros(2,1); zeros(2,1), Bu];
        argumentos.B4 = B4;
        %Vetor posição de equilíbrio:
        uTilEqMod1 = [0; 0; u2_Eq; u3_Eq];
        argumentos.uTilEqMod1 = uTilEqMod1;
        %Funções:
        fteta3 = @(upTil)(Id4(2,:)*upTil);
        fu3 = @(uTil)(Id4(4,:)*uTil);
        Tc = @(uTilb,uTilbp)([0; 
            Tatr3c(fteta3(uTilbp),fu3(uTilb)); 
            0; 
            Normalb(fu3(uTilb))]);
        repositorio.Tc = Tc;
        Tdc = @(uTilb,uTilbp)([0;
            Tatr3dc(fteta3(uTilbp),fu3(uTilb));
            0;
            Normalb(fu3(uTilb))]);
        repositorio.Tdc = Tdc;
        %Entradas:
        uRef = u11Ref;
        %uRef = u12Ref;
        Tal = argumentos.Tal; %s
        wRef = argumentos.wRef;
        fteta1pUc = @(t)(wRef*tanh(pi*t/Tal));
        repositorio.fteta1pUc = fteta1pUc;
        fteta1Uc = @(t)(wRef*(pi/Tal)^-1*log(cosh(pi*t/Tal)));
        repositorio.fteta1Uc = fteta1Uc;
        fU = argumentos.fU;
        fu1Uc = @(t)(fU*uRef*tanh(pi*t/Tal));
        ftetaPUc = @(t)(betaS*fu1Uc(t));
        fUc = @(t)([fteta1Uc(t); ftetaPUc(t)]);
        repositorio.fUc = fUc;
        %Derivadas e outras expressões (para integração):
        dNormalbduTilbT = @(uTilb)([0, 0, 0, dNormalbdu3b(fu3(uTilb))]);
        dTatr3cduTilbT = @(uTilb,uTilbp)(mi3c(fteta3(uTilbp))*dNormalbduTilbT(uTilb)*R3);
        dTcduTilbT = @(uTilb,uTilbp)([zeros(1,4); 
            dTatr3cduTilbT(uTilb,uTilbp); 
            zeros(1,4); 
            dNormalbduTilbT(uTilb)]);
        repositorio.dTcduTilbT = dTcduTilbT;
        dmi3cduTilbpT = @(uTilbp)([0, dmi3cdw3(fteta3(uTilbp)), 0, 0]);
        dTatr3cduTilbpT = @(uTilb,uTilbp)(dmi3cduTilbpT(uTilbp)*Normalb(fu3(uTilb))*R3);
        dTcduTilbpT = @(uTilb,uTilbp)([zeros(1,4); 
            dTatr3cduTilbpT(uTilb,uTilbp); 
            zeros(1,4); 
            zeros(1,4)]);
        repositorio.dTcduTilbpT = dTcduTilbpT;
    case 3
        M1 = M1ac;
        argumentos.M1 = M1;
        K = Kac;
        argumentos.K = K;
        C = Cac;
        argumentos.C = C;
        B6 = @(uTilb)(BbAc(uTilb));
        repositorio.B6 = B6;
        fTilc = @(uTilb,uTilbp,uTilbpp)(fTilbc(uTilb,uTilbp,uTilbpp));
        repositorio.fTilc = fTilc;
        fTildc = @(uTilb,uTilbp,uTilbpp)(fTilbdc(uTilb,uTilbp,uTilbpp));
        repositorio.fTildc = fTildc;
        %Entradas:
        uRef = u11Ref;
        %uRef = u12Ref;
        Tal = argumentos.Tal; %s
        wRef = argumentos.wRef;
        fteta1pUc = @(t)(wRef*tanh(pi*t/Tal));
        repositorio.fteta1pUc = fteta1pUc;
        fteta1Uc = @(t)(wRef*(pi/Tal)^-1*log(cosh(pi*t/Tal)));
        repositorio.fteta1Uc = fteta1Uc;
        fU = argumentos.fU;
        fu1Uc = @(t)(fU*uRef*tanh(pi*t/Tal));
        ftetaPUc = @(t)(betaS*fu1Uc(t));
        fUc = @(t)([fteta1Uc(t); ftetaPUc(t)]);
        repositorio.fUc = fUc;
        %Derivadas e outras expressções (para integração):
        dfTilcdubT = @(uTilb,uTilbp,uTilbpp)(dfTilbcdubTAc(uTilb,uTilbp,uTilbpp));
        repositorio.dfTilcdubT = dfTilcdubT;
        dfTilcdubpT = @(uTilb,uTilbp)(dfTilbcdubpTAc(uTilb,uTilbp));
        repositorio.dfTilcdubpT = dfTilcdubpT;
        dfTildubppT = @(uTilb)(dfTilbdubppTAc(uTilb));
        repositorio.dfTildubppT = dfTildubppT;
        dfTildcdubT = @(uTilb,uTilbp,uTilbpp)(dfTilbdcdubTAc(uTilb,uTilbp,uTilbpp));
        repositorio.dfTildcdubT = dfTildcdubT;
        dfTildcdubpT = @(uTilb,uTilbp)(dfTilbdcdubpTAc(uTilb,uTilbp));
        repositorio.dfTildcdubpT = dfTildcdubpT;
        argumentos.dB61duT = dBb1duTAc;
        argumentos.dB62duT = dBb2duTAc;
    case 4
        %Obs.:
        %(1) Simulação da trajetória desejada e da camada limite dinâmica
        %(2) Não há necessidade de se definir elementos da planta
    case 5
        M = Mt;
        argumentos.M = M;
        K = Kt;
        argumentos.K = K;
        C = Ct;
        argumentos.C = C;
        B2 = Bt;
        argumentos.B2 = B2;
        Normalb0 = repositorio.Normalb0;
        Tatr3c0 = @(teta3p)(mi3c(teta3p)*Normalb0(teta3p)*R3);
        fteta3 = @(upTil)(Id2(2,:)*upTil);
        Tc0 = @(uTilbp)([0; Tatr3c0(fteta3(uTilbp))]);
        repositorio.Tc0 = Tc0;
        
        %Derivadas e outras expressões (para integração):
        dNormalb0dw3 = repositorio.dNormalb0dw3;
        dw3duTilpT = Id2(2,:);
        dTatr3c0dw3 = @(teta3p)((dmi3cdw3(teta3p)*Normalb0(teta3p) + mi3c(teta3p)*dNormalb0dw3(teta3p))*R3);
        repositorio.dTatr3c0dw3 = dTatr3c0dw3;
        dTatr3c0duTilpT = @(uTilbp)(dTatr3c0dw3(fteta3(uTilbp))*dw3duTilpT);
        dTc0duTilpT = @(uTilbp)([zeros(1,2); dTatr3c0duTilpT(uTilbp)]);
        repositorio.dTc0duTilpT = dTc0duTilpT;
        
        %Formas vetorizadas:
        Normalb0vt = repositorio.Normalb0vt;
        Tatr3c0vt = @(teta3p)(mi3cVt(teta3p).*Normalb0vt(teta3p)*R3);
        repositorio.Tatr3c0vt = Tatr3c0vt;
    case 6
        M = [Mt, zeros(2); zeros(2), Mu];
        argumentos.M = M;
        K = [Kt, zeros(2); zeros(2), Ku];
        argumentos.K = K;
        C = [Ct, zeros(2); zeros(2), Cu];
        argumentos.C = C;
        B4 = [Bt, zeros(2,1); zeros(2,1), Bu];
        argumentos.B4 = B4;
        fteta3 = @(upTil)(Id4(2,:)*upTil);
        fu3 = @(uTil)(Id4(4,:)*uTil);
        Tc = @(uTilb,uTilbp)([0; 
            Tatr3c(fteta3(uTilbp),fu3(uTilb)); 
            0; 
            Normalb(fu3(uTilb))]);
        repositorio.Tc = Tc;
        Tdc = @(uTilb,uTilbp)([0;
            Tatr3dc(fteta3(uTilbp),fu3(uTilb));
            0;
            Normalb(fu3(uTilb))]);
        repositorio.Tdc = Tdc;
        
        %Função tetaP:
        Auc = Id4(:,1:2);
        argumentos.Auc = Auc;
        fuTilc = @(uTilb)(Auc.'*uTilb);
        repositorio.fuTilc = fuTilc;
        tetaPuc = @(uTilcp)(tetaP(fw3(uTilcp)));
        repositorio.tetaPuc = tetaPuc;
        tetaPub = @(uTilbp)(tetaPuc(fuTilc(uTilbp)));
        repositorio.tetaPub = tetaPub;
        
        %Derivadas e outras expressÃµes (para integraÃ§Ã£o):
        dNormalbduTilbT = @(uTilb)([0, 0, 0, dNormalbdu3b(fu3(uTilb))]);
        dTatr3cduTilbT = @(uTilb,uTilbp)(mi3c(fteta3(uTilbp))*dNormalbduTilbT(uTilb)*R3);
        dTcduTilbT = @(uTilb,uTilbp)([zeros(1,4); 
            dTatr3cduTilbT(uTilb,uTilbp); 
            zeros(1,4); 
            dNormalbduTilbT(uTilb)]);
        repositorio.dTcduTilbT = dTcduTilbT;
        dmi3cduTilbpT = @(uTilbp)([0, dmi3cdw3(fteta3(uTilbp)), 0, 0]);
        dTatr3cduTilbpT = @(uTilb,uTilbp)(dmi3cduTilbpT(uTilbp)*Normalb(fu3(uTilb))*R3);
        dTcduTilbpT = @(uTilb,uTilbp)([zeros(1,4); 
            dTatr3cduTilbpT(uTilb,uTilbp); 
            zeros(1,4); 
            zeros(1,4)]);
        repositorio.dTcduTilbpT = dTcduTilbpT;
        
        %Função dtetaPubduTilbT:
        dw3duTilcpT = Id2(2,:);
        duTilcpduTilbpT = Auc.';
        fteta3 = @(upTil)(Id2(2,:)*upTil);
        dtetaPucduTilcpT = @(uTilcp)(dtetaPdw3(fteta3(uTilcp))*dw3duTilcpT);
        dtetaPubduTilbpT = @(uTilbp)(dtetaPucduTilcpT(fuTilc(uTilbp))*duTilcpduTilbpT);
        repositorio.dtetaPubduTilbpT = dtetaPubduTilbpT;
    case 7
        M1 = M1ac;
        argumentos.M1 = M1;
        K = Kac;
        argumentos.K = K;
        C = Cac;
        argumentos.C = C;
        B6 = @(uTilb)(BbAc(uTilb));
        repositorio.B6 = B6;
        fteta3 = @(upTil)(Id2(2,:)*upTil);
        repositorio.fteta3 = fteta3;
        fTilc = @(uTilb,uTilbp,uTilbpp)(fTilbc(uTilb,uTilbp,uTilbpp));
        repositorio.fTilc = fTilc;
        fTildc = @(uTilb,uTilbp,uTilbpp)(fTilbdc(uTilb,uTilbp,uTilbpp));
        repositorio.fTildc = fTildc;
        %Derivadas e outras expressões (para integração):
        dfTilbcdubT = @(uTilb,uTilbp,uTilbpp)(dfTilbcdubTAc(uTilb,uTilbp,uTilbpp));
        repositorio.dfTilcdubT = dfTilbcdubT;
        dfTilbcdubpT = @(uTilb,uTilbp)(dfTilbcdubpTAc(uTilb,uTilbp));
        repositorio.dfTilcdubpT = dfTilbcdubpT;
        dfTilbdubppT = @(uTilb)(dfTilbdubppTAc(uTilb));
        repositorio.dfTildubppT = dfTilbdubppT;
        dfTilbdcdubT = @(uTilb,uTilbp,uTilbpp)(dfTilbdcdubTAc(uTilb,uTilbp,uTilbpp));
        repositorio.dfTildcdubT = dfTilbdcdubT;
        dfTilbdcdubpT = @(uTilb,uTilbp)(dfTilbdcdubpTAc(uTilb,uTilbp));
        repositorio.dfTildcdubpT = dfTilbdcdubpT;
        dfTilbdcdubpT = @(uTilb,uTilbp)(dfTilbdcdubppTAc(uTilb,uTilbp));
        repositorio.dfTildcdubpT = dfTilbdcdubpT;
        
        %Função tetaP:
        Auc = [Id6(:,1), Id6(:,6)];
        argumentos.Auc = Auc;
        fuTilc = @(uTilb)(Auc.'*uTilb);
        repositorio.fuTilc = fuTilc;
        tetaPuc = @(uTilcp)(tetaP(fw3(uTilcp)));
        repositorio.tetaPuc = tetaPuc;
        tetaPub = @(uTilbp)(tetaPuc(fuTilc(uTilbp)));
        repositorio.tetaPub = tetaPub;
        
        %Função dtetaPubduTilbT:
        dw3duTilcpT = Id2(2,:);
        duTilcpduTilbpT = Auc.';
        dtetaPucduTilcpT = @(uTilcp)(dtetaPdw3(fteta3(uTilcp))*dw3duTilcpT);
        dtetaPubduTilbpT = @(uTilbp)(dtetaPucduTilcpT(fuTilc(uTilbp))*duTilcpduTilbpT);
        repositorio.dtetaPubduTilbpT = dtetaPubduTilbpT;
        
        dB61duTilbT = dBb1duTAc;
        dB62duTilbT = dBb2duTAc;
        dZduTilcT = Au;
        dZduTilcpT = Av;
        duTilcduTilbT = Auc.';
        dB6UccTetaPduTilbT = @(t,uTilb,uTilbp,fiBd,mid,Zd0)(...
            Uccu(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Zd0)*dB61duTilbT +...
            tetaPub(uTilbp)*dB62duTilbT +...
            B6(uTilb)*[DUccuDZT(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Zd0)*dZduTilcT*duTilcduTilbT; zeros(size(uTilb.'))]);
        repositorio.dB6UccTetaPduTilbT = dB6UccTetaPduTilbT;
        dB6UccTetaPduTilbpT = @(t,uTilb,uTilbp,fiBd,mid,Zd0)(...
            B6(uTilb)*[DUccuDZT(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Zd0)*dZduTilcpT*duTilcpduTilbpT; dtetaPubduTilbpT(uTilbp)]);
        repositorio.dB6UccTetaPduTilbpT = dB6UccTetaPduTilbpT;
        
        dB6UcdcTetaPduTilbT = @(t,uTilb,uTilbp,fiBd,mid,Zd0)(...
            Ucdcu(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Zd0)*dB61duTilbT +...
            tetaPub(uTilbp)*dB62duTilbT +...
            B6(uTilb)*[DUcdcuDZT(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Z0)*dZduTilcT*duTilcduTilbT; zeros(size(uTilb.'))]);
        repositorio.dB6UcdcTetaPduTilbT = dB6UcdcTetaPduTilbT;
        dB6UcdcTetaPduTilbpT = @(t,uTilb,uTilbp,fiBd,mid,Zd0)(...
            B6(uTilb)*[DUcdcuDZT(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Z0)*dZduTilcpT*duTilcpduTilbpT; dtetaPubduTilbpT(uTilbp.')]);
        repositorio.dB6UcdcTetaPduTilbpT = dB6UcdcTetaPduTilbpT;
    case 8
        %Entradas:
        uRef = u11Ref;
        %uRef = u12Ref;
        Tal = argumentos.Tal; %s
        %Funções adimensionais:
        fi = @(nu)(tanh(pi*nu));
        IntFi = @(nu)(pi^-1*log(cosh(pi*nu)));
        dfidnu = @(nu)(pi*(ones(size(nu)) - tanh(pi*nu).^2));
        d2fidnu2 = @(nu)(2*pi^2*tanh(pi*nu).*(tanh(pi*nu).^2 - ones(size(nu))));
        %Funções auxiliares:
        wRef = argumentos.wRef;
        fteta1pUc = @(t)(wRef*fi(t/Tal));
        repositorio.fteta1pUc = fteta1pUc;
        fteta1Uc = @(t)(wRef*Tal*IntFi(t/Tal));
        repositorio.fteta1Uc = fteta1Uc;
        fteta1ppUc = @(t)((wRef/Tal)*dfidnu(t/Tal));
        fU = argumentos.fU;
        fu1Uc = @(t)(fU*uRef*fi(t/Tal));
        ftetaPUc = @(t)(betaS*fu1Uc(t));
        fu1pUc = @(t)(fU*(uRef/Tal)*dfidnu(t/Tal));
        ftetaPpUc = @(t)(betaS*fu1pUc(t));
        fu1ppUc = @(t)(fU*(uRef/Tal^2)*d2fidnu2(t/Tal));
        ftetaPppUc = @(t)(betaS*fu1ppUc(t));
        %Entradas e suas derivadas temporais:
        fUc = @(t)([fteta1Uc(t); ftetaPUc(t)]);
        repositorio.fUc = fUc;
        fUpc = @(t)([fteta1pUc(t); ftetaPpUc(t)]);
        repositorio.fUpc = fUpc;
        fUppc = @(t)([fteta1ppUc(t); ftetaPppUc(t)]);
        repositorio.fUppc = fUppc;
    case 9
        %Entradas:
        uRef = u11Ref;
        %uRef = u12Ref;
        Tal = argumentos.Tal; %s
        %Funções adimensionais:
        fi = @(nu)(tanh(pi*nu));
        IntFi = @(nu)(pi^-1*log(cosh(pi*nu)));
        dfidnu = @(nu)(pi*(ones(size(nu)) - tanh(pi*nu).^2));
        d2fidnu2 = @(nu)(2*pi^2*tanh(pi*nu).*(tanh(pi*nu).^2 - ones(size(nu))));
        %Funções auxiliares:
        wRef = argumentos.wRef;
        fteta1pUc = @(t)(wRef*fi(t/Tal));
        repositorio.fteta1pUc = fteta1pUc;
        fteta1Uc = @(t)(wRef*Tal*IntFi(t/Tal));
        repositorio.fteta1Uc = fteta1Uc;
        fteta1ppUc = @(t)((wRef/Tal)*dfidnu(t/Tal));
        fU = argumentos.fU;
        fu1Uc = @(t)(fU*uRef*fi(t/Tal));
        ftetaPUc = @(t)(betaS*fu1Uc(t));
        fu1pUc = @(t)(fU*(uRef/Tal)*dfidnu(t/Tal));
        ftetaPpUc = @(t)(betaS*fu1pUc(t));
        fu1ppUc = @(t)(fU*(uRef/Tal^2)*d2fidnu2(t/Tal));
        ftetaPppUc = @(t)(betaS*fu1ppUc(t));
        %Entradas e suas derivadas temporais:
        fUc = @(t)([fteta1Uc(t); ftetaPUc(t)]);
        repositorio.fUc = fUc;
        fUpc = @(t)([fteta1pUc(t); ftetaPpUc(t)]);
        repositorio.fUpc = fUpc;
        fUppc = @(t)([fteta1ppUc(t); ftetaPppUc(t)]);
        repositorio.fUppc = fUppc;
    case 10
        M1 = M1ac;
        argumentos.M1 = M1;
        K = Kac;
        argumentos.K = K;
        C = Cac;
        argumentos.C = C;
        B6 = @(uTilb)(BbAc(uTilb));
        repositorio.B6 = B6;
        fteta3 = @(upTil)(Id2(2,:)*upTil);
        repositorio.fteta3 = fteta3;
        fTilc = @(uTilb,uTilbp,uTilbpp)(fTilbc(uTilb,uTilbp,uTilbpp));
        repositorio.fTilc = fTilc;
        fTildc = @(uTilb,uTilbp,uTilbpp)(fTilbdc(uTilb,uTilbp,uTilbpp));
        repositorio.fTildc = fTildc;
        %Derivadas e outras expressÃµes (para integração):
        dfTilbcdubT = @(uTilb,uTilbp,uTilbpp)(dfTilbcdubTAc(uTilb,uTilbp,uTilbpp));
        repositorio.dfTilcdubT = dfTilbcdubT;
        dfTilbcdubpT = @(uTilb,uTilbp)(dfTilbcdubpTAc(uTilb,uTilbp));
        repositorio.dfTilcdubpT = dfTilbcdubpT;
        dfTilbdubppT = @(uTilb)(dfTilbdubppTAc(uTilb));
        repositorio.dfTildubppT = dfTilbdubppT;
        dfTilbdcdubT = @(uTilb,uTilbp,uTilbpp)(dfTilbdcdubTAc(uTilb,uTilbp,uTilbpp));
        repositorio.dfTildcdubT = dfTilbdcdubT;
        dfTilbdcdubpT = @(uTilb,uTilbp)(dfTilbdcdubpTAc(uTilb,uTilbp));
        repositorio.dfTildcdubpT = dfTilbdcdubpT;
        dfTilbdcdubpT = @(uTilb,uTilbp)(dfTilbdcdubppTAc(uTilb,uTilbp));
        repositorio.dfTildcdubpT = dfTilbdcdubpT;
        
        %Função tetaP:
        Auc = [Id6(:,1), Id6(:,6)];
        argumentos.Auc = Auc;
        fuTilc = @(uTilb)(Auc.'*uTilb);
        repositorio.fuTilc = fuTilc;
        tetaPuc = @(uTilcp)(tetaP(fw3(uTilcp)));
        repositorio.tetaPuc = tetaPuc;
        tetaPub = @(uTilbp)(tetaPuc(fuTilc(uTilbp)));
        repositorio.tetaPub = tetaPub;
        
        %Função dtetaPubduTilbT:
        dw3duTilcpT = Id2(2,:);
        duTilcpduTilbpT = Auc.';
        dtetaPucduTilcpT = @(uTilcp)(dtetaPdw3(fteta3(uTilcp))*dw3duTilcpT);
        dtetaPubduTilbpT = @(uTilbp)(dtetaPucduTilcpT(fuTilc(uTilbp))*duTilcpduTilbpT);
        repositorio.dtetaPubduTilbpT = dtetaPubduTilbpT;
        
        dB61duTilbT = dBb1duTAc;
        dB62duTilbT = dBb2duTAc;
        dZduTilcT = Au;
        dZduTilcpT = Av;
        duTilcduTilbT = Auc.';
        dB6UccTetaPduTilbT = @(t,uTilb,uTilbp,fiBd,mid,Zd0)(...
            Uccu(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Zd0)*dB61duTilbT +...
            tetaPub(uTilbp)*dB62duTilbT +...
            B6(uTilb)*[DUccuDZT(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Zd0)*dZduTilcT*duTilcduTilbT; zeros(size(uTilb.'))]);
        repositorio.dB6UccTetaPduTilbT = dB6UccTetaPduTilbT;
        dB6UccTetaPduTilbpT = @(t,uTilb,uTilbp,fiBd,mid,Zd0)(...
            B6(uTilb)*[DUccuDZT(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Zd0)*dZduTilcpT*duTilcpduTilbpT; dtetaPubduTilbpT(uTilbp)]);
        repositorio.dB6UccTetaPduTilbpT = dB6UccTetaPduTilbpT;
        
        dB6UcdcTetaPduTilbT = @(t,uTilb,uTilbp,fiBd,mid,Zd0)(...
            Ucdcu(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Zd0)*dB61duTilbT +...
            tetaPub(uTilbp)*dB62duTilbT +...
            B6(uTilb)*[DUcdcuDZT(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Z0)*dZduTilcT*duTilcduTilbT; zeros(size(uTilb.'))]);
        repositorio.dB6UcdcTetaPduTilbT = dB6UcdcTetaPduTilbT;
        dB6UcdcTetaPduTilbpT = @(t,uTilb,uTilbp,fiBd,mid,Zd0)(...
            B6(uTilb)*[DUcdcuDZT(t,fuTilc(uTilb),fuTilc(uTilbp),fiBd,mid,Z0)*dZduTilcpT*duTilcpduTilbpT; dtetaPubduTilbpT(uTilbp.')]);
        repositorio.dB6UcdcTetaPduTilbpT = dB6UcdcTetaPduTilbpT;
    case 11
    case 12
    case 13
        %Obs.: (Ref.: simulação 4)
        %(1) Simulação apenas da trajetória desejada
        %(2) Não há necessidade de se definir elementos da planta
    otherwise
        disp('other value')
end
end