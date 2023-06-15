function [ argumentos ] = zzpreprocess2( simulacao, opcao_m2, opcao_sign,...
    inc_dsignS, opcao_atr, ind_LVC )
%UNTITLED32 Summary of this function goes here
%   simulacao:
%   opcao (Rotor 2, excentricidade): (1) Aço; (2) Alumínio; (3) Vazio
argumentos.simulacao = simulacao;
argumentos.opcao_m2 = opcao_m2;
argumentos.opcao_sign = opcao_sign;
argumentos.inc_dsignS = inc_dsignS;
argumentos.opcao_atr = opcao_atr;
argumentos.ind_LVC = ind_LVC;
%(1) ad=0: ausência de parâmetros concentrados
%(2) ad=1: inclusão plena dos parâmetros concentrados
ad = 1;
%Número de nós para a quadratura de Gauss-Legendre
nGL = 5;
argumentos.nGL = nGL;

percent_geral = 0.1;
argumentos.percent_geral = percent_geral;
%Espessura do eixo metálico:
%Obs.: Mudanças em d levam a mudanças na ordem dos modos de vibração: 
%(1) d = 3.175*10^-3: tX1, v, w, tX2, u1, u2
%(2) d = 1.500*10^-3: tX1, tX2, v, w, u1, u2
d_hat = 3.175*10^-3;
%d_hat = 1.5*10^-3;
dd = percent_geral*d_hat;
d = 3.175*10^-3; %m
%d = 1.5*10^-3; %m
argumentos.d = d;
rEixo = d/2;
argumentos.raio = rEixo;
A_hat = pi*d_hat^2/4;
argumentos.Area_hat = A_hat;
dAdd_hat = (pi*d_hat)/2;
argumentos.dAdd_hat = dAdd_hat;
A = pi*rEixo^2;
argumentos.Area = A;
%Momento polar de inércia do eixo:
Jp_hat = pi*d_hat^4/32; %m^4
argumentos.Jp_hat = Jp_hat;
dJpdd_hat = (pi*d_hat^3)/8;
argumentos.dJpdd_hat = dJpdd_hat;
Jp = pi*d^4/32; %m^4
argumentos.Jp = Jp;
%Momento da seção transversal do eixo:
J = Jp/2; %m^4
argumentos.J = J;
%Densidade do aço:
rho_hat = 7850; %kg/m^3
argumentos.rho_hat = rho_hat;
drho = percent_geral*rho_hat;
rho = 7850; %kg/m^3
argumentos.rho = rho;
%Shearing correct factor (tese Thiago Ritto, p.128 pdf):
ks = 6/7;
argumentos.ks = ks;
%Módulo de elasticidade (ou módulo de Young):
E_hat = 200*10^9; %Pa
argumentos.E_hat = E_hat;
dE = percent_geral*E_hat;
E = 200*10^9; %Pa
argumentos.E = E;
%Coeficiente de Poisson:
nu_hat = 0.3;
argumentos.nu_hat = nu_hat;
dnu = percent_geral*nu_hat;
nu = 0.3;
argumentos.nu = nu;
%Módulo de cisalhamento (Pa):
G_hat = E_hat/(2*(1 + nu_hat));
argumentos.G_hat = G_hat;
G = E/(2*(1 + nu));
argumentos.G = G;
dGdE_hat = 1/(2*nu_hat + 2);
argumentos.dGdE_hat = dGdE_hat;
dGdnu_hat = -(2*E_hat)/(2*nu_hat + 2)^2;
argumentos.dGdnu_hat = dGdnu_hat;
%Rotor 1:
L1b_hat = 0.026; %m
argumentos.L1b_hat = L1b_hat;
dL1b = percent_geral*L1b_hat;
L1b = 0.026; %m
Dr1_hat = ad*0.25; %m
dDr1 = percent_geral*Dr1_hat;
Dr1 = ad*0.25; %m
A1_hat = pi*Dr1_hat^2/4; %m^2
argumentos.A1_hat = A1_hat;
dAr1dDr1_hat = (pi*Dr1_hat)/2;
argumentos.dAr1dDr1_hat = dAr1dDr1_hat;
A1 = pi*Dr1^2/4; %m^2
Jp1_hat = pi*Dr1_hat^4/32; %m^4
argumentos.Jp1_hat = Jp1_hat;
dJp1dDr1_hat = (pi*Dr1_hat^3)/8;
argumentos.dJp1dDr1_hat = dJp1dDr1_hat;
Jp1 = pi*Dr1^4/32; %m^4
M1_hat = rho_hat*A1_hat*L1b_hat; %Massa do rotor 1
argumentos.M1_hat = M1_hat;
M1 = rho*A1*L1b; %Massa do rotor 1
argumentos.M1rt = M1;
%argumentos.M1 = M1;
Ip1_hat = rho_hat*Jp1_hat*L1b_hat; %Momento polar do rotor 1
argumentos.Ip1_hat = Ip1_hat;
alfaIp1 = 1;
Ip1 = alfaIp1*rho*Jp1*L1b; %Momento polar do rotor 1
argumentos.Ip1 = Ip1;

%Rotor 2:
alfaL2b = 1;
L2b_hat = alfaL2b*0.06; %m
argumentos.L2b_hat = L2b_hat;
dL2b = percent_geral*L2b_hat;
L2b = alfaL2b*0.06; %m
Dr2_hat = ad*0.14; %m
dDr2 = percent_geral*Dr2_hat;
Dr2 = ad*0.14; %m
%Dr2 = ad*0.02; %m
argumentos.Dr2 = Dr2;
R2 = Dr2/2;
argumentos.R2 = R2;
A2_hat = pi*Dr2_hat^2/4; %m^2
argumentos.A2_hat = A2_hat;
dAr2dDr2_hat = (pi*Dr2_hat)/2;
argumentos.dAr2dDr2_hat = dAr2dDr2_hat;
A2 = pi*Dr2^2/4; %m^2

switch opcao_m2
    case 1
        rhom2_hat = rho_hat;
        rhom2 = rho;
    case 2
        rhom2_hat = 2700; %kg/m3 (alumínio)
        rhom2 = 2700; %kg/m3 (alumínio)
    case 3
        rhom2_hat = 0; %furo vazio
        rhom2 = 0; %furo vazio
    otherwise
        disp('other value')
end
argumentos.rhom2_hat = rhom2_hat;
drhom2 = percent_geral*rhom2_hat;

nd = 0.37; %uLmax = 5.7mm; 120RPM; d = 3mm; rhom2 = 0; mic0 = 0 (AD, 17/07/2021)
Dm_hat = nd*((1/2)*Dr2_hat);
dDm = percent_geral*Dm_hat; 
Am2_hat = pi*Dm_hat^2/4;
argumentos.Am2_hat = Am2_hat;
dAm2dDm_hat = (pi*Dm_hat)/2;
argumentos.dAm2dDm_hat = dAm2dDm_hat;
m20_hat = rho_hat*Am2_hat*L2b_hat;
M2_hat = rho_hat*A2_hat*L2b_hat - m20_hat;
argumentos.M2_hat = M2_hat;
m2_hat = rhom2_hat*Am2_hat*L2b_hat;
argumentos.m2_hat = m2_hat;
M2b_hat = rho_hat*A2_hat*L2b_hat;
argumentos.M2b_hat = M2b_hat;
dm_hat = (rhom2_hat - rho_hat)*Am2_hat*L2b_hat;
argumentos.dm_hat = dm_hat;
n_l2 = 0.5; %0 < n_l2 < 1
l2_min_hat = Dm_hat/2;
l2_max_hat = (Dr2_hat - Dm_hat)/2;
l2_hat_ref = (1 - n_l2)*l2_min_hat + n_l2*l2_max_hat;
argumentos.l2_hat_ref = l2_hat_ref;
l2_hat = 0.035; %m
argumentos.l2_hat = l2_hat;
dl2 = percent_geral*l2_hat;
rCG_hat = l2_hat*(m2_hat - m20_hat)/(M2_hat + m2_hat);
argumentos.rCG_hat = rCG_hat;

Dm = nd*((1/2)*Dr2);
Am2 = pi*Dm^2/4;
m20 = rho*Am2*L2b;
M2 = rho*A2*L2b - m20;
argumentos.M2rt = M2;
m2 = rhom2*Am2*L2b;
argumentos.m2 = m2;
M2b = rho*A2*L2b;
argumentos.M2b = M2b;
dm = (rhom2 - rho)*Am2*L2b;
argumentos.dm = dm;
l2_min = Dm/2;
l2_max = (Dr2 - Dm)/2;
l2_ref = (1 - n_l2)*l2_min + n_l2*l2_max;
argumentos.l2_ref = l2_ref;
l2 = 0.035; %m
argumentos.l2 = l2;
rCG = l2*(m2 - m20)/(M2 + m2);
argumentos.rCG = rCG;
alfa_m2 = 30*(pi/180); %Ângulo de posição de m2 no plano YZ, 1o quadrante (rad)
argumentos.alfa_m2 = alfa_m2;

Jp_m2_hat = pi*Dm_hat^4/32 + l2_hat^2*Am2_hat;
Jp2_hat = pi*Dr2_hat^4/32 - Jp_m2_hat;
Jp2b_hat = pi*Dr2_hat^4/32;
argumentos.Jp2b_hat = Jp2b_hat;
dJp2bdDr2_hat = (pi*Dr2_hat^3)/8;
argumentos.dJp2bdDr2_hat = dJp2bdDr2_hat;
Jp_m2 = pi*Dm^4/32 + l2^2*Am2;
Jp2 = pi*Dr2^4/32 - Jp_m2;
Jp2b = pi*Dr2^4/32;

Ip2_hat = rho_hat*Jp2_hat*L2b_hat; %Momento polar do rotor 2
argumentos.Ip2_hat = Ip2_hat;
Ip2b_hat = rho_hat*Jp2b_hat*L2b_hat; %Momento polar do rotor 2
argumentos.Ip2b_hat = Ip2b_hat;
Ip2 = rho*Jp2*L2b; %Momento polar do rotor 2
argumentos.Ip2 = Ip2;
Ip2b = rho*Jp2b*L2b;
argumentos.Ip2b = Ip2b;

%Parâmetros do motor de potência:
%(1) Redução na saída do motor de potência:
nM = 7;
argumentos.nM = nM;
%(2) Inércia da bobina
Im = 4.698*10^-4; %kg*m^2
argumentos.Im = Im;
%(3) Inércia total na saída do motor:
Ip3 = 0;
Ip3b = Ip3 + nM^2*Im;
argumentos.Ip3bMt = Ip3b;
%(4) Amortecimento interno do motor:
cM = 1.784*10^-4; %Nm*s/rad
argumentos.cM = cM;

%Parâmetros do motor de passo:
%(1) Redução na saída do motor de passo:
nMs = 16;
argumentos.nMs = nMs;
%(2) Inércia da bobina:
Ims = 4.698*10^-4; %kg*m^2
argumentos.Ims = Ims;
%(4) Amortecimento interno do motor:
cMs = 1.784*10^-4; %Nm*s/rad
argumentos.cMs = cMs;

%Plataforma móvel:
l_pl_hat = ad*0.45; %comprimento lateral (m)
argumentos.l_pl_hat = l_pl_hat;
%dl_pl = percent_geral*l_pl_hat;
l_pl = ad*0.45;
e_pl_hat = 0.02; %espessura da plataforma (m)
argumentos.e_pl_hat = e_pl_hat;
%de_pl = percent_geral*e_pl_hat;
e_pl = 0.02;
M3_pl_hat = rho_hat*l_pl_hat^2*e_pl_hat;
argumentos.M3_pl_hat = M3_pl_hat;
M3_pl = rho*l_pl^2*e_pl;
argumentos.M3_pl = M3_pl;
I3_pl_hat = (1/12)*M3_pl_hat*l_pl_hat^2;
argumentos.I3_pl_hat = I3_pl_hat;
I3_pl = (1/12)*M3_pl*l_pl^2;
argumentos.I3_pl = I3_pl;
Ip3_pl = 2*I3_pl;
argumentos.Ip3_pl = Ip3_pl;

%Passo do fuso (movimento vertical):
pf_hat = 0.005; %m
argumentos.pf_hat = pf_hat;
dpf = percent_geral*pf_hat;
pf = 0.005; %m
betaS_hat = 2*pi/pf_hat;
argumentos.betaS_hat = betaS_hat;
dbetaSdpf_hat = -(2*pi)/pf_hat^2;
argumentos.dbetaSdpf_hat = dbetaSdpf_hat;
betaS = 2*pi/pf;
argumentos.betaS = betaS;

%Gravidade:
g_hat = 9.8;
argumentos.g_hat = g_hat;
%dg = percent_geral*g_hat;
g = 9.8; %m/s2
argumentos.g = g;

R_NM1_hat = 0.045;
argumentos.R_NM1_hat = R_NM1_hat; %m
R3_hat = R_NM1_hat;
argumentos.R3_hat = R3_hat;
dR = percent_geral*R_NM1_hat;
dR3 = dR;
R_NM1 = 0.045; %m
R3 = R_NM1;
argumentos.R_NM1 = R_NM1; %m
argumentos.R3 = R3; %m

%Parâmetros do modelo de atrito descontínuo na origem:
mi0_hat = 0.1;
argumentos.mi0_hat = mi0_hat;
dmi0 = percent_geral*mi0_hat;
mi0 = 0.1;
argumentos.mi0 = mi0;
mi1_hat = 0.1;
argumentos.mi1_hat = mi1_hat;
dmi1 = percent_geral*mi1_hat;
mi1 = 0.1;
argumentos.mi1 = mi1;
mi2_hat = 0.05;
argumentos.mi2_hat = mi2_hat;
dmi2 = percent_geral*mi2_hat;
mi2 = 0.05;
argumentos.mi2 = mi2;
beta1mi_hat = 1.5;
argumentos.beta1mi_hat = beta1mi_hat;
dbeta1mi = percent_geral*beta1mi_hat;
beta1mi = 1.5;
argumentos.beta1mi = beta1mi;
beta2mi_hat = 1.0;
argumentos.beta2mi_hat = beta2mi_hat;
dbeta2mi = percent_geral*beta2mi_hat;
beta2mi = 1.0;
argumentos.beta2mi = beta2mi;

%Parâmetros do modelo de atrito contínuo na origem:
switch opcao_atr
    case 1
        %Atrito aço-aço a seco:
        mic0_hat = 0.6;
        argumentos.mic0_hat = mic0_hat;
        mic0 = 0.6;
        argumentos.mic0 = mic0;
        lambdami_hat = (0.8 - 0.6)/0.6;
        argumentos.lambdami_hat = lambdami_hat;
        lambdami = (0.8 - 0.6)/0.6;
        argumentos.lambdami = lambdami;
        indPsiAtr = 2; %com atrito na inversão de PsiDif
        argumentos.indPsiAtr = indPsiAtr;
    case 2
        %Atrito aço-aço com lubrificação:
        mic0_hat = 0.05;
        argumentos.mic0_hat = mic0_hat;
        mic0 = 0.05;
        argumentos.mic0 = mic0;
        lambdami_hat = (0.1 - 0.05)/0.05;
        argumentos.lambdami_hat = lambdami_hat;
        lambdami = (0.1 - 0.05)/0.05;
        argumentos.lambdami = lambdami;
        indPsiAtr = 2; %com atrito na inversão de PsiDif
        argumentos.indPsiAtr = indPsiAtr;
    case 3
        %Sem atrito:
        mic0_hat = 0;
        argumentos.mic0_hat = mic0_hat;
        mic0 = 0;
        argumentos.mic0 = mic0;
        lambdami_hat = (0.8 - 0.6)/0.6;
        argumentos.lambdami_hat = lambdami_hat;
        lambdami = (0.8 - 0.6)/0.6;
        argumentos.lambdami = lambdami;
        indPsiAtr = 1; %sem atrito na inversão de PsiDif
        argumentos.indPsiAtr = indPsiAtr;
    otherwise
        disp('other value')
end
dmic0 = percent_geral*mic0_hat;
dlambdami = percent_geral*lambdami_hat;
w0_hat = 10^-2*pi; %velodade referência da função de atrito (rad/s)
argumentos.w0_hat = w0_hat;
dw0 = percent_geral*w0_hat;
w0 = 10^-2*pi;
argumentos.w0 = w0;
%nmi = 2;
%Obs.: 
%(1) d5fiAdmdnu5 apresentou o termo nu^(2*nMi - 5) no numerador
%(2) para nu=0, nu^(2*nMi - 5) é indefinido para nMi=2
nmi = 3;
argumentos.nmi = nmi;

%Elementos Finitos:
%1. Número de elementos:
%N = 12;
N = 12;
%Obs.: O valor mínimo de N é 3
argumentos.N = N;
%Comprimento total do eixo:
L_hat = 2.5; %m
argumentos.L_hat = L_hat;
dL = percent_geral*L_hat;
%Comprimento entre rotor intermediário e motor:
L1_hat = 0.8*L_hat; %m
argumentos.L1_hat = L1_hat;
dL1 = percent_geral*L1_hat;
%Comprimento entre rotor intermediário e rotor da extremidade:
L2_hat = L_hat - L1_hat; %m
%
%Comprimento entre a saída do motor e o ponto de amortecimento lateral:
deltaVW_hat = 0.5; %m
%Número de elementos em deltaVW:
N0 = floor(N*deltaVW_hat/L_hat)*double(floor(N*deltaVW_hat/L_hat)>=1) +...
    (1)*double(floor(N*deltaVW_hat/L_hat)<1);
argumentos.N0 = N0;
%N1 = ceil(N*L1/L); %número de elementos em L1
N1 = N0 +...
    floor(N*(L1_hat-deltaVW_hat)/L_hat)*double(floor(N*(L1_hat-deltaVW_hat)/L_hat)>=1) +...
    (1)*double(floor(N*(L1_hat-deltaVW_hat)/L_hat)<1);
%Número de elementos em L2:
N2 = N - N1;
if N2<1
    N1 = floor(N*L1_hat/L_hat);
    N2 = N - N1;
end
argumentos.N1 = N1;
%2. Malha:
X_hat = zeros(1,N+1);
for i = 1:N+1
    if i<N0+2
        X_hat(i) = (i-1)*deltaVW_hat/N0;
    elseif i<N1+2
        X_hat(i) = X_hat(N0+1) + (i-N0-1)*(L1_hat-deltaVW_hat)/(N1-N0);
    else
        X_hat(i) = X_hat(N1+1) + (i-N1-1)*L2_hat/N2;
    end
end
argumentos.X_hat = X_hat;
%}
%{
%Número de elementos em L1:
N1 = ceil(N*L1_hat/L_hat);
%Número de elementos em L2:
N2 = N - N1;
if N2<1
    N1 = floor(N*L1_hat/L_hat);
    N2 = N - N1;
end
argumentos.N1 = N1;
%2. Malha:
X_hat = zeros(1,N+1);
for i = 1:N+1
    if i<N1+2
        X_hat(i) = (i-1)*L1_hat/N1;
    else
        X_hat(i) = X_hat(N1+1) + (i-N1-1)*L2_hat/N2;
    end
end
argumentos.X_hat = X_hat;
%}

%Comprimento total do eixo:
L = 2.5; %m
argumentos.L = L;
%Comprimento entre rotor intermediário e motor:
L1 = 0.8*L; %m
argumentos.L1 = L1;
L2 = L - L1; %m (cumprimento entre rotor intermediário e rotor da extremidade)
N1 = ceil(N*L1/L); %número de elementos em L1
N2 = N - N1; %número de elementos em L2
if N2<1
    N1 = floor(N*L1/L);
    N2 = N - N1;
end
argumentos.N1 = N1;
%2. Malha:
X = zeros(1,N+1);
for i = 1:N+1
    if i<N1+2
        X(i) = (i-1)*L1/N1;
    else
        X(i) = X(N1+1) + (i-N1-1)*L2/N2;
    end
end
argumentos.X = X;
%}
%{
%Número de elementos em L1:
N1 = ceil(N*L1_hat/L_hat);
%Número de elementos em L2:
N2 = N - N1;
if N2<1
    N1 = floor(N*L1_hat/L_hat);
    N2 = N - N1;
end
argumentos.N1 = N1;
%2. Malha:
X = zeros(1,N+1);
for i = 1:N+1
    if i<N1+2
        X(i) = (i-1)*L1_hat/N1;
    else
        X(i) = X(N1+1) + (i-N1-1)*L2_hat/N2;
    end
end
argumentos.X = X;
%}

%Flambagem:
%kf = 0.7; %rotulado-engastado
%kf = 2; %livre-engastado
kf = 0.5; %engastado-engastado
fS = 2.0; %fator de segurança
P_cr = fS*(M2+m2)*g;
d_cr = (64*P_cr*(kf*L2)^2/(pi^3*E))^(1/4); %diâmetro crítico de flambagem
argumentos.d_cr = d_cr;

%Frequência natural lateral (rotor intermediário):
kW0 = 3*E*J*(1/L2^3 + 1/L1^3);
wW = sqrt(kW0/(M2+m2)); %rad/s
fw = (wW/(2*pi))*60; %RPM
argumentos.fw = fw;

%%
%Frequências naturais (parâmetros concentrados):
%Obs.: Parâmetros apenas de referência
%(1) Dinâmica longitudinal:
M2b0 = M2b + dm;
k1 = E*A/L1;
argumentos.k10 = k1;
k2 = E*A/(L - L1);
argumentos.k20 = k2;
lambda_u1 = (1/2)*((k1 + k2)/M2b0 + k2/M1 - sqrt(((k1 + k2)/M2b0 - k2/M1)^2 + 4*k2^2/(M1*M2b0)));
lambda_u2 = (1/2)*((k1 + k2)/M2b0 + k2/M1 + sqrt(((k1 + k2)/M2b0 - k2/M1)^2 + 4*k2^2/(M1*M2b0)));
wu1 = sqrt(lambda_u1);
wu2 = sqrt(lambda_u2);
freqNu1 = wu1/(2*pi);
argumentos.freqNu1 = freqNu1;
freqNu2 = wu2/(2*pi);
argumentos.freqNu2 = freqNu2;
%(2) Dinâmica torsional:
K10 = G*J/L1;
argumentos.K10 = K10;
K20 = G*J/(L - L1);
argumentos.K20 = K20;
I2b = Ip2;
I2b0 = I2b + dm*l2^2;
lambda_teta1 = (1/2)*((K10 + K20)/I2b0 + K20/Ip1 - sqrt(((K10 + K20)/I2b0 - K20/Ip1)^2 + 4*K20^2/(Ip1*I2b0)));
lambda_teta2 = (1/2)*((K10 + K20)/I2b0 + K20/Ip1 + sqrt(((K10 + K20)/I2b0 - K20/Ip1)^2 + 4*K20^2/(Ip1*I2b0)));
wTeta1 = sqrt(lambda_teta1); %rad/s
wTeta2 = sqrt(lambda_teta2); %rad/s
freqNteta1 = wTeta1/(2*pi); %Hz
argumentos.freqNteta1 = freqNteta1;
freqNteta2 = wTeta2/(2*pi); %Hz
argumentos.freqNteta2 = freqNteta2;

%%
%Parâmetros da matriz de amortecimento:

%psic2 = 0;
psic3 = 0.1;
%psib2 = 0;
psib3 = 0.1;
%Ib20_hat = Ip2b_hat + dm_hat*l2_hat^2;
%K10_hat = G_hat*Jp_hat/L1_hat;
K20_hat = G_hat*Jp_hat/(L_hat - L1_hat);
%K1_hat = K10_hat;
K2_hat = K20_hat;
c1 = nM^2*cM;
%c2_hat = 2*psic2*sqrt((K1_hat + K2_hat)/Ib20_hat);
%c2 = c2_hat;
c3_hat = 2*psic3*sqrt(K2_hat/Ip1_hat);
c3 = c3_hat;

%(2) Parâmetros de contato:
nk = 2;
alfa_kNM1 = 10^-nk;
E_Aluminio = 70*10^9; %Pa
A_celulaS = (8*10^-3)*(55*10^-3); %m2
L_celulaS = 80*10^-3; %m
kNM1_hat = alfa_kNM1*E_Aluminio*A_celulaS/L_celulaS;%%%%%%%%%%%%%%%%%%%%%%%
k3_hat = kNM1_hat;
pkNM1 = 5;
dkNM1 = percent_geral*pkNM1*kNM1_hat;
dk3 = dkNM1;
argumentos.dk3 = dk3;
argumentos.kNM1_hat = kNM1_hat;
argumentos.k3_hat = k3_hat;
kNM1 = alfa_kNM1*E_Aluminio*A_celulaS/L_celulaS;%%%%%%%%%%%%%%%%%%%%%%%%%%%
argumentos.kNM1 = kNM1;
k3 = kNM1;
argumentos.k3 = k3;
b1 = nMs^2*cMs;
%b2 = 2*psib2*sqrt((k1+k2)/(M2b + dm));
b3 = 2*psib3*sqrt((k2+k3)/M1);

kV = 5*kW0; %N/m
kW = 5*kW0; %N/M
psiV = 0.5;
psiW = 0.5;
bU1 = b1;
bTx1 = c1;
bV = 2*psiV*kV;
bW = 2*psiW*kW;
bUNM1 = b3;
bTxNM1 = c3;
argumentos.bU1 = bU1;
argumentos.bTx1 = bTx1;
argumentos.bV = bV;
argumentos.bW = bW;
argumentos.bUNM1 = bUNM1;
argumentos.bTxNM1 = bTxNM1;

%%
% Parâmetros do Rubbing:
fg = 0.006; %Folga (m)
argumentos.fg = fg;
%lbd = 1 + 10^-2; %adm (Apresentou melhor resultado, 01/09/2022)
lbd = 1 + 5*10^-3; %adm 
delta0 = lbd*fg;
argumentos.delta0 = delta0; %Avanço por deformação durante o impacto [+ folga] (m)
Af = 50*10^-3*5*10^-3;
Lf = 10*10^-3;
Kip = E*Af/Lf;
argumentos.Kip = Kip; %Rigidez de impacto (N/m)

%(a) Parâmetros do modelo de atrito contínuo na origem:
miIp0c = mic0;
argumentos.miIp0c = miIp0c;
v0Ipc = w0;
argumentos.v0Ipc = v0Ipc;
nmiIpc = nmi;
argumentos.nmiIpc = nmiIpc;
lambdamiIpc = lambdami;
argumentos.lambdamiIpc = lambdamiIpc;

%(b) Parâmetros do modelo de atrito descontínuo na origem:
mi0ip = mi0;
argumentos.mi0ip = mi0ip;
mi1ip = mi1;
argumentos.mi1ip = mi1ip;
mi2ip = mi2;
argumentos.mi2ip = mi2ip;
beta1miIp = beta1mi;
argumentos.beta1miIp = beta1miIp;
beta2miIp = beta2mi;
argumentos.beta2miIp = beta2miIp;

%%
%Configuração inicial:
dLmax = (M1*g/(E*A))*L + (M2*g/(E*A))*L1;
argumentos.dLmax = dLmax;
%betaL = 10;
%{
%(2) Parâmetros de contato:
nk = 2;
alfa_kNM1 = 10^-nk;
E_Aluminio = 70*10^9; %Pa
A_celulaS = (8*10^-3)*(55*10^-3); %m2
L_celulaS = 80*10^-3; %m
kNM1_hat = alfa_kNM1*E_Aluminio*A_celulaS/L_celulaS;%%%%%%%%%%%%%%%%%%%%%%%
k3_hat = kNM1_hat;
pkNM1 = 5;
dkNM1 = percent_geral*pkNM1*kNM1_hat;
dk3 = dkNM1;
argumentos.dk3 = dk3;
argumentos.kNM1_hat = kNM1_hat;
argumentos.k3_hat = k3_hat;
kNM1 = alfa_kNM1*E_Aluminio*A_celulaS/L_celulaS;%%%%%%%%%%%%%%%%%%%%%%%%%%%
argumentos.kNM1 = kNM1;
k3 = kNM1;
argumentos.k3 = k3;
%}

%Avanço máximo durante o contato:
%(a) Todo o eixo principal em regime de tração:
u31Ref = M1*g/k3;
argumentos.u31Ref = u31Ref;
%(b) Apenas a parte do eixo de comprimento L1 em regime de tração:
u32Ref = (M1 + M2b + dm)*g/k3;
argumentos.u32Ref = u32Ref;

xi_amortecimento = 0;
wC_NM1 = sqrt(kNM1/M1);
c_NM1 = 2*xi_amortecimento*M1*wC_NM1;
argumentos.c_NM1 = c_NM1;

wC_NM1_Hat = sqrt(kNM1_hat/M1);
c_NM1_Hat = 2*xi_amortecimento*M1*wC_NM1_Hat;
argumentos.c_NM1_Hat = c_NM1_Hat;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Equação de restrições "Gama*u = h":
%(1) Sem a restrição do movimento longitudinal no atrito (u_{N+1}):
Ng = 6*(N + 1);
GamaR = zeros(8,Ng);
GamaR(1,2) = 1;
GamaR(2,3) = 1;
GamaR(3,4) = 1;
GamaR(4,5) = 1;
GamaR(5,6*(N+1)-4) = 1;
GamaR(6,6*(N+1)-3) = 1;
GamaR(7,6*(N+1)-2) = 1;
GamaR(8,6*(N+1)-1) = 1;
argumentos.GamaR = GamaR;
hRestr = zeros(8,1);
argumentos.hRestr = hRestr;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Número de restrições no modelo em Elementos Finitos:
rEstr = 8;
%rEstr = 12;
argumentos.rEstr = rEstr;

%%
%Trajetória desejada: y1d = teta3pd(t) = wRef*fiAdm(t/tRef)
%Parâmetros da função adimensional
%Descrição: fiAdm(nu) = tanh(alfa*nu) - alfa*nu/(beta^2*nu^2 + 1)
%(i) fi(nu=0) = 0;
%(ii) dfidnu(nu=0) = 0;
%(iii) d2fidnu2(nu=0) = 0;
%(iv) d3fidnu3(nu=0) = 0;
%(v) d4fidnu4(nu=0) = 0;
fiOpt = @(nu,alfa)(tanh(alfa*nu) - (alfa*nu)./((alfa^2*nu.^2)/3 + ones(size(nu))));
CUSTO = @(alfa)((fiOpt(1,alfa) - 0.9)^2);
alfa0 = 28; %fi(1) == 0.9
alfa = fminsearch(CUSTO, alfa0);
argumentos.alfa_trajetoria_desejada = alfa;
t0_9wRef = 5; %s (instante em que teta3pd = 0.9*wRef)
argumentos.t0_9wRef = t0_9wRef;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARÂMETROS PRIMÁRIOS (ESTIMATIVAS E VARIABILIDADE)

dqsiTil = [dd; drho; dE; dnu;...
    dL1b; dDr1; dL2b; dDr2;...
    drhom2; dDm; dl2;...
    dmi0; dmi1; dmi2; dbeta1mi; dbeta2mi;...
    dmic0; dw0; dlambdami;...
    dR3;...
    dL; dL1;
    dpf];

qsiTil = [d_hat; rho_hat; E_hat; nu_hat; 
    L1b_hat; Dr1_hat; L2b_hat; Dr2_hat; 
    rhom2_hat; Dm_hat; l2_hat;
    mi0_hat; mi1_hat; mi2_hat; beta1mi_hat; beta2mi_hat; 
    mic0_hat; w0_hat; lambdami_hat; 
    R3_hat; 
    L_hat; L1_hat;
    pf_hat];

argumentos.dqsiTil = dqsiTil;
argumentos.qsiTil = qsiTil;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%OUTROS PARÂMETROS:

alfam2 = 45; %posição angular inicial de dm (em graus) no rotor 2
alfam2 = alfam2*pi/180; %posição angular inicial de dm (em rad)
argumentos.alfam2 = alfam2;

%Fator de suavização da função sinal (sign):
%Obs.: signS(nu) = tanh(nS*pi*nu);
nS = 10;
argumentos.nS = nS;

%Fator de suavização da função degrau:
%Obs.: f(x) = (x>=x0) ~ (1/2)*(tanh(nf*pi*(x - x0)) + 1);
nf = 10^8;
argumentos.nf = nf;

%Velocidade de operação:
%wRef = 120; %RPM
wRef = 30; %RPM
wRef = wRef*(pi/30); %rad/s
argumentos.wRef = wRef;

%Carga de operação:
%fUref = 0.6;
fUref = 0.95;
argumentos.fUref = fUref;
Nref = fUref*M1*g; %N
argumentos.Nref = Nref;

%%
%SIMULAÇÕES

switch simulacao
    case 1
        
    case 2
        
    case 3
        
    case 4
        
    case 5
        
    case 6
        
    case 7
        
    case 8
        
    case 9
        
    case 10
        
    case 11
        
    case 12
        
    case 13
        
    case 14
        
    case 15
        
    case 16    
    case 17    
    case 18    
    case 19    
    case 20    
    case 21    
    case 22    
    case 23    
    case 24    
    case 25    
    case 26    
    case 27    
    case 28    
    case 29
    case 30
    case 31
    case 32
    case 33
    case 34
    case 35
    case 36
    case 37
    case 38
    case 39
    case 40
    case 41
    otherwise
        disp('other value')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fator multiplicativo de FcLC e FdcLC: FcLC = sigmaF*abs(daC)
sigmaF = 1;
argumentos.sigmaF = sigmaF;
%Fator multiplicativo de Dlc: Dlc = sigmaD*(1/bHat)*abs(dbHat)
sigmaD = 1;
argumentos.sigmaD = sigmaD;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Aproximação de zero:
epsZero = 10^-15;
argumentos.epsZero = epsZero;

%%
%Newmark - parâmetros de integração:
%epsilon = 10^-15;
%epsilon = 10^-12;
epsilon = 10^-9;
%epsilon = 10^-3;
argumentos.epsNewtRaph = epsilon;
Nmax = 100;
argumentos.Nmax = Nmax;
%(1) Passo e tempo de integração:
%dt = 1*10^-5; %s
dt = 1*10^-4; %s
%dt = 1*10^-3; %s
argumentos.dt_newmark = dt;
t0 = 0; %s
argumentos.t0 = t0;
tf = 1.0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tf = 25.0; %s
argumentos.tf = tf;
T = (t0:dt:tf).';
argumentos.T = T;
nT = length(T);
argumentos.nT = nT;
%(2) Outros parâmetros:
gama_newmark = 1/2;
argumentos.gama_newmark = gama_newmark;
beta_newmark = 1/4;
argumentos.beta_newmark = beta_newmark;
k_newmark = (beta_newmark*dt^2)^-1;
argumentos.k_newmark = k_newmark;
end