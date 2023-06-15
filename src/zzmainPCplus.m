clear
clc
%%
%Simulações:
%(1) Dinâmica de torção em malha aberta
%    (a) Sem atrito
%(2) Dinâmica longitudinal em malha aberta
%    (a) Aplicação de PID na entrada
%(3) Dinâmica interna e camada limite para torção
%    (a) Sem atrito
%(4) Camada limite para dinâmica longitudinal
%(5) Dinâmica de torção com controlador
%    (a) Sem atrito
%(6) Dinâmica longitudinal com controlador
%    Obs.: Não deu certo
%    - Problemas de convergência com Newton-Raphson
%    - Alta frequência de chaveamento na entrada
%(7) Dinâmica longitudinal com controlador e perturbações singulares
%    (a) Trajetória desejada nula: regulador
%    Obs.: Não deu certo
%    - Problemas de inversão com mal condicionamento de matrizes
%(8) Dinâmica longitudinal com controlador e perturbações singulares
%    (a) Trajetória desejada não nula
%    Obs.: Não deu certo
%    - Problemas de convergência com Newton-Raphson
%    - Alta frequência de chaveamento na entrada
%(9) Dinâmica longitudinal em malha aberta (ref.: simulação 2)
%    (a) Aplicação de trajetória desejada diretamente na entrada
%    (b) Fu = y2d(t)
%    (c) Sem PID
%    Obs.: Bom resultado
%    - O resultado obtido por ode45 apresentou maior oscilação de y2
%    - A amplitude de oscilação foi bem baixa
%(10) Lei de controle para sistema sem acoplamento lateral (ref.: 5 e 9)
%    (a) Fu = y2d(t)
%    (b) Tx = Uc(t,Zt,mid,fiBd)
%    (c) Sem atrito
%(11) Lei de controle para sistema sem acoplamento lateral
%    (a) Fu = y2d(t)
%    (b) Tx = Uc(t,Zt,y2,...,y4p2,mid,fiBd,Z0)
%    (c) Com atrito
%(12) Dinâmica de torção com controlador 
%    (a) Planta: modelo com acoplamento de dinâmicas lateral e
%    longitudinal, com atrito
%    (b) Metamodelo: modelo sem acoplamento, sem atrito
%(13) Dinâmica de torção com controlador 
%    (a) Planta: modelo com acoplamento de dinâmicas lateral e longitudinal
%    (b) Metamodelo: modelo sem acoplamento na torção e no movimento
%    longitudinal
%    (c) Planta e metamodelo com atrito
%    (d) Entrada: FuUc(t,u1) = y2d(t) - Kf*(u1 - y2d(t)/kNM1)
%    (e) Suavização da curva de impácto no Rubbing
%(14) Dinâmica de torção com controlador (ref.: 12)
%    (a) Planta: modelo com acoplamento de dinâmicas lateral e longitudinal
%    (b) Metamodelo: modelo sem acoplamento na torção e movimento longitud.
%    (c) Planta e metamodelo sem atrito
%    (d) Entrada: FuUc(t,u1) = y2d(t) - Kf*(u1 - y2d(t)/kNM1)
%    (e) Suavização da curva de impácto no Rubbing
%(15) Dinâmica do sistema com acoplamento de dinâmicas:
%    (a) Sem atrito
%    (b) Malha aberta
%    (c) Entrada: FuUc(t) = y2d(t)
%(16) Dinâmica do sistema com acoplamento de dinâmicas:
%    (a) Sem atrito
%    (b) Malha aberta
%    (c) Entrada: FuUc(t,u1) = y2d(t) - Kf*(u1 - y2d(t)/kNM1)
%(17) Lei de controle para sistema sem acoplamento lateral (ref.: 10)
%    (a) Entrada: FuUc(t,u1) = y2d(t) - Kf*(u1 - y2d(t)/kNM1)
%    (b) Tx = Uc(t,Zt,mid,fiBd)
%    (c) Sem atrito
%(18) Dinâmica do sistema sem acoplamento (ref.: 1, 2):
%    (a) Entrada: FuUc(t,u1) = y2d(t) - Kf*(u1 - y2d(t)/kNM1)
%    (b) Com atrito
%    (c) Malha aberta em Tx
%(19) Modelo em Elementos Finitos (ref.: 18):
%    (a) Malha aberta em Tx
%    (b) Entrada: FuUc(u1) = -Kf*(ub1 - ub1Ref)
%    (c) Com ou sem atrito, a depender de indTatr (usar buscador)
%    (d) Restrições aplicadas de modo INTRÍNSECO
%(20) Dinâmica de torção com controlador (ref.: 14)
%    (a) Planta: modelo em Elementos Finitos
%    (b) Metamodelo: modelo sem acoplamento na torção e movimento longitud.
%    (c) Planta e metamodelo sem atrito
%    (d) Entrada: FuUc(t,u1) = y2d(t) - Kf*(u1 - y2d(t)/kNM1)
%    (e) Suavização da curva de impácto no Rubbing
%    (f) Restrições aplicadas de modo INTRÍNSECO
%(21) Modelo em Elementos Finitos (ref.: 19):
%    (a) Malha aberta em Tx
%    (b) Entrada: FuUc(u1) = -Kf*(ub1 - ub1Ref)
%    (c) Com ou sem atrito, a depender de indTatr (usar buscador)
%    (d) Restrições aplicadas de modo EXTRÍNSECO
%(22) Dinâmica de torção com controlador (ref.: 20, 21)
%    (a) Planta: modelo em Elementos Finitos
%    (b) Metamodelo: modelo sem acoplamento na torção e movimento longitud.
%    (c) Planta e metamodelo sem atrito
%    (d) Entrada: FuUc(t,u1) = y2d(t) - Kf*(u1 - y2d(t)/kNM1)
%    (e) Suavização da curva de impácto no Rubbing
%    (f) Restrições aplicadas de modo EXTRÍNSECO
%(23) Dinâmica de torção com controlador e restrição variavel em t (ref.: 22)
%    (a) Planta: modelo em Elementos Finitos
%    (b) Metamodelo: modelo sem acoplamento na torção e movimento longitud.
%    (c) Planta e metamodelo SEM atrito
%    (d) Suavização da curva de impácto no Rubbing
%    (e) Restrições aplicadas de modo EXTRÍNSECO
%    (f) Restrição variável no tempo em ub1: fhRestr(t) = [y2d(t)/kNM1; hRestr]
%(24) Lei de controle para sistema sem acoplamento lateral (ref.: 11)
%    (a) Entrada: FuUc(t,u1) = y2d(t) - Kf*(u1 - y2d(t)/kNM1)
%    (b) Tx = Uc(t,Zt,mid,fiBd)
%    (c) Com atrito
%(25) Modelo em Elementos Finitos (ref.: 21):
%    (a) Malha aberta via restrições variáveis no tempo
%    (b) Com ou sem atrito, a depender de indTatr (usar buscador)
%    (c) Restrições aplicadas de modo EXTRÍNSECO
%    (d) Restrição variável no tempo em ub1: fhU(t) = [hU(t); hRestr]
%    (e) Restrição variável no tempo em tetaX1: fhTeta(t) = [hTeta(t); hRestr]
%(26) Dinâmica de torção com controlador e restrição variavel em t (ref.: 13, 23)
%    (a) Planta: modelo em Elementos Finitos
%    (b) Metamodelo: modelo sem acoplamento na torção e movimento longitud.
%    (c) Planta e metamodelo COM atrito
%    (d) Suavização da curva de impácto no Rubbing
%    (e) Restrições aplicadas de modo EXTRÍNSECO
%    (f) Restrição variável no tempo em ub1: fhRestr(t) = [y2d(t)/kNM1; hRestr]
%(27) Modelo em Elementos Finitos (ref.: 25):
%    (a) Malha aberta via restrições variáveis no tempo
%    (b) Com ou sem atrito, a depender de indTatr (usar buscador)
%    (c) Restrições aplicadas de modo EXTRÍNSECO
%    (d) Restrição variável no tempo em ub1: fhU(t) = y2d(t)/kNM1
%    (e) tetaX1 como entrada do sistema: tetaX1(t), tetapX1(t), tetappX1(t)
%(28) Modelo em Elementos Finitos (ref.: 27): [UTILIZADO NA TESE]
%    (a) Malha aberta via restrições variáveis no tempo
%    (b) Com ou sem atrito, a depender de indTatr (usar buscador)
%    (c) Restrições aplicadas de modo EXTRÍNSECO
%    (d) Restrição variável no tempo em ub1: fhU(t) = y2d(t)/kNM1
%    (e) tetaX1 como entrada do sistema: tetaX1(t), tetapX1(t), tetappX1(t)
%    (f) avanço inicial sem rotação, até um deslocamento desejado
%(29) Dinâmica longitudinal e completa em modelos EF e PC: (ref.: 26)
%    (a) Modelo em Parâmetros Concentrados sem acoplamento
%    (b) Modelo em Parâmetros Concentrados com acoplamento
%    (c) Modelo em Elementos Finitos
%    (d) AUSÊNCIA DE ROTAÇÃO
%    (e) Restrições aplicadas de modo EXTRÍNSECO
%    (f) Restrição variável no tempo em ub1: fhRestr(t) = [y2d(t)/kNM1; 0; hRestr]
%    (g) Ausência de controle
%(30) Dinâmica em modelos EF e PC: (ref.: 29, 32)
%    (a) Modelo em Parâmetros Concentrados com acoplamento
%    (b) Modelo em Elementos Finitos
%    (c) COM ROTAÇÃO
%    (d) Restrições aplicadas de modo EXTRÍNSECO
%    (e) Restrição variável no tempo em ub1: fhRestr(t) = [y2d(t)/kNM1; inty1d(t); hRestr]
%    (f) Ausência de controle
%(31) Modelos reduzidos em malha aberta (Pk):
%    (a) Modelo linear
%    (b) Aplicação extrínseca de restrição (ub1 = hU(t))
%(32) Modelos reduzidos em malha aberta (Pk):
%    (a) Modelo linear
%    (b) Aplicação intrínseca de restrição (ub1 = hU(t), ubp1 = hUp(t), ubpp1 = hUpp(t))
%(33) Modelo alternativo (Tat aprox Tr):
%    (a) Entrada: Tx
%    Obs.: Não deu certo
%(34) Modelo alternativo (Tat aprox Tr):
%    (a) Entrada: teta1, tetap1, tetapp1
%    Obs.: Não deu certo
%(35) Dinâmica de torção com controlador e restrição variavel em t (ref.: 26)
%    (a) Planta: modelo em Elementos Finitos
%    (b) Metamodelo: modelo SEM acoplamento na dinâmica de torção
%    (c) Movimento longitudinal (y2): estimado a partir do modelo em EF
%    (c) Planta e metamodelo COM atrito
%    (d) Suavização da curva de impácto no Rubbing
%    (e) Restrições aplicadas de modo EXTRÍNSECO
%    (f) Restrição variável no tempo em ub1: fhRestr(t) = [y2d(t)/kNM1; hRestr]
%(36) Dinâmica de torção com controlador e restrição variavel em t (ref.: 35)
%    (a) Planta: modelo em Elementos Finitos
%    (b) Metamodelo: modelo SEM acoplamento na dinâmica de torção
%    (c) Movimento longitudinal (y2): estimado a partir do modelo em EF
%    (c) Planta e metamodelo COM atrito
%    (d) Suavização da curva de impácto no Rubbing
%    (e) Restrições aplicadas de modo EXTRÍNSECO
%    (f) Restrição variável no tempo em ub1: fhRestr(t) = [y2d(t)/kNM1; hRestr]
%    (g) Algoritmo para convergência de y2
%(37) Dinâmica de torção com controlador 
%    (a) Planta: modelo com acoplamento de dinâmicas lateral e longitudinal
%    (b) Metamodelo: modelo sem acoplamento na torção e no movimento
%    longitudinal
%    (c) Planta e metamodelo com atrito
%    (d) Restrição variável no tempo em ub1: hlU(t) ou hU(t)
%    (e) Suavização da curva de impácto no Rubbing
%(38) Dinâmica de torção com controlador (ref. 37 e 26) 
%    (a) Planta: modelo em EF
%    (b) Metamodelo: modelo sem acoplamento na torção e no movimento
%    longitudinal
%    (c) Planta e metamodelo com atrito
%    (d) Restrição variável no tempo em ub1: hlU(t) ou hU(t)
%    (e) Suavização da curva de impácto no Rubbing
%(39) Dinâmica de torção com controlador (ref. 38): [UTILIZADO NA TESE]
%    (a) Planta: modelo em EF
%    (b) Metamodelo: modelo sem acoplamento na torção e no movimento
%    longitudinal
%    (c) Planta e metamodelo com atrito
%    (d) Restrição variável no tempo em ub1: hlU(t) ou hU(t)
%    (e) Suavização da curva de impácto no Rubbing
%    (f) Melhoramentos na sim. 38: ângulos de Euler, correção de restrição
%(40) Dinâmica de torção com controlador (ref. 39) 
%    (a) Planta: modelo em EF
%    (b) Metamodelo: modelo sem acoplamento na torção e no movimento
%    longitudinal
%    (c) Planta e metamodelo com atrito
%    (d) Restrição variável no tempo em ub1 (metamodelo): hlU(t,yb1)
%    (e) Restrição variável no tempo em ub1 (sistema): hU(t)
%    (f) Suavização da curva de impácto no Rubbing
%    (g) Melhoramentos na sim. 38: ângulos de Euler, correção de restrição
%    (h) Implementação de media móvel yb1
%    OBS.: Mesmo com diferentes carregamentos longitudinal, não funcionou!!


%%
%Preprocessamento:
simulacao = 39;
opcao_m2 = 3;
%Obs.:
%(1) Valor preferencial: "opcao_m2 = 3"
%(2) Valor para identificação de parâmetros: "opcao_m2 = 1"
opcao_sign = 2;
opcaoSignAt = 1;
inc_dsignS = 1;
opcao_atr = 1;
ind_LVC = 0;
[ argumentos ] = zzpreprocess( simulacao, opcao_m2, opcao_sign, inc_dsignS,...
    opcao_atr, ind_LVC );
%Tamanho do vetor para cálculo da média (dinâmicas de alta frequência):
MEANativo = true;
nPD = 10;
ODEativo = true;
%Inclui (indTatr = 1) ou não (indTatr <> 1) o torque de atrito no mod EF:
indTatr = 1;
argumentos.indTatr = indTatr;
argumentos.opcaoSignAt = opcaoSignAt;
%Indicador de ajuste do metamodelo:
IndAjustMod = 5;
%Obs.:
%(1) IndAjustMod = 'sa': sem ajuste
%(2) IndAjustMod = 1: ajuste 1 (apenas rigidezes)
%(3) IndAjustMod = 2: ajuste 2 (apenas rigidezes)
%(4) IndAjustMod = 3: ajuste 3 (apenas rigidezes)
%(5) IndAjustMod = 4: ajuste 4 (rigidezes e inércias)
%(6) IndAjustMod = 5: ajuste 5 (rigidezes e inércias, IpTotal e Mtotal fixos)
%(7) IndAjustMod = 6: ajuste 6 (mmod's separados a partir de modEF: tX e u)
%    Obs.: Ainda em implementação


%%
%Parâmetros:
Ib2 = argumentos.Ip2b;
I3 = argumentos.Ip3bMt;
M1 = argumentos.M1rt;
Mb2 = argumentos.M2b;
M3 = argumentos.M3_pl;
g = argumentos.g;
G = argumentos.G;
Jp = argumentos.Jp;
E = argumentos.E;
A = argumentos.Area;
L = argumentos.L;
L1 = argumentos.L1;
l2 = argumentos.l2;
dm = argumentos.dm;
b1 = argumentos.bU1;
c1 = argumentos.bTx1;
b3 = argumentos.bUNM1;
c3 = argumentos.bTxNM1;
R3 = argumentos.R3;
alfaAc = argumentos.alfam2;
rho = argumentos.rho;
J = Jp/2;
nGL = argumentos.nGL;
GamaR = argumentos.GamaR;
hRestr = argumentos.hRestr;
rEstr = argumentos.rEstr;
kNmk = argumentos.k_newmark;
Kip = argumentos.Kip;
Dr2 = argumentos.Dr2;
fg = argumentos.fg;

%Acrescimo de deslocamento longitudinal na traj des acopl y2d:
Nref = argumentos.Nref;
kNM1 = argumentos.kNM1;

nUaj = 35;
uNM1aj = nUaj*(Nref/kNM1); %m
nFg = 1;
uNM1aj0 = L - (sqrt(L1^2 - (fg/nFg)^2) + sqrt((L-L1)^2 - (fg/nFg)^2));
%uNM1aj = uNM1aj0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modelo em Elementos Finitos
Ip1 = argumentos.Ip1;
I1 = Ip1/2;
M2b = Mb2;
Ip2b = Ib2;
I2 = Ip2b/2;
M3pl = M3;
I3pl = argumentos.I3_pl;
Ip3pl = argumentos.Ip3_pl;
Ip3bMt = argumentos.Ip3bMt;
N = argumentos.N;
N1 = argumentos.N1;
R = argumentos.raio;
ks = argumentos.ks;
bU1 = argumentos.bU1;
bUNM1 = argumentos.bUNM1;
bTx1 = argumentos.bTx1;
bTxNM1 = argumentos.bTxNM1;
R_NM1 = argumentos.R_NM1;

n1 = @(xi, le)(1 - xi);
n2 = @(xi, le)(xi);
H1 = @(xi, le)(1 - 3*xi^2 + 2*xi^3);
H2 = @(xi, le)(le*(xi - 2*xi^2 + xi^3));
H3 = @(xi, le)(3*xi^2 - 2*xi^3);
H4 = @(xi, le)(le*(-xi^2 + xi^3));
%ok, 27/02/2020
dn1dx = @(xi, le)(-1/le);
dn2dx = @(xi, le)(1/le);
dH1dx = @(xi, le)((6*xi*(xi - 1))/le);
dH2dx = @(xi, le)(3*xi^2 - 4*xi + 1);
dH3dx = @(xi, le)(-(6*xi*(xi - 1))/le);
dH4dx = @(xi, le)(xi*(3*xi - 2));
%ok, 27/02/2020
d2H1dx2 = @(xi, le)((6*(2*xi - 1))/le^2);
d2H2dx2 = @(xi, le)((6*xi - 4)/le);
d2H3dx2 = @(xi, le)(-(6*(2*xi - 1))/le^2);
d2H4dx2 = @(xi, le)((6*xi - 2)/le);
%ok, 27/02/2020
Nu = @(xi, le)([n1(xi, le), 0, 0, 0, 0, 0,...
    n2(xi, le), 0, 0, 0, 0, 0]);
repositorio.Nu = Nu;
Nv = @(xi, le)([0, H1(xi, le), H2(xi, le), 0, 0, 0,...
    0, H3(xi, le), H4(xi, le), 0, 0, 0]);
repositorio.Nv = Nv;
NtetaZ = @(xi, le)([0, dH1dx(xi, le), dH2dx(xi, le), 0, 0, 0,...
    0, dH3dx(xi, le), dH4dx(xi, le), 0, 0, 0]);
repositorio.NtetaZ = NtetaZ;
Nw = @(xi, le)([0, 0, 0, H1(xi, le), -H2(xi, le), 0,...
    0, 0, 0, H3(xi, le), -H4(xi, le), 0]);
repositorio.Nw = Nw;
NtetaY = @(xi, le)([0, 0, 0, -dH1dx(xi, le), dH2dx(xi, le), 0,...
    0, 0, 0, -dH3dx(xi, le), dH4dx(xi, le), 0]);
repositorio.NtetaY = NtetaY;
NtetaX = @(xi, le)([0, 0, 0, 0, 0, n1(xi, le),...
    0, 0, 0, 0, 0, n2(xi, le)]);
repositorio.NtetaX = NtetaX;
%ok, 27/02/2020
dNu_dx = @(xi, le)([dn1dx(xi, le), 0, 0, 0, 0, 0,...
    dn2dx(xi, le), 0, 0, 0, 0, 0]);
repositorio.dNu_dx = dNu_dx;
dNv_dx = @(xi, le)([0, dH1dx(xi, le), dH2dx(xi, le), 0, 0, 0,...
    0, dH3dx(xi, le), dH4dx(xi, le), 0, 0, 0]);
repositorio.dNv_dx = dNv_dx;
dNtetaZ_dx = @(xi, le)([0, d2H1dx2(xi, le), d2H2dx2(xi, le), 0, 0, 0,...
    0, d2H3dx2(xi, le), d2H4dx2(xi, le), 0, 0, 0]);
repositorio.dNtetaZ_dx = dNtetaZ_dx;
dNw_dx = @(xi, le)([0, 0, 0, dH1dx(xi, le), -dH2dx(xi, le), 0,...
    0, 0, 0, dH3dx(xi, le), -dH4dx(xi, le), 0]);
repositorio.dNw_dx = dNw_dx;
dNtetaY_dx = @(xi, le)([0, 0, 0, -d2H1dx2(xi, le), d2H2dx2(xi, le), 0,...
    0, 0, 0, -d2H3dx2(xi, le), d2H4dx2(xi, le), 0]);
repositorio.dNtetaY_dx = dNtetaY_dx;
dNtetaX_dx = @(xi, le)([0, 0, 0, 0, 0, dn1dx(xi, le),...
    0, 0, 0, 0, 0, dn2dx(xi, le)]);
repositorio.dNtetaX_dx = dNtetaX_dx;
%ok, 27/02/2020

%u_aprox = @(xi, le, uj)(Nu(xi, le)*uj);
%v_aprox = @(xi, le, uj)(Nv(xi, le)*uj);
tetaZ = @(xi, le, uj)(NtetaZ(xi, le)*uj);
%w_aprox = @(xi, le, uj)(Nw(xi, le)*uj);
tetaY = @(xi, le, uj)(NtetaY(xi, le)*uj);
tetaX = @(xi, le, uj)(NtetaX(xi, le)*uj);
%ok, 27/02/2020
up = @(xi, le, upj)(Nu(xi, le)*upj);
vp = @(xi, le, upj)(Nv(xi, le)*upj);
tetaZp = @(xi, le, upj)(NtetaZ(xi, le)*upj);
wp = @(xi, le, upj)(Nw(xi, le)*upj);
tetaYp = @(xi, le, upj)(NtetaY(xi, le)*upj);
tetaXp = @(xi, le, upj)(NtetaX(xi, le)*upj);
%ok, 27/02/2020
upp = @(xi, le, uppj)(Nu(xi, le)*uppj);
vpp = @(xi, le, uppj)(Nv(xi, le)*uppj);
tetaZpp = @(xi, le, uppj)(NtetaZ(xi, le)*uppj);
wpp = @(xi, le, uppj)(Nw(xi, le)*uppj);
tetaYpp = @(xi, le, uppj)(NtetaY(xi, le)*uppj);
tetaXpp = @(xi, le, uppj)(NtetaX(xi, le)*uppj);
%ok, 05/03/2020
du_dx = @(xi, le, uj)(dNu_dx(xi, le)*uj);
dv_dx = @(xi, le, uj)(dNv_dx(xi, le)*uj);
dtetaZ_dx = @(xi, le, uj)(dNtetaZ_dx(xi, le)*uj);
dw_dx = @(xi, le, uj)(dNw_dx(xi, le)*uj);
dtetaY_dx = @(xi, le, uj)(dNtetaY_dx(xi, le)*uj);
dtetaX_dx = @(xi, le, uj)(dNtetaX_dx(xi, le)*uj);
%ok, 05/03/2020

%%
%(3) Matriz M1 (para termos elementares e concentrados):
%(3.1) Eixo e rotores 1 e 2:
T_m1j_Aux = @(xi, le, rhoA_ou_M, rhoJp_ou_Ip, rhoJ_ou_I)(2*rhoA_ou_M*...
    (Nu(xi, le).'*Nu(xi, le) + Nv(xi, le).'*Nv(xi, le) + Nw(xi, le).'*Nw(xi, le)) +...
    2*rhoJp_ou_Ip*...
    (NtetaX(xi, le).'*NtetaX(xi, le)) +...
    2*rhoJ_ou_I*...
    (NtetaY(xi,le).'*NtetaY(xi,le) + NtetaZ(xi,le).'*NtetaZ(xi,le)));
%repositorio.T_m1j_Aux = T_m1j_Aux;
%ok! 21/12/2021

%(3.2) Excentricidade:
%Termo a ser subtraído de M2 e acrescido a M1, de modo a que M2(0)=0:
T_dmj_dm_Aux = @(xi, le)(2*dm*...
    (l2*((NtetaX(xi, le).')*(Nw(xi, le)*cos(alfaAc)) +...
           (Nw(xi, le).'*cos(alfaAc))*(NtetaX(xi, le)))) +...
    dm*...
    (l2*(((Nu(xi, le).')*(NtetaY(xi, le)*sin(alfaAc) -...
                    NtetaZ(xi, le)*cos(alfaAc))) +...
           ((NtetaY(xi, le).'*sin(alfaAc) - NtetaZ(xi, le).'*cos(alfaAc))*(Nu(xi, le))))) +...
    dm*...
    (l2^2*((NtetaY(xi, le).'*sin(alfaAc) - NtetaZ(xi, le).'*cos(alfaAc))*(NtetaY(xi, le)*sin(alfaAc) - NtetaZ(xi, le)*cos(alfaAc))) -...
     l2*(Nv(xi, le).'*(NtetaX(xi, le)*sin(alfaAc)) +...
        (NtetaX(xi, le).'*sin(alfaAc))*Nv(xi, le))));
%Parcela da matriz M1 referente à massa m2:
T_m1j_m2_Aux = @(xi, le)(2*dm*...
    (Nu(xi, le).'*Nu(xi, le) + Nv(xi, le).'*Nv(xi, le) + Nw(xi, le).'*Nw(xi, le) +...
     l2^2*(NtetaX(xi, le).'*NtetaX(xi, le))) +...
     T_dmj_dm_Aux(xi, le));
%repositorio.T_m1j_m2_Aux = T_m1j_m2_Aux;
%ok! 22/12/2021

%(3.3) Parcela da matriz M1 referente à plataforma móvel:
T_m1j_M3_Aux = @(xi, le)(2*...
    (M3pl*((Nu(xi, le).'*Nu(xi, le)) + (Nv(xi, le).'*Nv(xi, le)) + (Nw(xi, le).'*Nw(xi, le))) +...
     I3pl*(NtetaY(xi, le).'*NtetaY(xi, le)) +...
     I3pl*(NtetaZ(xi, le).'*NtetaZ(xi, le))));
%ok! 24/12/2021

%(3.4) Parcela da matriz M1 referente à inércia interna do motor:
T_m1j_Motor_Aux = @(xi, le)(2*Ip3bMt*(NtetaX(xi, le).'*NtetaX(xi, le)));
%ok! 24/12/2021

%%
%(4) Matriz M2upp (para termos elementares e concentrados):
%(4.1) Eixo e rotores 1 e 2:
T_m2j_Aux = @(xi, le, uj, rhoJp_ou_Ip)(2*rhoJp_ou_Ip*...
    (tetaZ(xi,le,uj)*...
    (NtetaX(xi,le).'*NtetaY(xi,le) + NtetaY(xi,le).'*NtetaX(xi,le)) +...
    tetaZ(xi,le,uj)^2*...
    (NtetaY(xi,le).'*NtetaY(xi,le))));
%repositorio.T_m2j_Aux = T_m2j_Aux;
%ok! 21/12/2021
T_m2juppj_Aux = @(xi, le, uj, uppj, rhoJp_ou_Ip)(2*rhoJp_ou_Ip*...
    (tetaZ(xi,le,uj)*...
    (NtetaX(xi,le).'*tetaYpp(xi,le,uppj) + NtetaY(xi,le).'*tetaXpp(xi,le,uppj)) +...
    tetaZ(xi,le,uj)^2*...
    (NtetaY(xi,le).'*tetaYpp(xi,le,uppj))));
%repositorio.T_m2juppj_Aux = T_m2juppj_Aux;
%ok! 21/12/2021
T_dm2juppj_dujT_Aux = @(xi, le, uj, uppj, rhoJp_ou_Ip)(2*rhoJp_ou_Ip*...
    (tetaYpp(xi,le,uppj)*...
     (NtetaX(xi,le).'*NtetaZ(xi,le) +...
      2*tetaZ(xi,le,uj)*NtetaY(xi,le).'*NtetaZ(xi,le)) +...
     tetaXpp(xi,le,uppj)*...
     (NtetaY(xi,le).'*NtetaZ(xi,le))));
%repositorio.T_m2juppj_Aux = T_dm2juppj_dujT_Aux;
%ok! 21/12/2021

%(4.2) Excentricidade:
T_m2j_m2_Aux = @(xi, le, uj)(dm*...
    (l2^2*(2*tetaZ(xi, le, uj)*(NtetaX(xi, le).'*NtetaY(xi, le) + NtetaY(xi, le).'*NtetaX(xi, le)) +...
           2*tetaZ(xi, le, uj)^2*(NtetaY(xi, le).'*NtetaY(xi, le))) +...
     (tetaZ(xi, le, uj)^2*(2*(Nv(xi, le).'*Nv(xi, le)) + 2*(Nu(xi, le).'*Nu(xi, le)))) +...
     (tetaY(xi, le, uj)^2*(2*(Nw(xi, le).'*Nw(xi, le)) + 2*(Nu(xi, le).'*Nu(xi, le))))) +...
    dm*...
    (2*l2*(((NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*(((tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc))*Nu(xi, le)) +...
                                         (cos(tetaX(xi, le, uj) + alfaAc)*Nw(xi, le)))) +...
           ((Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) + Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)))*(NtetaX(xi, le) + (NtetaY(xi, le)*tetaZ(xi, le, uj)))))) +...
    dm*...
    (2*l2*(((Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*((NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc)) -...
                                                             (NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc)))) +...
           ((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))*(Nu(xi, le) - (Nw(xi, le)*tetaY(xi, le, uj)))))) +...
    dm*...
    (l2^2*2*((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc))) -...
     2*(Nv(xi, le).'*((Nw(xi, le)*tetaZ(xi, le, uj)*tetaY(xi, le, uj)) +...
              (tetaZ(xi, le, uj)*NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc))*l2 +...
              (NtetaX(xi, le)*sin(tetaX(xi, le, uj) + alfaAc))*l2) +...
        ((Nw(xi, le).'*tetaZ(xi, le, uj)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*l2*NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) + l2*NtetaX(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc))*Nv(xi, le)))) -...
    T_dmj_dm_Aux(xi, le));
%repositorio.T_m2j_m2_Aux = T_m2j_m2_Aux;
%ok! 22/12/2021
T_m2juppj_m2_Aux = @(xi, le, uj, uppj)(...
    dm*...
    (l2^2*(2*tetaZ(xi, le, uj)*(NtetaX(xi, le).'*tetaYpp(xi, le, uppj) + tetaXpp(xi, le, uppj)*NtetaY(xi, le).') +...
           2*tetaYpp(xi, le, uppj)*tetaZ(xi, le, uj)^2*NtetaY(xi, le).') +...
     (tetaZ(xi, le, uj)^2*(2*vpp(xi, le, uppj)*Nv(xi, le).' + 2*upp(xi, le, uppj)*Nu(xi, le).')) +...
     (tetaY(xi, le, uj)^2*(2*wpp(xi, le, uppj)*Nw(xi, le).' + 2*upp(xi, le, uppj)*Nu(xi, le).'))) +...
    dm*...
    (2*l2*(((NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((upp(xi, le, uppj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc))) +...
                                                     (wpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc)))) +...
           ((tetaXpp(xi, le, uppj) + (tetaYpp(xi, le, uppj)*tetaZ(xi, le, uj)))*(Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) +...
                                                                                 Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) +...
                                                                                               tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)))))) +...
    dm*...
    (2*l2*(((Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*((tetaYpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfaAc)) -...
                                                             (tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc)))) +...
           ((upp(xi, le, uppj) - (wpp(xi, le, uppj)*tetaY(xi, le, uj)))*(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) -...
                                                                         NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))))) +...
    dm*...
    (l2^2*2*(((tetaYpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfaAc)) - (tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc)))*(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))) -...
     2*(Nv(xi, le).'*((wpp(xi, le, uppj)*tetaZ(xi, le, uj)*tetaY(xi, le, uj)) +...
                      (tetaZ(xi, le, uj)*tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc))*l2 +...
                      (tetaXpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfaAc))*l2) +...
        (vpp(xi, le, uppj)*(Nw(xi, le).'*tetaZ(xi, le, uj)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*l2*NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) + l2*NtetaX(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc))))) -...
    T_dmj_dm_Aux(xi, le)*uppj);
%repositorio.T_m2juppj_m2_Aux = T_m2juppj_m2_Aux;
%ok! 22/12/2021
T_dm2juppjdujT_m2_Aux = @(xi, le, uj, uppj)(dm*...
    (l2^2*(2*(NtetaX(xi, le).'*tetaYpp(xi, le, uppj) + tetaXpp(xi, le, uppj)*NtetaY(xi, le).')*NtetaZ(xi, le) +...
           2*tetaYpp(xi, le, uppj)*NtetaY(xi, le).'*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)) +...
     ((2*vpp(xi, le, uppj)*Nv(xi, le).' + 2*upp(xi, le, uppj)*Nu(xi, le).')*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)) +...
     ((2*wpp(xi, le, uppj)*Nw(xi, le).' + 2*upp(xi, le, uppj)*Nu(xi, le).')*2*tetaY(xi, le, uj)*NtetaY(xi, le))) +...
    dm*...
    (2*l2*((((NtetaY(xi, le).'*NtetaZ(xi, le))*((upp(xi, le, uppj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc))) +...
                                                (wpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc))) +...
             (NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((upp(xi, le, uppj)*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) + (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)))) +...
                                                                      (-wpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))))) +...
           ((Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) + Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)))*(tetaYpp(xi, le, uppj)*NtetaZ(xi, le)) +...
            (tetaXpp(xi, le, uppj) + (tetaYpp(xi, le, uppj)*tetaZ(xi, le, uj)))*(-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + Nu(xi, le).'*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) + (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))))))) +...
    dm*...
    (2*l2*(((-Nw(xi, le).'*NtetaY(xi, le))*((tetaYpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfaAc)) -...
                                            (tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc))) +...
            (Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*(tetaYpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) +...
                                                             tetaZpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))) +...
           ((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))*(-wpp(xi, le, uppj)*NtetaY(xi, le)) + (upp(xi, le, uppj) - (wpp(xi, le, uppj)*tetaY(xi, le, uj)))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))))) +...
    dm*...
    (l2^2*2*((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))*(tetaYpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + tetaZpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
             (tetaYpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))) -...
     2*(Nv(xi, le).'*(wpp(xi, le, uppj)*(NtetaZ(xi, le)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le)) +...
                      tetaZpp(xi, le, uppj)*(NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*l2 +...
              tetaXpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*l2) +...
        (vpp(xi, le, uppj)*(Nw(xi, le).'*(NtetaZ(xi, le)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le)) +...
              NtetaZ(xi, le).'*(NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*l2 +...
              NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*l2)))));
%repositorio.T_dm2juppj_dujT_m2_Aux = T_dm2juppj_dujT_m2_Aux;
%ok! 22/12/2021

%(4.3) Plataforma móvel:
T_m2j_M3_Aux = @(xi, le, uj)(2*Ip3pl*...
    tetaZ(xi,le,uj)^2*(NtetaY(xi,le).'*NtetaY(xi,le)));
%repositorio.T_m2j_M3_Aux = T_m2j_M3_Aux;
%ok! 27/12/2021
T_m2uppj_M3_Aux = @(xi, le, uj, uppj)(2*Ip3pl*...
    tetaZ(xi,le,uj)^2*(NtetaY(xi,le).'*tetaYpp(xi,le,uppj)));
%repositorio.T_m2uppj_M3_Aux = T_m2uppj_M3_Aux;
%ok! 27/12/2021
T_dm2juppjdujT_M3_Aux = @(xi, le, uj, uppj)(2*Ip3pl*...
    (NtetaY(xi,le).'*tetaYpp(xi,le,uppj))*2*tetaZ(xi,le,uj)*NtetaZ(xi,le));
%repositorioT_dm2juppj_dujT_M3_Aux = T_dm2juppj_dujT_M3_Aux;
%ok! 27/12/2021

%%
%(5) Vetor giroscópico Gc (para termos elementares e concentrados):
%(5.1) Eixo e rotores 1 e 2:
T_gCj_Aux = @(xi, le, uj, upj, rhoJp_ou_Ip)(2*rhoJp_ou_Ip*tetaYp(xi, le, upj)*...
    (NtetaX(xi, le).'*tetaZp(xi, le, upj)*tetaYp(xi, le, upj) +...
     NtetaY(xi, le).'*tetaZp(xi, le, upj)*(tetaXp(xi, le, upj) + 2*tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) -...
     NtetaZ(xi, le).'*tetaYp(xi, le, upj)*(tetaXp(xi, le, upj) + tetaZ(xi, le, uj)*tetaYp(xi, le, upj))));
%repositorio.T_gCj_Aux = T_gCj_Aux;
%ok, 21/12/2021
T_dgCjdujT_Aux = @(xi, le, upj, rhoJp_ou_Ip)(2*rhoJp_ou_Ip*...
    (NtetaY(xi, le).'*tetaZp(xi, le, upj)*(2*NtetaZ(xi, le)*tetaYp(xi, le, upj)) -...
     NtetaZ(xi, le).'*tetaYp(xi, le, upj)*(NtetaZ(xi, le)*tetaYp(xi, le, upj))));
%repositorio.T_dgCj_dujT_Aux = T_dgCj_dujT_Aux;
%ok, 21/12/2021
T_dgCjdupjT_Aux = @(xi, le, uj, upj, rhoJp_ou_Ip)(2*rhoJp_ou_Ip*...
    (NtetaX(xi, le).'*(NtetaZ(xi, le)*tetaYp(xi, le, upj) + tetaZp(xi, le, upj)*NtetaY(xi, le)) +...
     NtetaY(xi, le).'*(NtetaZ(xi, le)*(tetaXp(xi, le, upj) + 2*tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) +...
                       tetaZp(xi, le, upj)*(NtetaX(xi, le) + 2*tetaZ(xi, le, uj)*NtetaY(xi, le))) -...
     NtetaZ(xi, le).'*(NtetaY(xi, le)*(tetaXp(xi, le, upj) + tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) +...
                       tetaYp(xi, le, upj)*(NtetaX(xi, le) + tetaZ(xi, le, uj)*NtetaY(xi, le)))));
%repositorio.T_dgCj_dupjT_Aux = T_dgCj_dupjT_Aux;
%ok, 21/12/2021

%(5.2) Excentricidade:
T_gCj_m2_Aux = @(xi, le, uj, upj)(dm*...
    (l2^2*(2*(tetaZp(xi, le, upj)*(NtetaX(xi, le).'*tetaYp(xi, le, upj) + tetaXp(xi, le, upj)*NtetaY(xi, le).')) +...
           2*(2*tetaZ(xi, le, uj)*tetaYp(xi, le, upj)*tetaZp(xi, le, upj))*NtetaY(xi, le).') +...
     (2*tetaZ(xi, le, uj)*tetaZp(xi, le, upj)*(2*vp(xi, le, upj)*Nv(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')) +...
     (2*tetaY(xi, le, uj)*tetaYp(xi, le, upj)*(2*wp(xi, le, upj)*Nw(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).'))) +...
    dm*...
    (2*l2*(((NtetaY(xi, le).'*tetaZp(xi, le, upj))*(up(xi, le, upj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)) + wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)) +...
            (NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((up(xi, le, upj)*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
                                              (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)))) +...
                                         (-wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)))) +...
           (((tetaYp(xi, le, upj)*tetaZp(xi, le, upj)))*(Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) + Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc))) +...
            (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) +...
                                     Nu(xi, le).'*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
                                           (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))))))) +...
    dm*...
    (2*l2*(((-Nw(xi, le).'*tetaYp(xi, le, upj))*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)) +...
            (Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*((tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) -...
                                 (-tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)))) +...
           ((-(wp(xi, le, upj)*tetaYp(xi, le, upj)))*(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)) +...
            (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))))) +...
    dm*...
    (l2^2*2*(((tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) - (-tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)))*(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)) +...
             (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))) -...
     2*(Nv(xi, le).'*((wp(xi, le, upj)*tetaZp(xi, le, upj)*tetaY(xi, le, uj) + wp(xi, le, upj)*tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) +...
              (tetaZp(xi, le, upj)^2*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))*l2 +...
              (tetaXp(xi, le, upj)^2*cos(tetaX(xi, le, uj) + alfaAc))*l2) +...
        (vp(xi, le, upj)*(Nw(xi, le).'*(tetaZp(xi, le, upj)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) +...
             l2*NtetaZ(xi, le).'*(tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
             l2*NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))))) -...
    dm*(l2^2*(2*tetaXp(xi, le, upj)*tetaYp(xi, le, upj) + tetaYp(xi, le, upj)^2*2*tetaZ(xi, le, uj))*NtetaZ(xi, le).' +...
    2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*(vp(xi, le, upj)^2 + up(xi, le, upj)^2) +...
    2*tetaY(xi, le, uj)*NtetaY(xi, le).'*(wp(xi, le, upj)^2 + up(xi, le, upj)^2) +...
    2*l2*((tetaYp(xi, le, upj)*NtetaZ(xi, le).')*(up(xi, le, upj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)) + wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)) +...
          (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(up(xi, le, upj)*((NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).') +...
                                       (NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')) - wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')) +...
    2*l2*((-wp(xi, le, upj)*NtetaY(xi, le).')*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)) +...
          (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc))*NtetaX(xi, le).') +...
    l2^2*2*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc))*NtetaX(xi, le).' -...
    2*vp(xi, le, upj)*(wp(xi, le, upj)*(NtetaZ(xi, le).'*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le).') +...
    l2*tetaZp(xi, le, upj)*(NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).') +...
    l2*tetaXp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')));
%repositorio.T_gCj_m2_Aux = T_gCj_m2_Aux;
T_dgCjdujT_m2_Aux = @(xi, le, uj, upj)(dm*...
    (l2^2*(2*(2*tetaYp(xi, le, upj)*tetaZp(xi, le, upj))*(NtetaY(xi, le).'*NtetaZ(xi, le))) +...
     (2*tetaZp(xi, le, upj)*(2*vp(xi, le, upj)*Nv(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')*NtetaZ(xi, le)) +...
     (2*tetaYp(xi, le, upj)*(2*wp(xi, le, upj)*Nw(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')*NtetaY(xi, le))) +...
    dm*...
    (2*l2*(((NtetaY(xi, le).'*tetaZp(xi, le, upj))*(up(xi, le, upj)*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
                                   (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))) -...
                               wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
            ((NtetaY(xi, le).'*NtetaZ(xi, le))*((up(xi, le, upj)*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
                                              (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)))) +...
                                (-wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))) +...
             (NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((up(xi, le, upj)*((tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + (NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*tetaXp(xi, le, upj)) +...
                                               (-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) - (NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*tetaXp(xi, le, upj)))) +...
                                          (-wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj))))) +...
           ((tetaYp(xi, le, upj)*tetaZp(xi, le, upj))*(-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) +...
                             Nu(xi, le).'*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
                                   (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)))) +...
            ((-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) +...
              Nu(xi, le).'*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
                    (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))))*(tetaYp(xi, le, upj)*NtetaZ(xi, le)) +...
             (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(-Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj) +...
                                      Nu(xi, le).'*((tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + (NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*tetaXp(xi, le, upj)) +...
                                            (-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) - (NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*tetaXp(xi, le, upj)))))))) +...
    dm*...
    (2*l2*(((-Nw(xi, le).'*tetaYp(xi, le, upj))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
            ((-Nw(xi, le).'*NtetaY(xi, le))*((tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
                             (tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))) +...
             (Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*((-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj)) -...
                                  (-tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj))))) +...
           ((-wp(xi, le, upj)*tetaYp(xi, le, upj))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
            ((NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))*(-wp(xi, le, upj)*NtetaY(xi, le)) +...
             (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(-NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj)))))) +...
    dm*...
    (l2^2*2*(((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))*(-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj)) +...
              (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))) +...
             ((NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
              (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc))*(-NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj)))) -...
     2*(Nv(xi, le).'*((wp(xi, le, upj)*tetaZp(xi, le, upj)*NtetaY(xi, le) + wp(xi, le, upj)*NtetaZ(xi, le)*tetaYp(xi, le, upj)) +...
              (-tetaZp(xi, le, upj)^2*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) - (NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*tetaZp(xi, le, upj)*tetaXp(xi, le, upj))*l2 +...
              (-tetaXp(xi, le, upj)^2*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*l2) +...
        vp(xi, le, upj)*(Nw(xi, le).'*(tetaZp(xi, le, upj)*NtetaY(xi, le) + NtetaZ(xi, le)*tetaYp(xi, le, upj)) +...
            l2*NtetaZ(xi, le).'*(-tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) - (NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*tetaXp(xi, le, upj)) -...
            l2*NtetaX(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)*tetaXp(xi, le, upj)))) -...
    dm*...
    (l2^2*NtetaZ(xi, le).'*(tetaYp(xi, le, upj)^2*2*NtetaZ(xi, le)) +...
     2*NtetaZ(xi, le).'*(vp(xi, le, upj)^2 + up(xi, le, upj)^2)*NtetaZ(xi, le) +...
     2*NtetaY(xi, le).'*(wp(xi, le, upj)^2 + up(xi, le, upj)^2)*NtetaY(xi, le) +...
     2*l2*((tetaYp(xi, le, upj)*NtetaZ(xi, le).')*(up(xi, le, upj)*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
                                  (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))) -...
                              wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
           ((up(xi, le, upj)*((NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) +...
                  NtetaX(xi, le).'*(NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))) +...
                 (-NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) -...
                  NtetaX(xi, le).'*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)))) -...
             wp(xi, le, upj)*NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*(tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj)) +...
            (up(xi, le, upj)*((NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).') +...
                 (NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')) -...
             wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')*(tetaYp(xi, le, upj)*NtetaZ(xi, le)))) +...
     2*l2*((-wp(xi, le, upj)*NtetaY(xi, le).')*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
           NtetaX(xi, le).'*((-wp(xi, le, upj)*NtetaY(xi, le))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)) +...
                     (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)))) +...
     l2^2*2*NtetaX(xi, le).'*((tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)) +...
                      (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc))*(-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))) -...
     2*vp(xi, le, upj)*(wp(xi, le, upj)*(NtetaZ(xi, le).'*NtetaY(xi, le) + NtetaY(xi, le).'*NtetaZ(xi, le)) +...
           l2*tetaZp(xi, le, upj)*(-NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) - NtetaX(xi, le).'*(NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))) -...
           l2*tetaXp(xi, le, upj)*NtetaX(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))));
%repositorio.T_dgCj_dujT_m2_Aux = T_dgCj_dujT_m2_Aux;
T_dgCjdupjT_m2_Aux = @(xi, le, uj, upj)(dm*...
    (l2^2*(2*(((NtetaX(xi, le).'*tetaYp(xi, le, upj) + tetaXp(xi, le, upj)*NtetaY(xi, le).')*NtetaZ(xi, le) +...
               tetaZp(xi, le, upj)*(NtetaX(xi, le).'*NtetaY(xi, le) + NtetaY(xi, le).'*NtetaX(xi, le)))) +...
           2*(2*tetaZ(xi, le, uj))*NtetaY(xi, le).'*(NtetaY(xi, le)*tetaZp(xi, le, upj) + tetaYp(xi, le, upj)*NtetaZ(xi, le))) +...
     (2*tetaZ(xi, le, uj)*((2*vp(xi, le, upj)*Nv(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')*NtetaZ(xi, le) + 2*tetaZp(xi, le, upj)*(Nv(xi, le).'*Nv(xi, le) + Nu(xi, le).'*Nu(xi, le)))) +...
     (2*tetaY(xi, le, uj)*((2*wp(xi, le, upj)*Nw(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')*NtetaY(xi, le) + 2*tetaYp(xi, le, upj)*(Nw(xi, le).'*Nw(xi, le) + Nu(xi, le).'*Nu(xi, le))))) +...
    dm*...
    (2*l2*((((NtetaY(xi, le).'*NtetaZ(xi, le))*(up(xi, le, upj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)) + wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)) +...
             (NtetaY(xi, le).'*tetaZp(xi, le, upj))*(Nu(xi, le)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)) + Nw(xi, le)*cos(tetaX(xi, le, uj) + alfaAc))) +...
            (NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((Nu(xi, le)*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
                                              (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))) +...
                                          up(xi, le, upj)*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
                                              (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)))) +...
                                         (-sin(tetaX(xi, le, uj) + alfaAc)*(Nw(xi, le)*tetaXp(xi, le, upj) + wp(xi, le, upj)*NtetaX(xi, le))))) +...
           ((Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) + Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)))*(NtetaY(xi, le)*tetaZp(xi, le, upj) + tetaYp(xi, le, upj)*NtetaZ(xi, le)) +...
            ((-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) +...
              Nu(xi, le).'*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
                    (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))))*(NtetaX(xi, le) + NtetaY(xi, le)*tetaZ(xi, le, uj)) +...
             (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) +...
                                      Nu(xi, le).'*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
                                            (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)))))))) +...
    dm*...
    (2*l2*((-Nw(xi, le).'*(NtetaY(xi, le)*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)) +...
                   tetaYp(xi, le, upj)*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc))) +...
            (Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*(cos(tetaX(xi, le, uj) + alfaAc)*(NtetaY(xi, le)*tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*NtetaX(xi, le)) +...
                                 sin(tetaX(xi, le, uj) + alfaAc)*(NtetaZ(xi, le)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*NtetaX(xi, le)))) +...
           (-(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))*(Nw(xi, le)*tetaYp(xi, le, upj) + wp(xi, le, upj)*NtetaY(xi, le)) +...
            ((NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))*(Nu(xi, le) - Nw(xi, le)*tetaY(xi, le, uj)) +...
             (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)))))) +...
    dm*...
    (l2^2*2*((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc))*(cos(tetaX(xi, le, uj) + alfaAc)*(NtetaY(xi, le)*tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*NtetaX(xi, le)) + sin(tetaX(xi, le, uj) + alfaAc)*(NtetaZ(xi, le)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*NtetaX(xi, le))) +...
             ((NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc)) +...
              (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)))) -...
     2*(Nv(xi, le).'*(((Nw(xi, le)*tetaZp(xi, le, upj) + wp(xi, le, upj)*NtetaZ(xi, le))*tetaY(xi, le, uj) + (Nw(xi, le)*tetaYp(xi, le, upj) + wp(xi, le, upj)*NtetaY(xi, le))*tetaZ(xi, le, uj)) +...
              (2*tetaZp(xi, le, upj)*NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*(NtetaZ(xi, le)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*NtetaX(xi, le)))*l2 +...
              (2*tetaXp(xi, le, upj)*NtetaX(xi, le)*cos(tetaX(xi, le, uj) + alfaAc))*l2) +...
        ((Nw(xi, le).'*(tetaZp(xi, le, upj)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) +...
          l2*NtetaZ(xi, le).'*(tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj)) +...
          l2*NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*tetaXp(xi, le, upj))*Nv(xi, le) +...
         vp(xi, le, upj)*(Nw(xi, le).'*(NtetaZ(xi, le)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le)) +...
             l2*NtetaZ(xi, le).'*(NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le)) +...
             l2*NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le))))) -...
    dm*...
    (2*l2^2*NtetaZ(xi, le).'*((NtetaX(xi, le)*tetaYp(xi, le, upj) + tetaXp(xi, le, upj)*NtetaY(xi, le)) + 2*tetaYp(xi, le, upj)*NtetaY(xi, le)*tetaZ(xi, le, uj)) +...
     2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*(2*vp(xi, le, upj)*Nv(xi, le) + 2*up(xi, le, upj)*Nu(xi, le)) +...
     2*tetaY(xi, le, uj)*NtetaY(xi, le).'*(2*wp(xi, le, upj)*Nw(xi, le) + 2*up(xi, le, upj)*Nu(xi, le)) +...
     2*l2*(((NtetaZ(xi, le).'*NtetaY(xi, le))*(up(xi, le, upj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)) + wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)) +...
            (tetaYp(xi, le, upj)*NtetaZ(xi, le).')*(Nu(xi, le)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)) + Nw(xi, le)*cos(tetaX(xi, le, uj) + alfaAc))) +...
           ((up(xi, le, upj)*((NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).') +...
                 (NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')) -...
             wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')*(NtetaX(xi, le) + NtetaY(xi, le)*tetaZ(xi, le, uj)) +...
            (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(((NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfaAc) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).') +...
                                         (NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).'))*Nu(xi, le) -...
                                     sin(tetaX(xi, le, uj) + alfaAc)*(NtetaX(xi, le).'*Nw(xi, le))))) +...
     2*l2*(((-NtetaY(xi, le).'*Nw(xi, le))*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)) +...
            (-wp(xi, le, upj)*NtetaY(xi, le).')*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc))) +...
           NtetaX(xi, le).'*((Nu(xi, le) - Nw(xi, le)*tetaY(xi, le, uj))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)) +...
                     (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) + NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc)))) +...
     l2^2*2*NtetaX(xi, le).'*((NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfaAc) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfaAc))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc)) +...
                      (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfaAc) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc))*(NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfaAc) + NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfaAc))) -...
     2*((wp(xi, le, upj)*(NtetaZ(xi, le).'*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le).') +...
         l2*tetaZp(xi, le, upj)*(NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).') +...
         l2*tetaXp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')*Nv(xi, le) +...
        vp(xi, le, upj)*((NtetaZ(xi, le).'*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le).')*Nw(xi, le) +...
            l2*(NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfaAc) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfaAc)*NtetaX(xi, le).')*NtetaZ(xi, le) +...
            l2*cos(tetaX(xi, le, uj) + alfaAc)*(NtetaX(xi, le).'*NtetaX(xi, le))))));
%repositorio.T_dgCj_dupjT_m2_Aux = T_dgCj_dupjT_m2_Aux;

%(5.3) Plataforma móvel:
T_gCj_M3_Aux = @(xi, le, uj, upj)(2*Ip3pl*...
    (NtetaY(xi, le).'*tetaZp(xi, le, upj)*(2*tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) -...
     NtetaZ(xi, le).'*tetaYp(xi, le, upj)*(tetaZ(xi, le, uj)*tetaYp(xi, le, upj))));
%27/12/2021
T_dgCjdujT_M3_Aux = @(xi, le, uj, upj)(2*Ip3pl*...
    (NtetaY(xi, le).'*tetaZp(xi, le, upj)*(2*NtetaZ(xi, le)*tetaYp(xi, le, upj)) -...
     NtetaZ(xi, le).'*tetaYp(xi, le, upj)*(NtetaZ(xi, le)*tetaYp(xi, le, upj))));
%27/12/2021
T_dgCjdupjT_M3_Aux = @(xi, le, uj, upj)(2*Ip3pl*...
    (NtetaY(xi, le).'*2*tetaZ(xi, le, uj)*(NtetaZ(xi, le)*tetaYp(xi, le, upj) + tetaZp(xi, le, upj)*NtetaY(xi, le)) -...
     NtetaZ(xi, le).'*tetaZ(xi, le, uj)*2*tetaYp(xi, le, upj)*NtetaY(xi, le)));
%27/12/2021

%%
%(5) Matriz K (para termos elementares):
%{
piR2E = @(xi, le)(2*(dNu_dx(xi, le).'*dNu_dx(xi, le)));
%ok! 01/03/2020
piR2ksG = @(xi, le)(2*...
    (dNv_dx(xi, le).'*dNv_dx(xi, le) +...
    dNw_dx(xi, le).'*dNw_dx(xi, le) +...
    NtetaY(xi, le).'*NtetaY(xi, le) +...
    NtetaZ(xi, le).'*NtetaZ(xi, le)));
%ok! 01/03/2020
piR4E = @(xi, le)((1/4)*(2*...
    (dNtetaY_dx(xi, le).'*dNtetaY_dx(xi, le) +...
    dNtetaZ_dx(xi, le).'* dNtetaZ_dx(xi, le))));
%ok! 01/03/2020
piR4ksG = @(xi, le)((1/2)*...
    (2*...
    (dNtetaX_dx(xi, le).'*dNtetaX_dx(xi, le))));
%ok! 01/03/2020
U_kj_Aux = @(xi, le)(...
    (pi*R^2*E)*piR2E(xi, le) +...
    (pi*R^2*ks*G)*piR2ksG(xi, le) +...
    (pi*R^4*E)*piR4E(xi, le) +...
    (pi*R^4*ks*G)*piR4ksG(xi, le));
%}
U_kj_Aux = @(xi, le)(pi*R^2*E*...
    2*(dNu_dx(xi, le).'*dNu_dx(xi, le)) +...
    pi*R^4*E*...
    ((1/4)*2*(dNtetaY_dx(xi, le).'*dNtetaY_dx(xi, le) + dNtetaZ_dx(xi, le).'*dNtetaZ_dx(xi, le))) +...
    pi*R^2*ks*G*...
    2*((NtetaZ(xi, le).' - dNv_dx(xi, le).')*(NtetaZ(xi, le) - dNv_dx(xi, le)) + (NtetaY(xi, le).' + dNw_dx(xi, le).')*(NtetaY(xi, le) + dNw_dx(xi, le))) +...
    pi*R^4*ks*G*...
    ((1/2)*2*...
    (dNtetaX_dx(xi, le).'*dNtetaX_dx(xi, le))));
repositorio.U_kj_Aux = U_kj_Aux;
%ok! 21/12/2021
%ad = 1;

%(6) Matriz Fu (para termos elementares):
U_piR2E_fuj = @(xi, le, uj)(3*du_dx(xi, le, uj)^2*dNu_dx(xi, le).' + du_dx(xi, le, uj)^3*dNu_dx(xi, le).' + dv_dx(xi, le, uj)^3*dNv_dx(xi, le).' + dw_dx(xi, le, uj)^3*dNw_dx(xi, le).' +...
     (dNu_dx(xi, le).'*(dv_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)^2) + (2*dv_dx(xi, le, uj)*dNv_dx(xi, le).' + 2*dw_dx(xi, le, uj)*dNw_dx(xi, le).')*du_dx(xi, le, uj)) +...
     (du_dx(xi, le, uj)*dNu_dx(xi, le).'*dv_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*dNv_dx(xi, le).'*du_dx(xi, le, uj)^2) +...
     (du_dx(xi, le, uj)*dNu_dx(xi, le).'*dw_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*dNw_dx(xi, le).'*du_dx(xi, le, uj)^2) +...
     (dv_dx(xi, le, uj)*dNv_dx(xi, le).'*dw_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*dNw_dx(xi, le).'*dv_dx(xi, le, uj)^2));
%ok! 22/12/2021
U_piR4E_fuj = @(xi, le, uj)(((tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2)/2 +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2)/2 +...
     (2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*du_dx(xi, le, uj) + dNu_dx(xi, le).'*dtetaX_dx(xi, le, uj)^2)/2 +...
     (2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*du_dx(xi, le, uj) + dNu_dx(xi, le).'*dtetaY_dx(xi, le, uj)^2)*(3/4) +...
     (2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*du_dx(xi, le, uj) + dNu_dx(xi, le).'*dtetaZ_dx(xi, le, uj)^2)*(3/4) +...
     (dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*du_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*dNu_dx(xi, le).'*dtetaX_dx(xi, le, uj)^2)/2 +...
     (dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*du_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*dNu_dx(xi, le).'*dtetaY_dx(xi, le, uj)^2)*(3/4) +...
     (dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*du_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*dNu_dx(xi, le).'*dtetaZ_dx(xi, le, uj)^2)*(3/4) +...
     (dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*dv_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*dNv_dx(xi, le).'*dtetaX_dx(xi, le, uj)^2) +...
     (dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*dv_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*dNv_dx(xi, le).'*dtetaY_dx(xi, le, uj)^2)/4 +...
     (dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*dv_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*dNv_dx(xi, le).'*dtetaZ_dx(xi, le, uj)^2)/4 +...
     (dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*dw_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*dNw_dx(xi, le).'*dtetaX_dx(xi, le, uj)^2) +...
     (dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*dw_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*dNw_dx(xi, le).'*dtetaY_dx(xi, le, uj)^2)/4 +...
     (dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*dw_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*dNw_dx(xi, le).'*dtetaZ_dx(xi, le, uj)^2)/4 -...
     (dNtetaX_dx(xi, le).'*dtetaY_dx(xi, le, uj)*dv_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*dtetaX_dx(xi, le, uj)*dv_dx(xi, le, uj) + dNv_dx(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj))/2 -...
     (dNtetaX_dx(xi, le).'*dtetaZ_dx(xi, le, uj)*dw_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*dtetaX_dx(xi, le, uj)*dw_dx(xi, le, uj) + dNw_dx(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj))/2 -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj))/2 +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj))/2 +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*du_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*dNu_dx(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)*(3/4) +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2*du_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*dNu_dx(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)*(3/4) +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*dv_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*dNv_dx(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)/4 +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2*dv_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*dNv_dx(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)/4 +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*dw_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*dNw_dx(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)/4 +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2*dw_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*dNw_dx(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)/4 -...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2*dv_dx(xi, le, uj) + 2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dv_dx(xi, le, uj) + dNv_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2)/2 +...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2*dw_dx(xi, le, uj) + 2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dw_dx(xi, le, uj) + dNw_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2)/2 +...
     (2*tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj) + 2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*du_dx(xi, le, uj) + dNu_dx(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)*(3/4) +...
     (2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj) + 2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2*du_dx(xi, le, uj) + dNu_dx(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)*(3/4) -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + 2*du_dx(xi, le, uj)*dNu_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj))*(3/4) +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + 2*du_dx(xi, le, uj)*dNu_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj))*(3/4) -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + 2*dv_dx(xi, le, uj)*dNv_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj))/4 +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + 2*dv_dx(xi, le, uj)*dNv_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj))/4 -...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + 2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dNu_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dv_dx(xi, le, uj) + dNv_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj))/2 -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + 2*dw_dx(xi, le, uj)*dNw_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj))/4 +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + 2*dw_dx(xi, le, uj)*dNw_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj))/4 +...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + 2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dNu_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dw_dx(xi, le, uj) + dNw_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj))/2 -...
     (dNtetaX_dx(xi, le).'*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dNu_dx(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*dv_dx(xi, le, uj) + dNv_dx(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj))/2 -...
     (dNtetaX_dx(xi, le).'*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dNu_dx(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dw_dx(xi, le, uj) + dNw_dx(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj))/2 -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj) + dNu_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj))*(3/2) +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj) + dNu_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj))*(3/2)));
%ok! 22/12/2021
U_piR6E_fuj = @(xi, le, uj)((dNtetaX_dx(xi, le).'*dtetaX_dx(xi, le, uj)^3)/3 + (dNtetaY_dx(xi, le).'*dtetaY_dx(xi, le, uj)^3)/8 + (dNtetaZ_dx(xi, le).'*dtetaZ_dx(xi, le, uj)^3)/8 +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^4 + 2*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2)/6 +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^4 + 2*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2)/6 +...
     (tetaY(xi, le, uj)^3*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^4 + dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^4)/8 +...
     (tetaZ(xi, le, uj)^3*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^4 + dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^4)/8 +...
     (dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*dtetaY_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*dtetaX_dx(xi, le, uj)^2)/6 +...
     (dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*dtetaZ_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*dtetaX_dx(xi, le, uj)^2)/6 +...
     (dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*dtetaZ_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*dtetaY_dx(xi, le, uj)^2)/8 +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*dtetaY_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)/8 +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*dtetaZ_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)*(3/8) +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2*dtetaY_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)*(3/8) +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2*dtetaZ_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2)/8 +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^4 + tetaZ(xi, le, uj)*NtetaZ(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^4 + 2*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)^2)/8 -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + 3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^3)/6 -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)^3 + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)^3 + 3*dtetaZ_dx(xi, le, uj)^2*dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj))/8 +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + 3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^3)/6 +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)^3 + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)^3 + 3*dtetaY_dx(xi, le, uj)^2*dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj))/8 -...
     (3*tetaY(xi, le, uj)^2*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + 3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)^3*dtetaX_dx(xi, le, uj)^3)/8 +...
     (3*tetaZ(xi, le, uj)^2*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + 3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)^3*dtetaX_dx(xi, le, uj)^3)/8 -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaY_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + 2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)^2)/8 +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)^2 + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)^2 + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)^2 + 2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj))/8 +...
     (2*tetaY(xi, le, uj)*NtetaY(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + NtetaZ(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + 3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^3)/8 -...
     (NtetaY(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + 2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + 3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^3)/8 -...
     (NtetaY(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + NtetaZ(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + 2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj))/4);
%ok! 22/12/2021
U_piR2ksG_fuj = @(xi, le, uj)((2*tetaY(xi, le, uj)*NtetaY(xi, le).'*du_dx(xi, le, uj) + tetaY(xi, le, uj)^2*dNu_dx(xi, le).')*2 +...
     (2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*du_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*dNu_dx(xi, le).')*2 +...
     (2*tetaY(xi, le, uj)*NtetaY(xi, le).'*du_dx(xi, le, uj)^2 + 2*du_dx(xi, le, uj)*dNu_dx(xi, le).'*tetaY(xi, le, uj)^2) +...
     (2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*du_dx(xi, le, uj)^2 + 2*du_dx(xi, le, uj)*dNu_dx(xi, le).'*tetaZ(xi, le, uj)^2) -...
     (NtetaZ(xi, le).'*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dNu_dx(xi, le).'*tetaZ(xi, le, uj)*dv_dx(xi, le, uj) + dNv_dx(xi, le).'*tetaZ(xi, le, uj)*du_dx(xi, le, uj))*2 +...
     (NtetaY(xi, le).'*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dNu_dx(xi, le).'*tetaY(xi, le, uj)*dw_dx(xi, le, uj) + dNw_dx(xi, le).'*tetaY(xi, le, uj)*du_dx(xi, le, uj))*2);
%ok! 22/12/2021
U_piR4ksG_fuj = @(xi, le, uj)((tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2) +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaY_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*tetaY(xi, le, uj)^2)/2 +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^2) +...
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*dtetaZ_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)^2)/2 +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaY_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)^2)/2 +...%
     (2*tetaY(xi, le, uj)^3*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^4)/2 +...
     (tetaZ(xi, le, uj)*NtetaZ(xi, le).'*dtetaZ_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le).'*tetaZ(xi, le, uj)^2)/2 +...
     (2*tetaZ(xi, le, uj)^3*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^4)/2 -...
     (NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj))/2 +...
     (NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj))/2 +...%
     (tetaY(xi, le, uj)*NtetaY(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*NtetaZ(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)^2) -...
     (3*tetaY(xi, le, uj)^2*NtetaY(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)^3*dtetaX_dx(xi, le, uj))/2 +...
     (3*tetaZ(xi, le, uj)^2*NtetaZ(xi, le).'*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaZ(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaZ(xi, le, uj)^3*dtetaX_dx(xi, le, uj))/2 +...
     (2*tetaY(xi, le, uj)*NtetaY(xi, le).'*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + NtetaZ(xi, le).'*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + dNtetaY_dx(xi, le).'*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj))/2 -...
     (NtetaY(xi, le).'*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + 2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dNtetaX_dx(xi, le).'*tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + dNtetaZ_dx(xi, le).'*tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj))/2);
%ok! 22/12/2021
U_fuj_Aux = @(xi, le, uj)(pi*R^2*E*U_piR2E_fuj(xi, le, uj) +...
    pi*R^4*E*U_piR4E_fuj(xi, le, uj) +...
    pi*R^6*E*U_piR6E_fuj(xi, le, uj) +...
    pi*R^2*ks*G*U_piR2ksG_fuj(xi, le, uj) +...
    pi*R^4*ks*G*U_piR4ksG_fuj(xi, le, uj));

%(6.1) Matriz dFuduT (para termos elementares):
U_piR2E_dfuj = @(xi, le, uj)(3*dNu_dx(xi, le).'*2*du_dx(xi, le, uj)*dNu_dx(xi, le) + dNu_dx(xi, le).'*3*du_dx(xi, le, uj)^2*dNu_dx(xi, le) + dNv_dx(xi, le).'*3*dv_dx(xi, le, uj)^2*dNv_dx(xi, le) + dNw_dx(xi, le).'*3*dw_dx(xi, le, uj)^2*dNw_dx(xi, le) +...
     (dNu_dx(xi, le).'*(2*dv_dx(xi, le, uj)*dNv_dx(xi, le) + 2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     (2*(dNv_dx(xi, le).'*dNv_dx(xi, le)) + 2*(dNw_dx(xi, le).'*dNw_dx(xi, le)))*du_dx(xi, le, uj) + (2*dv_dx(xi, le, uj)*dNv_dx(xi, le).' + 2*dw_dx(xi, le, uj)*dNw_dx(xi, le).')*dNu_dx(xi, le)) +...
     (dNu_dx(xi, le).'*(dNu_dx(xi, le)*dv_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
      dNv_dx(xi, le).'*(dNv_dx(xi, le)*du_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le))) +...
     (dNu_dx(xi, le).'*(dNu_dx(xi, le)*dw_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
      dNw_dx(xi, le).'*(dNw_dx(xi, le)*du_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le))) +...
     (dNv_dx(xi, le).'*(dNv_dx(xi, le)*dw_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
      dNw_dx(xi, le).'*(dNw_dx(xi, le)*dv_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le))));
%ok! 22/12/2021
U_piR4E_dfuj = @(xi, le, uj)((NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)) +...
    dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)))/2 +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)) +...
    dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)))/2 +...
    (dNtetaX_dx(xi, le).'*2*(dNtetaX_dx(xi, le)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNu_dx(xi, le)) +...
    dNu_dx(xi, le).'*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le))/2 +...
    (dNtetaY_dx(xi, le).'*2*(dNtetaY_dx(xi, le)*du_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)*dNu_dx(xi, le)) +...
    dNu_dx(xi, le).'*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le))*(3/4) +...
    (dNtetaZ_dx(xi, le).'*2*(dNtetaZ_dx(xi, le)*du_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)*dNu_dx(xi, le)) +...
    dNu_dx(xi, le).'*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le))*(3/4) +...
    (dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
    dNu_dx(xi, le).'*(dNu_dx(xi, le)*dtetaX_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/2 +...
    (dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*du_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
    dNu_dx(xi, le).'*(dNu_dx(xi, le)*dtetaY_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)))*(3/4) +...
    (dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*du_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
    dNu_dx(xi, le).'*(dNu_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)))*(3/4) +...
    (dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
    dNv_dx(xi, le).'*(dNv_dx(xi, le)*dtetaX_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le))) +...
    (dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*dv_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
    dNv_dx(xi, le).'*(dNv_dx(xi, le)*dtetaY_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)))/4 +... 
    (dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*dv_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
    dNv_dx(xi, le).'*(dNv_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)))/4 +...
    (dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
    dNw_dx(xi, le).'*(dNw_dx(xi, le)*dtetaX_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le))) +...
    (dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*dw_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
    dNw_dx(xi, le).'*(dNw_dx(xi, le)*dtetaY_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)))/4 +...
    (dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*dw_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
    dNw_dx(xi, le).'*(dNw_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)))/4 -...
    (dNtetaX_dx(xi, le).'*(dNtetaY_dx(xi, le)*dv_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)*dNv_dx(xi, le)) +...
    dNtetaY_dx(xi, le).'*(dNtetaX_dx(xi, le)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNv_dx(xi, le)) +...
    dNv_dx(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)))/2 -...
    (dNtetaX_dx(xi, le).'*(dNtetaZ_dx(xi, le)*dw_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)*dNw_dx(xi, le)) +...
    dNtetaZ_dx(xi, le).'*(dNtetaX_dx(xi, le)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNw_dx(xi, le)) +...
    dNw_dx(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)))/2 -...
    (NtetaY(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
    dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
    dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj) + tetaY(xi, le, uj)*dNtetaX_dx(xi, le)))/2 +...
    (NtetaZ(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
    dNtetaX_dx(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)) +...
    dNtetaY_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)))/2 +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*(dNu_dx(xi, le)*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))*(3/4) +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*(dNu_dx(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + du_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))*(3/4) +... 
     (NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*dv_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dv_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNv_dx(xi, le).'*(dNv_dx(xi, le)*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/4 +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*dv_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dv_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNv_dx(xi, le).'*(dNv_dx(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + dv_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/4 +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*dw_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dw_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNw_dx(xi, le).'*(dNw_dx(xi, le)*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/4 +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*dw_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dw_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNw_dx(xi, le).'*(dNw_dx(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + dw_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/4 -...
    (NtetaZ(xi, le).'*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*dNv_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*NtetaZ(xi, le)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)*dNv_dx(xi, le)) +...
     dNv_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/2 +...
     (NtetaY(xi, le).'*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*dNw_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*NtetaY(xi, le)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)*dNw_dx(xi, le)) +...
     dNw_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/2 +...
    (NtetaY(xi, le).'*2*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj) + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNu_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)^2*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*(2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaY(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))*(3/4) +...
    (NtetaZ(xi, le).'*2*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNu_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*(2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))*(3/4) -...
    (NtetaY(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dNtetaZ_dx(xi, le)*du_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*2*(dNu_dx(xi, le)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + du_dx(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + du_dx(xi, le, uj)*tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + du_dx(xi, le, uj)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)))*(3/4) +...
     (NtetaZ(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)*du_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)*du_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*2*(dNu_dx(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + du_dx(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + du_dx(xi, le, uj)*tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + du_dx(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)))*(3/4) -...
    (NtetaY(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*dtetaZ_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dNtetaZ_dx(xi, le)*dv_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dv_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNv_dx(xi, le).'*2*(dNv_dx(xi, le)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dv_dx(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dv_dx(xi, le, uj)*tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dv_dx(xi, le, uj)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)))/4 +...
    (NtetaZ(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)*dv_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)*dv_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dv_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dv_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*2*dv_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNv_dx(xi, le).'*2*(dNv_dx(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dv_dx(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dv_dx(xi, le, uj)*tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + dv_dx(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)))/4 -...
    (NtetaZ(xi, le).'*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*dNu_dx(xi, le)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*NtetaZ(xi, le)*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)*dNu_dx(xi, le)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)*du_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNu_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*dv_dx(xi, le, uj) + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dv_dx(xi, le, uj) + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNv_dx(xi, le)) +...
     dNv_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNu_dx(xi, le)))/2 -...
    (NtetaY(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*dtetaZ_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dNtetaZ_dx(xi, le)*dw_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dw_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNw_dx(xi, le).'*2*(dNw_dx(xi, le)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dw_dx(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dw_dx(xi, le, uj)*tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dw_dx(xi, le, uj)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)))/4 +...
    (NtetaZ(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)*dw_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)*dw_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dw_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dw_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*2*dw_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNw_dx(xi, le).'*2*(dNw_dx(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dw_dx(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dw_dx(xi, le, uj)*tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + dw_dx(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)))/4 +...
     (NtetaY(xi, le).'*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*dNu_dx(xi, le)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*NtetaY(xi, le)*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)*dNu_dx(xi, le)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)*du_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNu_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*dw_dx(xi, le, uj) + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dw_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNw_dx(xi, le)) +...
     dNw_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*du_dx(xi, le, uj) + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNu_dx(xi, le)))/2 -...
    (dNtetaX_dx(xi, le).'*(dNtetaY_dx(xi, le)*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)*dNu_dx(xi, le)*dv_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(dNtetaX_dx(xi, le)*du_dx(xi, le, uj)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNu_dx(xi, le)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNu_dx(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)*dv_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNv_dx(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*dNu_dx(xi, le)))/2 -...
    (dNtetaX_dx(xi, le).'*(dNtetaZ_dx(xi, le)*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)*dNu_dx(xi, le)*dw_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(dNtetaX_dx(xi, le)*du_dx(xi, le, uj)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNu_dx(xi, le)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNu_dx(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)*dw_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNw_dx(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dNu_dx(xi, le)))/2 -...
    (NtetaY(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*dtetaZ_dx(xi, le, uj)*du_dx(xi, le, uj) + tetaY(xi, le, uj)*dNtetaZ_dx(xi, le)*du_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj) + tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)))*(3/2) +...
    (NtetaZ(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)*du_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)))*(3/2));
%ok! 22/12/2021
U_piR6E_dfuj = @(xi, le, uj)((dNtetaX_dx(xi, le).'*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le))/3 +...
    (dNtetaY_dx(xi, le).'*3*dtetaY_dx(xi, le, uj)^2*dNtetaY_dx(xi, le))/8 +...
    (dNtetaZ_dx(xi, le).'*3*dtetaZ_dx(xi, le, uj)^2*dNtetaZ_dx(xi, le))/8 +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^4 + tetaY(xi, le, uj)*4*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)^3*2*tetaY(xi, le, uj)*NtetaY(xi, le)))/6 +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^4 + tetaZ(xi, le, uj)*4*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)^3*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)))/6 +...
    (NtetaY(xi, le).'*(3*tetaY(xi, le, uj)^2*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^4 + tetaY(xi, le, uj)^3*4*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^4 + dtetaX_dx(xi, le, uj)^3*4*tetaY(xi, le, uj)^3*NtetaY(xi, le)))/8 +...
    (NtetaZ(xi, le).'*(3*tetaZ(xi, le, uj)^2*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^4 + tetaZ(xi, le, uj)^3*4*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^4 + dtetaX_dx(xi, le, uj)^3*4*tetaZ(xi, le, uj)^3*NtetaZ(xi, le)))/8 +...
    (dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*dtetaX_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/6 +...
    (dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*dtetaX_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/6 +...
    (dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*dtetaY_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)))/8 +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*dtetaY_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaY_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/8 +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))*(3/8) +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2*dtetaY_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))*(3/8) +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)))/8 +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^4 + tetaY(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^4 + tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*4*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le)) +...
     NtetaZ(xi, le).'*(NtetaZ(xi, le)*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^4 + tetaZ(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^4 + tetaZ(xi, le, uj)*tetaY(xi, le, uj)^2*4*dtetaX_dx(xi, le, uj)^3*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)^3*2*tetaY(xi, le, uj)*NtetaY(xi, le)*tetaZ(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)^3*tetaY(xi, le, uj)^2*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)))/8 -...
    (NtetaY(xi, le).'*(3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^3*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*3*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*tetaY(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*NtetaY(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*tetaY(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^3 + tetaY(xi, le, uj)*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)))/6 -...
    (NtetaY(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)^3 + dtetaX_dx(xi, le, uj)*3*dtetaZ_dx(xi, le, uj)^2*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*dtetaZ_dx(xi, le, uj)^3 + tetaY(xi, le, uj)*3*dtetaZ_dx(xi, le, uj)^2*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*3*(2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)^2*NtetaY(xi, le)*dtetaX_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)^2*tetaY(xi, le, uj)*dNtetaX_dx(xi, le)))/8 +...
    (NtetaZ(xi, le).'*(3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^3*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*3*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*NtetaZ(xi, le)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^3 + tetaZ(xi, le, uj)*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)))/6 +...
    (NtetaZ(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)^3 + dtetaX_dx(xi, le, uj)*3*dtetaY_dx(xi, le, uj)^2*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)^3 + tetaZ(xi, le, uj)*3*dtetaY_dx(xi, le, uj)^2*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*3*(2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)^2*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)^2*tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)))/8 -...
    (NtetaY(xi, le).'*3*(2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)^2*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^3*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*3*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*3*tetaY(xi, le, uj)^2*NtetaY(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*tetaY(xi, le, uj)^3*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(3*tetaY(xi, le, uj)^2*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^3 + tetaY(xi, le, uj)^3*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)))/8 +...
    (NtetaZ(xi, le).'*3*(2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^3*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*3*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*3*tetaZ(xi, le, uj)^2*NtetaZ(xi, le)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*tetaZ(xi, le, uj)^3*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(3*tetaZ(xi, le, uj)^2*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^3 + tetaZ(xi, le, uj)^3*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)))/8 -...
    (NtetaY(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)^2*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*dtetaY_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaY_dx(xi, le, uj)^2*dNtetaZ_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*2*(dNtetaY_dx(xi, le)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)*tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaY_dx(xi, le, uj)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)))/8 +...
    (NtetaZ(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*2*(dNtetaZ_dx(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)*tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + dtetaZ_dx(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)))/8 +...
    (NtetaY(xi, le).'*2*(NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^3*dNtetaY_dx(xi, le)) +...
     NtetaZ(xi, le).'*(2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^3*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)^2*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^3*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*3*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*2*tetaY(xi, le, uj)*NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*tetaY(xi, le, uj)^2*NtetaZ(xi, le)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(2*tetaY(xi, le, uj)*NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^3 + tetaY(xi, le, uj)^2*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^3 + tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)))/8 -...
    (NtetaY(xi, le).'*(2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^3*dNtetaZ_dx(xi, le)) +...
     NtetaZ(xi, le).'*2*(NtetaZ(xi, le)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^3*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*tetaY(xi, le, uj)*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^3*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*3*(2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*NtetaY(xi, le)*tetaZ(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*tetaY(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)^2*tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^3 + tetaY(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^3 + tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*3*dtetaX_dx(xi, le, uj)^2*dNtetaX_dx(xi, le)))/8 -...
    (NtetaY(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNtetaY_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     NtetaZ(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNtetaY_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*2*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)*NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)^2*dNtetaY_dx(xi, le)))/4);
%ok! 22/12/2021
U_piR2ksG_dfuj = @(xi, le, uj)((NtetaY(xi, le).'*2*(NtetaY(xi, le)*du_dx(xi, le, uj) + tetaY(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*2*tetaY(xi, le, uj)*NtetaY(xi, le))*2 +...
    (NtetaZ(xi, le).'*2*(NtetaZ(xi, le)*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*2*tetaZ(xi, le, uj)*NtetaZ(xi, le))*2 +...
    (NtetaY(xi, le).'*2*(NtetaY(xi, le)*du_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*2*(dNu_dx(xi, le)*tetaY(xi, le, uj)^2 + du_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le))) +...
    (NtetaZ(xi, le).'*2*(NtetaZ(xi, le)*du_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*du_dx(xi, le, uj)*dNu_dx(xi, le)) +...
     dNu_dx(xi, le).'*2*(dNu_dx(xi, le)*tetaZ(xi, le, uj)^2 + du_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le))) -...
    (NtetaZ(xi, le).'*(dNu_dx(xi, le)*dv_dx(xi, le, uj) + du_dx(xi, le, uj)*dNv_dx(xi, le)) +...
     dNu_dx(xi, le).'*(NtetaZ(xi, le)*dv_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNv_dx(xi, le)) +...
     dNv_dx(xi, le).'*(NtetaZ(xi, le)*du_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNu_dx(xi, le)))*2 +...
    (NtetaY(xi, le).'*(dNu_dx(xi, le)*dw_dx(xi, le, uj) + du_dx(xi, le, uj)*dNw_dx(xi, le)) +...
     dNu_dx(xi, le).'*(NtetaY(xi, le)*dw_dx(xi, le, uj) + tetaY(xi, le, uj)*dNw_dx(xi, le)) +...
     dNw_dx(xi, le).'*(NtetaY(xi, le)*du_dx(xi, le, uj) + tetaY(xi, le, uj)*dNu_dx(xi, le)))*2);
%ok! 22/12/2021
U_piR4ksG_dfuj = @(xi, le, uj)((NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le))) +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaY_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*tetaY(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)))/2 +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le))) +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*dtetaZ_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*tetaY(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)))/2 +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaY_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(dNtetaY_dx(xi, le)*tetaZ(xi, le, uj)^2 + dtetaY_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)))/2 +...
    (NtetaY(xi, le).'*2*(3*tetaY(xi, le, uj)^2*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaY(xi, le, uj)^3*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^4 + dtetaX_dx(xi, le, uj)*4*tetaY(xi, le, uj)^3*NtetaY(xi, le)))/2 +...
    (NtetaZ(xi, le).'*(NtetaZ(xi, le)*dtetaZ_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*dtetaZ_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(dNtetaZ_dx(xi, le)*tetaZ(xi, le, uj)^2 + dtetaZ_dx(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)))/2 +...
    (NtetaZ(xi, le).'*2*(3*tetaZ(xi, le, uj)^2*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)^3*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaZ(xi, le, uj)^4 + dtetaX_dx(xi, le, uj)*4*tetaZ(xi, le, uj)^3*NtetaZ(xi, le)))/2 -...
    (NtetaY(xi, le).'*(dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*dtetaX_dx(xi, le, uj) + tetaY(xi, le, uj)*dNtetaX_dx(xi, le)))/2 +...
    (NtetaZ(xi, le).'*(dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaZ(xi, le)*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(NtetaZ(xi, le)*dtetaX_dx(xi, le, uj) + tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)))/2 +...
    (NtetaY(xi, le).'*(NtetaY(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)) +...
     NtetaZ(xi, le).'*(NtetaZ(xi, le)*tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)^2 + tetaZ(xi, le, uj)*tetaY(xi, le, uj)^2*2*dtetaX_dx(xi, le, uj)*dNtetaX_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(dNtetaX_dx(xi, le)*tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*2*tetaY(xi, le, uj)*NtetaY(xi, le)*tetaZ(xi, le, uj)^2 + dtetaX_dx(xi, le, uj)*tetaY(xi, le, uj)^2*2*tetaZ(xi, le, uj)*NtetaZ(xi, le))) -...
    (NtetaY(xi, le).'*3*(2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(3*tetaY(xi, le, uj)^2*NtetaY(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)^3*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(3*tetaY(xi, le, uj)^2*NtetaY(xi, le)*dtetaX_dx(xi, le, uj) + tetaY(xi, le, uj)^3*dNtetaX_dx(xi, le)))/2 +...
    (NtetaZ(xi, le).'*3*(2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(3*tetaZ(xi, le, uj)^2*NtetaZ(xi, le)*dtetaY_dx(xi, le, uj) + tetaZ(xi, le, uj)^3*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(3*tetaZ(xi, le, uj)^2*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj) + tetaZ(xi, le, uj)^3*dNtetaX_dx(xi, le)))/2 +...
    (NtetaY(xi, le).'*2*(NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     NtetaZ(xi, le).'*(2*tetaY(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)^2*dtetaX_dx(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(2*tetaY(xi, le, uj)*NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)^2*NtetaZ(xi, le)*dtetaY_dx(xi, le, uj) + tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*dNtetaY_dx(xi, le)) +...
     dNtetaY_dx(xi, le).'*(2*tetaY(xi, le, uj)*NtetaY(xi, le)*tetaZ(xi, le, uj)*dtetaX_dx(xi, le, uj) + tetaY(xi, le, uj)^2*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj) + tetaY(xi, le, uj)^2*tetaZ(xi, le, uj)*dNtetaX_dx(xi, le)))/2 -...
    (NtetaY(xi, le).'*(2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     NtetaZ(xi, le).'*2*(NtetaZ(xi, le)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le)*dtetaX_dx(xi, le, uj)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*tetaY(xi, le, uj)*dNtetaX_dx(xi, le)*dtetaZ_dx(xi, le, uj) + tetaZ(xi, le, uj)*tetaY(xi, le, uj)*dtetaX_dx(xi, le, uj)*dNtetaZ_dx(xi, le)) +...
     dNtetaX_dx(xi, le).'*(NtetaY(xi, le)*tetaZ(xi, le, uj)^2*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaZ_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*dNtetaZ_dx(xi, le)) +...
     dNtetaZ_dx(xi, le).'*(NtetaY(xi, le)*tetaZ(xi, le, uj)^2*dtetaX_dx(xi, le, uj) + tetaY(xi, le, uj)*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)*dtetaX_dx(xi, le, uj) + tetaY(xi, le, uj)*tetaZ(xi, le, uj)^2*dNtetaX_dx(xi, le)))/2);
%ok! 22/12/2021
U_dfujdujT_Aux = @(xi, le, uj)(pi*R^2*E*U_piR2E_dfuj(xi, le, uj) +...
    pi*R^4*E*U_piR4E_dfuj(xi, le, uj) +...
    pi*R^6*E*U_piR6E_dfuj(xi, le, uj) +...
    pi*R^2*ks*G*U_piR2ksG_dfuj(xi, le, uj) +...
    pi*R^4*ks*G*U_piR4ksG_dfuj(xi, le, uj));

%(7) Matriz W (força peso, para termos elementares e concentrados):
wj_Aux = @(xi, le, rhoA_ou_M)(2*rhoA_ou_M*g*Nu(xi, le).');
%repositorio.wj_Aux = wj_Aux;
%ok! 01/03/2020

%(23) Matriz W referente à massa m2 (força peso):
wj_m2_Aux = @(psi, le)(2*dm*g*Nu(psi, le).');
%repositorio.wj_m2_Aux = wj_m2_Aux;

%(7) Matriz C (amortecimento, para termos elementares e concentrados):
cjbU_Aux = @(xi, le, bU)(bU*(Nu(xi, le).'*Nu(xi, le)));
cjbTx_Aux = @(xi, le, bTx)(bTx*(NtetaX(xi, le).'*NtetaX(xi, le)));
cjbV_Aux = @(xi, le, bV)(bV*(Nv(xi, le).'*Nv(xi, le)));
cjbW_Aux = @(xi, le, bW)(bW*(Nw(xi, le).'*Nw(xi, le)));
%repositorio.wj_Aux = wj_Aux;
%ok! 06/07/2022

%%
%Funções auxiliares, referentes a termos concentrados:
ad_1 = @(nElem)(double(nElem == 1));
repositorio.ad_1 = ad_1;
ad_N0 = @(nElem)(double(nElem == N0));
repositorio.ad_N0 = ad_N0;
ad_N0M1 = @(nElem)(double(nElem == N0+1));
repositorio.ad_N0M1 = ad_N0M1;
ad_N1 = @(nElem)(double(nElem == N1));
repositorio.ad_N1 = ad_N1;
ad_N1M1 = @(nElem)(double(nElem == N1+1));
repositorio.ad_N1M1 = ad_N1M1;
ad_N = @(nElem)(double(nElem == N));
repositorio.ad_N = ad_N;

%%
%Termos concentrados:
%(1) Matriz M1j_c:
%
%Inclusão da inércia interna do motor de potência:
T_M1j_c = @(nElem, le)(ad_1(nElem)*...
    (1/2)*(T_m1j_M3_Aux(0, le) + T_m1j_Motor_Aux(0, le)) +...
    ad_N1(nElem)*...
    (1/2)*(T_m1j_Aux(1, le, M2b, Ip2b, I2) + T_m1j_m2_Aux(1, le)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_m1j_Aux(0, le, M2b, Ip2b, I2) + T_m1j_m2_Aux(0, le)) +...
    ad_N(nElem)*...
    (1/2)*T_m1j_Aux(1, le, M1, Ip1, I1));
repositorio.T_M1j_c = T_M1j_c;
%}
%{
T_M1j_c = @(nElem, le)(ad_1(nElem)*...
    (1/2)*(T_m1j_M3_Aux(0, le)) +...
    ad_N1(nElem)*...
    (1/2)*(T_m1j_Aux(1, le, M2b, Ip2, I2) + T_m1j_m2_Aux(1, le)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_m1j_Aux(0, le, M2b, Ip2, I2) + T_m1j_m2_Aux(0, le)) +...
    ad_N(nElem)*...
    (1/2)*T_m1j_Aux(1, le, M1, Ip1, I1));
repositorio.T_M1j_c = T_M1j_c;
%}
%(2) Matriz M2j_c:
T_M2j_c = @(nElem, le, uj)(ad_1(nElem)*...
    (1/2)*(T_m2j_M3_Aux(0, le, uj)) +...
    ad_N1(nElem)*...
    (1/2)*(T_m2j_Aux(1, le, uj, Ip2b) + T_m2j_m2_Aux(1, le, uj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_m2j_Aux(0, le, uj, Ip2b) + T_m2j_m2_Aux(0, le, uj)) +...
    ad_N(nElem)*...
    (1/2)*T_m2j_Aux(1, le, uj, Ip1));
repositorio.T_M2j_c = T_M2j_c;

%(2) Matriz M2juppj_c:
T_M2juppj_c = @(nElem, le, uj, uppj)(ad_1(nElem)*...
    (1/2)*T_m2uppj_M3_Aux(0, le, uj, uppj) +...
    ad_N1(nElem)*...
    (1/2)*(T_m2juppj_Aux(1, le, uj, uppj, Ip2b) + T_m2juppj_m2_Aux(1, le, uj, uppj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_m2juppj_Aux(0, le, uj, uppj, Ip2b) + T_m2juppj_m2_Aux(0, le, uj, uppj)) +...
    ad_N(nElem)*...
    (1/2)*T_m2juppj_Aux(1, le, uj, uppj, Ip1));
repositorio.T_M2juppj_c = T_M2juppj_c;

%(3) Matriz dM2juppj_dujT_c:
T_dM2juppjdujT_c = @(nElem, le, uj, uppj)(ad_1(nElem)*...
    (1/2)*T_dm2juppjdujT_M3_Aux(0, le, uj, uppj) +...
    ad_N1(nElem)*...
    (1/2)*(T_dm2juppj_dujT_Aux(1, le, uj, uppj, Ip2b) + T_dm2juppjdujT_m2_Aux(1, le, uj, uppj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_dm2juppj_dujT_Aux(0, le, uj, uppj, Ip2b) + T_dm2juppjdujT_m2_Aux(0, le, uj, uppj)) +...
    ad_N(nElem)*...
    (1/2)*T_dm2juppj_dujT_Aux(1, le, uj, uppj, Ip1));
repositorio.T_dM2juppj_dujT_c = T_dM2juppjdujT_c;

%(7) Matriz Gj_c:
T_Gj_c = @(nElem, le, uj, upj)(ad_1(nElem)*...
    (1/2)*T_gCj_M3_Aux(0, le, uj, upj) +...
    ad_N1(nElem)*...
    (1/2)*(T_gCj_Aux(1, le, uj, upj, Ip2b) + T_gCj_m2_Aux(1, le, uj, upj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_gCj_Aux(0, le, uj, upj, Ip2b) + T_gCj_m2_Aux(0, le, uj, upj)) +...
    ad_N(nElem)*...
    (1/2)*T_gCj_Aux(1, le, uj, upj, Ip1));
repositorio.T_Gj_c = T_Gj_c;

%(8) Matriz dGj_dujT_c:
T_dGjdujT_c = @(nElem, le, uj, upj)(ad_1(nElem)*...
    (1/2)*T_dgCjdujT_M3_Aux(0, le, uj, upj) +...
    ad_N1(nElem)*...
    (1/2)*(T_dgCjdujT_Aux(1, le, upj, Ip2b) + T_dgCjdujT_m2_Aux(1, le, uj, upj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_dgCjdujT_Aux(0, le, upj, Ip2b) + T_dgCjdujT_m2_Aux(0, le, uj, upj)) +...
    ad_N(nElem)*...
    (1/2)*T_dgCjdujT_Aux(1, le, upj, Ip1));
repositorio.T_dGj_dujT_c = T_dGjdujT_c;

%(9) Matriz dGj_dupjT_c:
T_dGjdupjT_c = @(nElem, le, uj, upj)(ad_1(nElem)*...
    (1/2)*T_dgCjdupjT_M3_Aux(0, le, uj, upj) +...
    ad_N1(nElem)*...
    (1/2)*(T_dgCjdupjT_Aux(1, le, uj, upj, Ip2b) + T_dgCjdupjT_m2_Aux(1, le, uj, upj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_dgCjdupjT_Aux(0, le, uj, upj, Ip2b) + T_dgCjdupjT_m2_Aux(0, le, uj, upj)) +...
    ad_N(nElem)*...
    (1/2)*T_dgCjdupjT_Aux(1, le, uj, upj, Ip1));
repositorio.T_dGj_dupjT_c = T_dGjdupjT_c;

%(10) Matriz Wj_c (força peso, em termos concentrados):
Wj_c = @(nElem, le)(ad_1(nElem)*...
    (1/2)*wj_Aux(0, le, M3pl) +...
    ad_N1(nElem)*...
    (1/2)*(wj_Aux(1, le, M2b) + wj_m2_Aux(1, le)) +...
    ad_N1M1(nElem)*...
    (1/2)*(wj_Aux(0, le, M2b) + wj_m2_Aux(0, le)) +...
    ad_N(nElem)*...
    (1/2)*wj_Aux(1, le, M1));
repositorio.Wj_c = Wj_c;

%{
cjbU_Aux = @(xi, le, bU)(bU*Nu(xi, le).');
cjbTx_Aux = @(xi, le, bTx)(bTx*NtetaX(xi, le).');
cjbV_Aux = @(xi, le, bV)(bV*Nv(xi, le).');
cjbW_Aux = @(xi, le, bW)(bW*Nw(xi, le).');
%}

%(10) Matriz Wj_c (força peso, em termos concentrados):
%{
Cj_c = @(nElem, le)(...
    ad_1(nElem)*...
    (cjbU_Aux(0, le, bUl) + cjbTx_Aux(0, le, bTx1)) +...
    ad_N0(nElem)*...
    (cjbV_Aux(1, le, bV) + cjbW_Aux(1, le, bW)) +...
    ad_N0M1(nElem)*...
    (cjbV_Aux(0, le, bV) + cjbW_Aux(0, le, bW)) +...
    ad_N(nElem)*...
    (cjbU_Aux(1, le, bUNM1) + cjbTx_Aux(1, le, bTxNM1)));
%}
Cj_c = @(nElem, le)(...
    ad_1(nElem)*...
    (cjbU_Aux(0, le, bU1) + cjbTx_Aux(0, le, bTx1)) +...
    ad_N(nElem)*...
    (cjbU_Aux(1, le, bUNM1) + cjbTx_Aux(1, le, bTxNM1)));
repositorio.Cj_c = Cj_c;

%%
%Funções referentes a matrizes elementares:
%(1) Matriz M1j_elem:
fAux = @(xi, le)(T_m1j_Aux(xi, le, rho*A, rho*Jp, rho*J));
T_M1j_elem = @(le)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le)), 0, 1, nGL ));
repositorio.T_M1j_elem = T_M1j_elem;

%(2) Matriz M2j_elem:
fAux = @(xi, le, uj)(T_m2j_Aux(xi, le, uj, rho*Jp));
T_M2j_elem = @(le, uj)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le, uj)), 0, 1, nGL ));
repositorio.T_M2j_elem = T_M2j_elem;

%(2) Matriz M2juppj_elem:
fAux = @(xi, le, uj, uppj)(T_m2juppj_Aux(xi, le, uj, uppj, rho*Jp));
T_M2juppj_elem = @(le, uj, uppj)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le, uj, uppj)), 0, 1, nGL ));
repositorio.T_M2juppj_elem = T_M2juppj_elem;

%(3) Matriz dM2juppj_dujT_elem:
fAux = @(xi, le, uj, uppj)(T_dm2juppj_dujT_Aux(xi, le, uj, uppj, rho*Jp));
T_dM2juppjdujT_elem = @(le, uj, uppj)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le, uj, uppj)), 0, 1, nGL ));
repositorio.T_dM2juppj_dujT_elem = T_dM2juppjdujT_elem;

%(7) Matriz Gj_elem:
fAux = @(xi, le, uj, upj)(T_gCj_Aux(xi, le, uj, upj, rho*Jp));
T_Gj_elem = @(le, uj, upj)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le, uj, upj)), 0, 1, nGL ));
repositorio.T_Gj_elem = T_Gj_elem;

%(8) Matriz dGj_dujT_elem:
fAux = @(xi, le, upj)(T_dgCjdujT_Aux(xi, le, upj, rho*Jp));
T_dGjdujT_elem = @(le, upj)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le, upj)), 0, 1, nGL ));
repositorio.T_dGj_dujT_elem = T_dGjdujT_elem;

%(9) Matriz dGj_dupjT_elem:
fAux = @(xi, le, uj, upj)(T_dgCjdupjT_Aux(xi, le, uj, upj, rho*Jp));
T_dGjdupjT_elem = @(le, uj, upj)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le, uj, upj)), 0, 1, nGL ));
repositorio.T_dGj_dupjT_elem = T_dGjdupjT_elem;

%(10) Matriz Kj_elem (para termos elementares):
fAux = @(xi, le)(U_kj_Aux(xi, le));
U_Kj_elem = @(le)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le)), 0, 1, nGL ));
repositorio.U_Kj_elem = U_Kj_elem;

%(11) Matriz Fuj_elem (para termos elementares):
fAux = @(xi, le, uj)(U_fuj_Aux(xi, le, uj));
U_Fuj_elem = @(le, uj)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le, uj)), 0, 1, nGL ));
repositorio.U_Fuj_elem = U_Fuj_elem;

%(12) Matriz dFuj_dujT_elem (para termos elementares):
fAux = @(xi, le, uj)(U_dfujdujT_Aux(xi, le, uj));
U_dFujdujT_elem = @(le, uj)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le, uj)), 0, 1, nGL ));
repositorio.U_dFuj_dujT_elem = U_dFujdujT_elem;

%(13) Matriz Wj_elem (força peso, para termos elementares):
fAux = @(xi, le)(wj_Aux(xi, le, rho*A));
Wj_elem = @(le)((1/2)*le*...
    faaGaussLegendre( @(xi)(fAux(xi, le)), 0, 1, nGL ));
repositorio.Wj_elem = Wj_elem;

%%
%Matrizes de inércia e rigidez:
[ M1gEF, KgEF ] = fGm1k( argumentos, repositorio );

%%
%Ponto de equilíbrio no modelo em EF:
N = argumentos.N;
Ng = 6*(N + 1);
IdNg = eye(Ng);
SigmaRMaU = zeros(Ng,N+1);
for i = 1:N+1
    nSigRMaU = 6*i-5;
    SigmaRMaU(:,i) = IdNg(:,nSigRMaU);
end
%Atribuição de u1Eq = 0:
SigmaRMauEq = SigmaRMaU(:,2:N+1);

KrM1 = (SigmaRMauEq.'*KgEF)*SigmaRMauEq;
fcEFeq = @(u)(fGfcEFeq( u, argumentos, repositorio ));
fcEFrM1eq = @(urMa)(SigmaRMauEq.'*fcEFeq(SigmaRMauEq*urMa));
%Derivadas:
dfcEFeqduT = @(u)(fGdfcEFduEq( u, argumentos, repositorio ));
dfcEFrM1eqdurM1T = @(urM1)(SigmaRMauEq.'*dfcEFeqduT(SigmaRMauEq*urM1)*SigmaRMauEq);

disp('busca ponto de equilibrio')
%Ponto de equilíbrio:
Eta0 = SigmaRMauEq.'*zeros(Ng,1);
Psi = @(eta)(KrM1*eta + fcEFrM1eq(eta));
dPsidEtaT = @(eta)(KrM1 + dfcEFrM1eqdurM1T(eta));
[ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
uTilrM1eq = EtaConv;
uTilEqEF = SigmaRMauEq*uTilrM1eq;
argumentos.uTilEqEF = uTilEqEF;
disp('feito!')

%%
%Atrito

%(3) Modelo de atrito contínuo na origem
mic0 = argumentos.mic0;
w0 = argumentos.w0;
lambdaMi = argumentos.lambdami;
nMi = argumentos.nmi;
%(3.a) Funções adimensionais e auxiliares:
fiAdm = @(nu)(tanh(pi*nu) +...
    lambdaMi*(2*nMi/(2*nMi - 1))./(1 + (1/(2*nMi - 1))*nu.^(2*nMi)));
%{
nu = 0;
ad = fiAdm(nu);
%ok!
%}
dfiAdmdnu = @(nu)(-pi*(tanh(pi*nu)^2 - 1) -...
    (4*lambdaMi*nMi^2*nu^(2*nMi - 1))/(2*nMi + nu^(2*nMi) - 1)^2);
%{
nu = 0;
ad = dfiAdmdnu(nu);
%ok!
%}
d2fiAdmdnu2 = @(nu)(2*pi^2*tanh(pi*nu)*(tanh(pi*nu)^2 - 1) +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 2))/(2*nMi + nu^(2*nMi) - 1)^3 -...
    (4*lambdaMi*nMi^2*nu^(2*nMi - 2)*(2*nMi - 1))/(2*nMi + nu^(2*nMi) - 1)^2);
%{
nu = 0;
ad = d2fiAdmdnu2(nu);
%ok!
%}
d3fiAdmdnu3 = @(nu)((16*lambdaMi*nMi^3*nu^(4*nMi - 3)*(2*nMi - 1))/(2*nMi + nu^(2*nMi) - 1)^3 -...
    4*pi^3*tanh(pi*nu)^2*(tanh(pi*nu)^2 - 1) -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^4 -...
    2*pi^3*(tanh(pi*nu)^2 - 1)^2 +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 3)*(4*nMi - 2))/(2*nMi + nu^(2*nMi) - 1)^3 -...
    (4*lambdaMi*nMi^2*nu^(2*nMi - 3)*(2*nMi - 1)*(2*nMi - 2))/(2*nMi + nu^(2*nMi) - 1)^2);
%{
nu = 0;
ad = d3fiAdmdnu3(nu);
%ok!
%}
d4fiAdmdnu4 = @(nu)(16*pi^4*tanh(pi*nu)*(tanh(pi*nu)^2 - 1)^2 +...
    8*pi^4*tanh(pi*nu)^3*(tanh(pi*nu)^2 - 1) +...
    (768*lambdaMi*nMi^5*nu^(8*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^5 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 4)*(2*nMi - 1))/(2*nMi + nu^(2*nMi) - 1)^4 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 4)*(4*nMi - 2))/(2*nMi + nu^(2*nMi) - 1)^4 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 4)*(6*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^4 +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 4)*(2*nMi - 1)*(2*nMi - 2))/(2*nMi + nu^(2*nMi) - 1)^3 +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 4)*(2*nMi - 1)*(4*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^3 +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 4)*(4*nMi - 2)*(4*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^3 -...
    (4*lambdaMi*nMi^2*nu^(2*nMi - 4)*(2*nMi - 1)*(2*nMi - 2)*(2*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^2);
%{
nu = 0;
ad = d4fiAdmdnu4(nu);
%ok!
%}
d5fiAdmdnu5 = @(nu)((768*lambdaMi*nMi^5*nu^(8*nMi - 5)*(2*nMi - 1))/(2*nMi + nu^(2*nMi) - 1)^5 -...
    16*pi^5*tanh(pi*nu)^4*(tanh(pi*nu)^2 - 1) -...
    88*pi^5*tanh(pi*nu)^2*(tanh(pi*nu)^2 - 1)^2 -...
    (7680*lambdaMi*nMi^6*nu^(10*nMi - 5))/(2*nMi + nu^(2*nMi) - 1)^6 -...
    16*pi^5*(tanh(pi*nu)^2 - 1)^3 +...
    (768*lambdaMi*nMi^5*nu^(8*nMi - 5)*(4*nMi - 2))/(2*nMi + nu^(2*nMi) - 1)^5 +...
    (768*lambdaMi*nMi^5*nu^(8*nMi - 5)*(6*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^5 +...
    (768*lambdaMi*nMi^5*nu^(8*nMi - 5)*(8*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^5 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 5)*(2*nMi - 1)*(2*nMi - 2))/(2*nMi + nu^(2*nMi) - 1)^4 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 5)*(2*nMi - 1)*(4*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^4 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 5)*(2*nMi - 1)*(6*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^4 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 5)*(4*nMi - 2)*(4*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^4 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 5)*(4*nMi - 2)*(6*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^4 -...
    (96*lambdaMi*nMi^4*nu^(6*nMi - 5)*(6*nMi - 3)*(6*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^4 +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 5)*(2*nMi - 1)*(2*nMi - 2)*(2*nMi - 3))/(2*nMi + nu^(2*nMi) - 1)^3 +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 5)*(2*nMi - 1)*(2*nMi - 2)*(4*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^3 +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 5)*(2*nMi - 1)*(4*nMi - 3)*(4*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^3 +...
    (16*lambdaMi*nMi^3*nu^(4*nMi - 5)*(4*nMi - 2)*(4*nMi - 3)*(4*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^3 -...
    (4*lambdaMi*nMi^2*nu^(2*nMi - 5)*(2*nMi - 1)*(2*nMi - 2)*(2*nMi - 3)*(2*nMi - 4))/(2*nMi + nu^(2*nMi) - 1)^2);
%{
nu = 0;
ad = d5fiAdmdnu5(nu);
%Obs.: denominador nulo na última parcela, para nu=0 e nMi=2
%}
miAtr = @(tp)(mic0*fiAdm(tp/w0));
repositorio.micAtr = miAtr;
dmiAtrdtp = @(tp)((mic0/w0)*dfiAdmdnu(tp/w0));
repositorio.dmiAtrdtp = dmiAtrdtp;
d2miAtrdtp2 = @(tp)((mic0/w0^2)*d2fiAdmdnu2(tp/w0));
d3miAtrdtp3 = @(tp)((mic0/w0^3)*d3fiAdmdnu3(tp/w0));
d4miAtrdtp4 = @(tp)((mic0/w0^4)*d4fiAdmdnu4(tp/w0));
d5miAtrdtp5 = @(tp)((mic0/w0^5)*d5fiAdmdnu5(tp/w0));
repositorio.miAtr = miAtr;
repositorio.dmiAtrdtp = dmiAtrdtp;
repositorio.d2miAtrdtp2 = d2miAtrdtp2;
repositorio.d3miAtrdtp3 = d3miAtrdtp3;
repositorio.d4miAtrdtp4 = d4miAtrdtp4;
repositorio.d5miAtrdtp5 = d5miAtrdtp5;

%(2) Modelo de atrito descontínuo na origem
%(2.a) Coeficientes de atrito:
mi0 = argumentos.mi0;
mi1 = argumentos.mi1;
mi2 = argumentos.mi2;
beta1dc = argumentos.beta1mi;
beta2dc = argumentos.beta2mi;
miDc = @(tetaXp)((mi0 +...
    mi1*(1 - 2/(1 + exp(beta1dc*abs(tetaXp)))) -...
    mi2*(1 - 2/(1 + exp(beta2dc*abs(tetaXp)))))*...
    sign(tetaXp));
repositorio.miDc = miDc;
dmidtXpDc = @(tetaXp)(mi1*...
    (2*beta1dc*exp(beta1dc*abs(tetaXp))/(1 + exp(beta1dc*abs(tetaXp)))^2) -...
    mi2*...
    (2*beta2dc*exp(beta2dc*abs(tetaXp))/(1 + exp(beta2dc*abs(tetaXp)))^2));
repositorio.dmidtXpDc = dmidtXpDc;

%Modelo de atrito:
%(1) Forças de interaçao no contato de atrito:
uNM1Eq = uTilEqEF(6*(N+1)-5);
uCef = uNM1Eq;
argumentos.uCef = uCef;

%(1.a)Condição de contato:
IipNM1 = @(uNM1)(double(uNM1 >= uCef));
repositorio.IipNM1 = IipNM1;

%(1.b)Força normal:
delta0 = argumentos.delta0;
nFnNM1 = (1/pi)*(1 - uCef/delta0)^(-1);
switch opcaoSignAt
    case 1
        FnNM1 = @(uNM1)(kNM1*(uNM1 - uCef).*double(uNM1 >= uCef));
        dFnNM1duNM1 = @(uNM1)(kNM1*double(uNM1 >= uCef));
    case 2
        FnNM1 = @(uNM1)(kNM1*(uNM1 - uCef).*(tanh(nFnNM1*pi*(uNM1/delta0 - 1)) + 1));
        dFnNM1duNM1 = @(uNM1)(kNM1.*(tanh(nFnNM1*pi*(uNM1/delta0 - 1)) + 1) +...
            kNM1*(uNM1 - uCef)*(-(pi*nFnNM1*(tanh(pi*nFnNM1*(uNM1/delta0 - 1))^2 - 1))/delta0));
    otherwise
        disp('other value')
end
repositorio.FnNM1 = FnNM1;
FnbNM1 = @(ubNM1)(FnNM1(ubNM1 + uNM1Eq));
repositorio.FnbNM1 = FnbNM1;
repositorio.dFnNM1duNM1 = dFnNM1duNM1;

%Funções de atrito no ponto de contato (nó N+1):
%Atrito descontínuo na origem:
TxNM1dc = @(nElem,tetaXpNM1,uNM1)(ad_N(nElem)*miDc(tetaXpNM1)*FnNM1(uNM1)*R_NM1);
repositorio.TxNM1nElemDc = TxNM1dc;
miNM1nElemDc = @(nElem, tetaXpNM1)(ad_N(nElem)*miDc(tetaXpNM1));
repositorio.miNM1nElemDc = miNM1nElemDc;
%Atrito contínuo na origem:
TxNM1c = @(nElem,tetaXpNM1,uNM1)(ad_N(nElem)*miAtr(tetaXpNM1)*FnNM1(uNM1)*R_NM1);
repositorio.TxNM1nElemC = TxNM1c;
miNM1nElemC = @(nElem, tetaXpNM1)(ad_N(nElem)*miAtr(tetaXpNM1));
repositorio.miNM1nElemC = miNM1nElemC;
dmidtXpNM1nElemC = @(nElem, tetaXpNM1)(ad_N(nElem)*...
    dmiAtrdtp(tetaXpNM1));
repositorio.dmidtXpNM1nElemC = dmidtXpNM1nElemC;

%%
%RUBBING
R2 = Dr2/2;
n = (1/pi)*(1 - fg/delta0)^-1;

%Órbita do rotor intermediário:
fuL = @(yL,zL)(sqrt(yL^2 + zL^2)*double(sqrt(yL^2 + zL^2)>eps) +...
    eps*double(sqrt(yL^2 + zL^2)<=eps));
repositorio.fuL = fuL;

%Força normal de contato:
switch opcao_sign
    case 1
        fFn = @(uL)(Kip*(uL - fg)*double(uL>=fg));
    case 2
        fFn = @(uL)(Kip*(uL - fg)*(tanh(n*pi*(uL/delta0 - 1)) + 1));
    otherwise
        disp('other value')
end
repositorio.fFn = fFn;

%Força tangencial de atrito no contato:
%(a) Modelo de atrito contínuo na origem:
miIp0c = argumentos.miIp0c;
v0Ipc = argumentos.v0Ipc;
nmiIpc = argumentos.nmiIpc;
lambdamiIpc = argumentos.lambdamiIpc;
fiIpc = @(nu)(tanh(pi*nu) +...
    lambdamiIpc*((2*nmiIpc)/(2*nmiIpc - 1))*...
    (1 + (1/(2*nmiIpc - 1))*nu^(2*nmiIpc))^-1);
miIpc = @(vCfi)(miIp0c*fiIpc(vCfi/v0Ipc));
repositorio.miIpc = miIpc;
fFtc = @(uL,vC)(miIpc(vC)*fFn(uL));
repositorio.fFtc = fFtc;
%(b) Modelo de atrito descontínuo na origem:
miIp0dc = argumentos.mi0ip;
miIp1dc = argumentos.mi1ip;
miIp2dc = argumentos.mi2ip;
beta1miIpdc = argumentos.beta1miIp;
beta2miIpdc = argumentos.beta2miIp;
nS = 5;
signS = @(nu)(tanh(nS*pi*nu));
dsignSdnu = @(nu)(-pi*nS*(tanh(pi*nS*nu)^2 - 1));
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
fFtdc = @(uL,vC)(miIpdc(vC)*fFn(uL));
repositorio.fFtdc = fFtdc;
repositorio.miIpdc = miIpdc;

%Velocidade parcial de contato:
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
%(a) Versor unitário na direção tangencial:
fuTilHatfiL = @(yL,zL)([-zL/fuL(yL,zL); yL/fuL(yL,zL)]);
%(b) Velocidade parcial do ponto de contato:
fvCfi = @(yL,zL,teta2p,yLp,zLp)([1, fuTilHatfiL(yL,zL).']*AvCfi(yL,zL)*[teta2p; yLp; zLp]);
repositorio.fvCfi = fvCfi;

%Forças e torque de interação no contato:
%(a) Modelo de atrito contínuo na origem:
fFyc = @(yL,zL,teta2p,yLp,zLp)(-fFn(fuL(yL,zL))*(yL/fuL(yL,zL)) +...
    fFtc(fuL(yL,zL),fvCfi(yL,zL,teta2p,yLp,zLp))*(zL/fuL(yL,zL)));
repositorio.fFyRubc = fFyc;
fFzc = @(yL,zL,teta2p,yLp,zLp)(-fFn(fuL(yL,zL))*(zL/fuL(yL,zL)) -...
    fFtc(fuL(yL,zL),fvCfi(yL,zL,teta2p,yLp,zLp))*(yL/fuL(yL,zL)));
repositorio.fFzRubc = fFzc;
fT2c = @(yL,zL,teta2p,yLp,zLp)(-miIpc(fvCfi(yL,zL,teta2p,yLp,zLp))*...
    fFn(fuL(yL,zL))*(R2 + fg));
repositorio.fT2cRub = fT2c;
%(b) Modelo de atrito descontínuo na origem:
fFydc = @(yL,zL,teta2p,yLp,zLp)(-fFn(fuL(yL,zL))*(yL/fuL(yL,zL)) +...
    fFtdc(fuL(yL,zL),fvCfi(yL,zL,teta2p,yLp,zLp))*(zL/fuL(yL,zL)));
repositorio.fFyRubdc = fFydc;
fFzdc = @(yL,zL,teta2p,yLp,zLp)(-fFn(fuL(yL,zL))*(zL/fuL(yL,zL)) -...
    fFtdc(fuL(yL,zL),fvCfi(yL,zL,teta2p,yLp,zLp))*(yL/fuL(yL,zL)));
repositorio.fFzRubdc = fFzdc;
fT2dc = @(yL,zL,teta2p,yLp,zLp)(-miIpdc(fvCfi(yL,zL,teta2p,yLp,zLp))*...
    fFn(fuL(yL,zL))*(R2 + fg));
repositorio.fT2dcRub = fT2dc;

%%
%Rubbing (DERIVADAS):
inc_dsignS = argumentos.inc_dsignS;

%Força normal de contato:
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

fduLdyL = @(yL,zL)(yL/fuL(yL,zL));
repositorio.fduLdyL = fduLdyL;
fduLdzL = @(yL,zL)(zL/fuL(yL,zL));
repositorio.fduLdzL = fduLdzL;
fdyLuLpwm1dyL = @(yL,zL)((1/fuL(yL,zL))*(1 - (yL/fuL(yL,zL))^2));
repositorio.fdyLuLpwm1dyL = fdyLuLpwm1dyL;
fdyLuLpwm1dzL = @(yL,zL)(-yL*zL/(fuL(yL,zL))^3);
repositorio.fdyLuLpwm1dzL = fdyLuLpwm1dzL;
fdzLuLpwm1dzL = @(yL,zL)((1/fuL(yL,zL))*(1 - (zL/fuL(yL,zL))^2));
repositorio.fdzLuLpwm1dzL = fdzLuLpwm1dzL;
fdzLuLpwm1dyL = @(yL,zL)(-zL*yL/(fuL(yL,zL))^3);
repositorio.fdzLuLpwm1dyL = fdzLuLpwm1dyL;

%Modelo de atrito contínuo na origem:
dfiIpcdvCfi = @(nu)(- pi*(tanh(pi*nu)^2 - 1) -...
    (4*lambdamiIpc*nmiIpc^2*nu^(2*nmiIpc - 1))/(2*nmiIpc + nu^(2*nmiIpc) - 1)^2);
dmiIpcdvCfi = @(vCfi)(miIp0c*(1/v0Ipc)*dfiIpcdvCfi(vCfi/v0Ipc));
repositorio.dmiIpcdvCfi = dmiIpcdvCfi;

%Modelo de atrito descontínuo na origem:
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
            inc_dsignS*dsignSdnu(vCfi)*...
            (miIp0dc +...
            miIp1dc*(1 - 2*(1 + exp(beta1miIpdc*abs(vCfi)))^-1) +...
            miIp2dc*(1 - 2*(1 + exp(beta2miIpdc*abs(vCfi)))^-1)));
    otherwise
        disp('other value')
end
repositorio.dmiIpdcdvCfi = dmiIpdcdvCfi;

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

fdvCfidteta2p = @(yL,zL)([1, fuTilHatfiL(yL,zL).']*AvCfi(yL,zL)*[1; 0; 0]);
repositorio.fdvCfidteta2p = fdvCfidteta2p;
fdvCfidyLp = @(yL,zL)([1, fuTilHatfiL(yL,zL).']*AvCfi(yL,zL)*[0; 1; 0]);
repositorio.fdvCfidyLp = fdvCfidyLp;
fdvCfidzLp = @(yL,zL)([1, fuTilHatfiL(yL,zL).']*AvCfi(yL,zL)*[0; 0; 1]);
repositorio.fdvCfidzLp = fdvCfidzLp;

%%
%Rigidez em torno do ponto de equilíbrio:
KgEFb = KgEF + fGdfuduT( uTilEqEF, argumentos, repositorio );
SigmaReF = [IdNg(:,1),IdNg(:,6),IdNg(:,7:6*N),IdNg(:,6*(N+1)-5),IdNg(:,6*(N+1))];
argumentos.SigmaReF = SigmaReF;
M1restr = SigmaReF.'*M1gEF*SigmaReF;
Krestr = SigmaReF.'*KgEF*SigmaReF;
KbRestr = SigmaReF.'*KgEFb*SigmaReF;

[ PtX, sigTiltX, Pu, sigTilu, Pvw, sigTilvw, Abvt, Abvl, info ] =...
    fModMMod( M1restr, KbRestr, argumentos );
%Primeiros modos torsionais:
argumentos.PtX = PtX;
%Primeiras frequências torsionais:
argumentos.sigTiltX = sigTiltX;
%Primeiros modos longitudinais:
argumentos.Pu = Pu;
%Primeiras frequências longitudinais:
argumentos.sigTilu = sigTilu;
%Modos reversos:
argumentos.Pvw = Pvw;
%Frequências dos modos reversos:
argumentos.sigTilvw = sigTilvw;
%Posição do modo PtX0 em Abvt (matriz de modos já ordenados):
jPtX0 = info.jPtX0;
argumentos.jPtX0 = jPtX0;
%Posição do modo PtX1 em Abvt:
jPtX1 = info.jPtX1;
argumentos.jPtX1 = jPtX1;
%Posição do modo PtX2 em Abvt:
jPtX2 = info.jPtX2;
argumentos.jPtX2 = jPtX2;
%Posição do modo Pu0 em Abvt:
jPu0 = info.jPu0;
argumentos.jPu0 = jPu0;
%Posição do modo Pu1 em Abvt:
jPu1 = info.jPu1;
argumentos.jPu1 = jPu1;
%Posição do modo Pu2 em Abvt:
jPu2 = info.jPu2;
argumentos.jPu2 = jPu2;
%Posição do modo Pv em Abvt:
jPv = info.jPv;
argumentos.jPv = jPv;
%Posição do modo Pw em Abvt:
jPw = info.jPw;
argumentos.jPw = jPw;
%Frequências do sistema em Hz:
wEF = sqrt(Abvl);
freqEF = wEF/(2*pi);

%%
disp('Matriz de amortecimento')
%Matriz de amortecimento:
%(a) Amortecimento estrutural:
nSig = length(sigTiltX);
wRmin = wEF(jPtX1,1); %rad/s
wRmax = max(wEF); %rad/s
qsiMax = 1; %Valor máximo dos coeficientes de amortecimento
if nSig<3
    %Matriz de amortecimento:
    %Obs: 
    %{
    (1) C = alfa*M1 + beta*K
    (2) alfa e beta foram escolhidos de modo que todos os modos fossem subamort:
       (a) alfa = 2*(wMin*wMax)/(wMin + wMax)
       (b) beta = 2/(wMin + wMax)
       (c) wMin e wMax são tais que f(wMin) = f(wMax) = 1
       (d) f(w) = (1/2)*(alfa/w + beta*w)
       (e) qsiMin: menor coeficiente de amortecimento nominal
    (3) Parâmetros para determinação de alfa e beta: qsiMin, wRmin, wRmax
       (a) wRmin: menor frequência natural do sistema
       (b) wRmax: maior frequência natural do sistema
       (c) wMin < wRmin <= w <= wRmax < wMax
    (4) Sequência de cálculo:
       (a) determinação de qsiMin, tal que 0 <= qsiMin <= qsiMinMax
       (b) deltaQsi = qsi*(1 + sqrt(1 - qsi^2))^-1
       (c) deltaQsi0 = wRmin/wRmax
       (d) deltaQsiMax = sqrt(wRmin/wRmax)
       (e) wQsi: wRmin <= wQsi <= wRmax, se 0 <= qsiMin < qsiMin0
           wQsi: wRmax*deltaQsi <= wQsi <= wRmin/deltaQsi, se qsiMin0 <= qsiMin <= qsiMinMax
    %}
    argumentos.qsiMax = qsiMax;
    deltaQsi0 = wRmin/wRmax;
    deltaQsiMax = sqrt(wRmin/wRmax);
    qsiMin0 = 2*qsiMax*(deltaQsi0 + deltaQsi0^-1)^-1;
    qsiMinMax = 2*qsiMax*(deltaQsiMax + deltaQsiMax^-1)^-1;
    sigQsi = 0.5; %Condição: 0 <= sigQsi <= 1
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
    alfaEF = 2*qsiMax*(wMin*wMax)/(wMin + wMax);
    betaEF = 2*qsiMax/(wMin + wMax);
    deltaMin = log10(wRmin) - log10(wMin);
    argumentos.deltaMin = deltaMin;
    deltaMax = log10(wMax) - log10(wRmax);
    argumentos.deltaMax = deltaMax;
    %(a) Amortecimento estrutural:
    CgEFestr = alfaEF*M1gEF + betaEF*KgEFb;
    wNatTil = freqNatKbw;
    psiAmort = (1/2)*(alfaEF./wNatTil + betaEF*wNatTil);
    argumentos.psiAmort = psiAmort;
else
    %{
    deltaMax = 0.1;
    wMax = wRmax*10^deltaMax;
    betaEF = 2*qsiMax/wMax;
    %}
    betaEF = 0;
    %betaEF = 0.5;
    %(a) Amortecimento estrutural:
    CgEFestr = betaEF*KgEFb;
    wNatTil = wEF;
    psiAmort = (1/2)*(betaEF*wNatTil);
    argumentos.psiAmort = psiAmort;
end
%(b) Amortecimento concentrado:
[ CgEFconc ] = fGcAmort( argumentos, repositorio );
%(c) Amortecimento total:
adCconc = 1; %(i) adCx == 0: amort nulo; (ii) adCx == 1: amort não nulo
adCestr = 1;
CgEF = adCestr*CgEFestr + adCconc*CgEFconc;
argumentos.CgEF = CgEF;
Crestr = SigmaReF.'*CgEF*SigmaReF;

Bef = [IdNg(:,6*1),IdNg(:,6*1-5)];
Br = SigmaReF.'*Bef;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vetor fEF(u,up,upp):
%(1) Modelo contínuo de atrito:
fcEF = @(u,up,upp)(fGfcEF( u, up, upp, argumentos, repositorio ));
repositorio.fcEF = fcEF;
%(2) Modelo descontínuo de atrito:
fdcEF = @(u,up,upp)(fGfdcEF( u, up, upp, argumentos, repositorio ));
repositorio.fdcEF = fdcEF;

%Força peso sobressalente:
[ Wef ] = fGw( argumentos, repositorio );
[ Fu ] = fGfu( uTilEqEF, argumentos, repositorio );
Wsb = Wef - KgEF*uTilEqEF - Fu;
Wsbr = SigmaReF.'*Wsb;

%Definição de variáveis e funções, em torno do ponto de equilíbrio:
ub0 = zeros(Ng,1); ubp0 = zeros(Ng,1); ubpp0 = zeros(Ng,1);
%(1) Modelo contínuo de atrito:
fcEFb = @(ub,ubp,ubpp)(KgEF*uTilEqEF + fcEF(ub + uTilEqEF,ubp,ubpp) + Wsb);
fc0 = fcEFb(ub0,ubp0,ubpp0);
fcEFb = @(ub,ubp,ubpp)(fcEFb(ub,ubp,ubpp) - fc0);
repositorio.fcEFb = fcEFb;
fcEFbr = @(ubr,ubpr,ubppr)(SigmaReF.'*fcEFb(SigmaReF*ubr,SigmaReF*ubpr,SigmaReF*ubppr));
%(2) Modelo descontínuo de atrito:
fdcEFb = @(ub,ubp,ubpp)(KgEF*uTilEqEF + fdcEF(ub + uTilEqEF,ubp,ubpp) + Wsb);
fdc0 = fdcEFb(ub0,ubp0,ubpp0);
fdcEFb = @(ub,ubp,ubpp)(fdcEFb(ub,ubp,ubpp) - fdc0);
repositorio.fdcEFb = fdcEFb;
fdcEFbr = @(ubr,ubpr,ubppr)(SigmaReF.'*fdcEFb(SigmaReF*ubr,SigmaReF*ubpr,SigmaReF*ubppr));

%Derivadas:
dudubT = IdNg;
dupdubpT = IdNg;
duppdubppT = IdNg;
%(1) Modelo contínuo de atrito:
dfcEFduT = @(u,up,upp)(fGdfcEFdu( u, up, upp, argumentos, repositorio ));
dfcEFbdubT = @(ub,ubp,ubpp)(dfcEFduT(ub + uTilEqEF,ubp,ubpp)*dudubT);
repositorio.dfcEFbdubT = dfcEFbdubT;
%{
ub0 = zeros(Ng,1);
ubp0 = zeros(Ng,1);
ubpp0 = zeros(Ng,1);
ad = dfcEFbdubT(ub0,ubp0,ubpp0);
%}
dfcEFbrdubrT = @(ubr,ubpr,ubppr)(...
    SigmaReF.'*dfcEFbdubT(SigmaReF*ubr,SigmaReF*ubpr,SigmaReF*ubppr)*SigmaReF);
dfcEFdupT = @(u,up)(fGdfcEFdup( u, up, argumentos, repositorio ));
dfcEFbdubpT = @(ub,ubp)(dfcEFdupT(ub + uTilEqEF,ubp)*dupdubpT);
repositorio.dfcEFbdubpT = dfcEFbdubpT;
dfcEFbrdubprT = @(ubr,ubpr)(...
    SigmaReF.'*dfcEFbdubpT(SigmaReF*ubr,SigmaReF*ubpr)*SigmaReF);
%(2) Modelo descontínuo de atrito:
dfdcEFduT = @(u,up,upp)(fGdfdcEFdu( u, up, upp, argumentos, repositorio ));
dfdcEFbdubT = @(ub,ubp,ubpp)(dfdcEFduT(ub + uTilEqEF,ubp,ubpp)*dudubT);
repositorio.dfdcEFbdubT = dfdcEFbdubT;
dfdcEFbrdubrT = @(ubr,ubpr,ubppr)(...
    SigmaReF.'*dfdcEFbdubT(SigmaReF*ubr,SigmaReF*ubpr,SigmaReF*ubppr)*SigmaReF);
dfdcEFdupT = @(u,up)(fGdfdcEFdup( u, up, argumentos, repositorio ));
dfdcEFbdubpT = @(ub,ubp)(dfdcEFdupT(ub + uTilEqEF,ubp)*dupdubpT);
repositorio.dfdcEFbdubpT = dfdcEFbdubpT;
dfdcEFbrdubprT = @(ubr,ubpr)(...
    SigmaReF.'*dfdcEFbdubpT(SigmaReF*ubr,SigmaReF*ubpr)*SigmaReF);
%(3) Derivada em relação à aceleração:
dfEFduppT = @(u)(fGdfEFdupp( u, argumentos, repositorio ));
dfEFbdubppT = @(ub)(dfEFduppT(ub + uTilEqEF)*duppdubppT);
repositorio.dfEFbdubppT = dfEFbdubppT;
dfEFbrdubpprT = @(ubr)(...
    SigmaReF.'*dfEFbdubppT(SigmaReF*ubr)*SigmaReF);

%%
%Vetor fEFaprox(u,up,upp):
indG = 0; 
argumentos.indG = indG;
%Obs.: 
%indG = 1: inclusão da função G em fc[d]EFaprox
%indG = 0: retirada da função G em fc[d]EFaprox
%(1) Modelo contínuo de atrito:
%Obs.: foram retirados os termos M2upp e FRubb
fcEFaprox = @(u,up)(fGfcEFaprox( u, up, argumentos, repositorio ));
repositorio.fcEFaprox = fcEFaprox;
%(2) Modelo descontínuo de atrito:
fdcEFaprox = @(u,up)(fGfdcEFaprox( u, up, argumentos, repositorio ));
repositorio.fdcEFaprox = fdcEFaprox;
%Definição de variáveis e funções, em torno do ponto de equilíbrio:
%(1) Modelo contínuo de atrito:
fcEFbAprox = @(ub,ubp)(KgEF*uTilEqEF + fcEFaprox(ub + uTilEqEF,ubp));
repositorio.fcEFbAprox = fcEFbAprox;
fcEFbAtr = @(ub,ubp)(fGfcEFbAtr( ub, ubp, argumentos, repositorio ));
repositorio.fcEFbAtr = fcEFbAtr;
%(2) Modelo descontínuo de atrito:
fdcEFbAprox = @(ub,ubp)(KgEF*uTilEqEF + fdcEFaprox(ub + uTilEqEF,ubp));
repositorio.fdcEFbAprox = fdcEFbAprox;
fdcEFbAtr = @(ub,ubp)(fGfdcEFbAtr( ub, ubp, argumentos, repositorio ));
repositorio.fdcEFbAtr = fdcEFbAtr;

%(1a) Aproximação:
dfcEFduTaprox = @(u,up)(fGdfcEFduAprox( u, up, argumentos, repositorio ));
dfcEFbdubTaprox = @(ub,ubp)(dfcEFduTaprox(ub + uTilEqEF,ubp)*dudubT);
repositorio.dfcEFbdubTaprox = dfcEFbdubTaprox;
dfcEFdupTaprox = @(u,up)(fGdfcEFdupAprox( u, up, argumentos, repositorio ));
dfcEFbdubpTaprox = @(ub,ubp)(dfcEFdupTaprox(ub + uTilEqEF,ubp)*dupdubpT);
repositorio.dfcEFbdubpTaprox = dfcEFbdubpTaprox;
dfcEFdubTAtr = @(ub,ubp)(fGdfcEFdubAtr( ub, ubp, argumentos, repositorio ));
repositorio.dfcEFdubTAtr = dfcEFdubTAtr;
dfcEFdubpTAtr = @(ub,ubp)(fGdfcEFdubpAtr( ub, ubp, argumentos, repositorio ));
repositorio.dfcEFdubpTAtr = dfcEFdubpTAtr;

%(2a) Aproximação:
dfdcEFduTaprox = @(u,up)(fGdfdcEFduAprox( u, up, argumentos, repositorio ));
dfdcEFbdubTaprox = @(ub,ubp)(dfdcEFduTaprox(ub + uTilEqEF,ubp)*dudubT);
repositorio.dfdcEFbdubTaprox = dfdcEFbdubTaprox;
dfdcEFdupTaprox = @(u,up)(fGdfdcEFdupAprox( u, up, argumentos, repositorio ));
dfdcEFbdubpTaprox = @(ub,ubp)(dfdcEFdupTaprox(ub + uTilEqEF,ubp)*dupdubpT);
repositorio.dfdcEFbdubpTaprox = dfdcEFbdubpTaprox;
dfdcEFdubTAtr = @(ub,ubp)(fGdfdcEFdubAtr( ub, ubp, argumentos, repositorio ));
repositorio.dfdcEFdubTAtr = dfdcEFdubTAtr;
dfdcEFdubpTAtr = @(ub,ubp)(fGdfdcEFdubpAtr( ub, ubp, argumentos, repositorio ));
repositorio.dfdcEFdubpTAtr = dfdcEFdubpTAtr;
%{
ad1 = fcEFb(ub0,ubp0,ubpp0);
ad2 = fdcEFb(ub0,ubp0,ubpp0);
ad3 = dfcEFbdubT(ub0,ubp0,ubpp0);
ad4 = dfdcEFbdubT(ub0,ubp0,ubpp0);
ad5 = dfdcEFbdubpT(ub0,ubp0);
%}

%%
%AJUSTE DE MODELO
disp('ajuste de modelo')

%Dinâmica de torção:
Id3 = eye(3);
argumentos.Id3 = Id3;
fI1 = @(alfa)(alfa*Ip1); 
fI2 = @(alfa)(alfa*(Ib2 + dm*l2^2)); 
fI3 = @(alfa)(alfa*I3);

fI12 = @(alfa1,alfa2)((fI1(alfa1)^-1 + (fI2(alfa2))^-1)^-1);
fI23 = @(alfa2,alfa3)(((fI2(alfa2))^-1 + fI3(alfa3)^-1)^-1);
fI123 = @(alfa1,alfa2,alfa3)(...
    (1/(fI3(alfa3)*fI1(alfa1)) + 1/((fI2(alfa2))*fI1(alfa1)) + 1/((fI2(alfa2))*fI3(alfa3)))^-1);
fsigTalfa = @(alfa1,alfa2,alfa3)(1 - fI12(alfa1,alfa2)*fI23(alfa2,alfa3)/fI123(alfa1,alfa2,alfa3));
falfa1Aj = @(alfaTil)(Id3(1,:)*alfaTil);
falfa2Aj = @(alfaTil)(Id3(2,:)*alfaTil);
falfa3Aj = @(alfaTil)(Id3(3,:)*alfaTil);
fsigTalfaTil = @(alfaTil)(fsigTalfa(falfa1Aj(alfaTil),falfa2Aj(alfaTil),falfa3Aj(alfaTil)));

sigT1 = sigTiltX(2,1);
sigT2 = sigTiltX(3,1);

%Dinâmica longitudinal:
fM1 = @(alfa)(alfa*M1); 
fM2 = @(alfa)(alfa*(Mb2 + dm)); 
fM3 = @(alfa)(alfa*M3);

fM12 = @(alfa1,alfa2)((fM1(alfa1)^-1 + (fM2(alfa2))^-1)^-1);
fM23 = @(alfa2,alfa3)(((fM2(alfa2))^-1 + fM3(alfa3)^-1)^-1);
fM123 = @(alfa1,alfa2,alfa3)(...
    (1/(fM3(alfa3)*fM1(alfa1)) + 1/((fM2(alfa2))*fM1(alfa1)) + 1/((fM2(alfa2))*fM3(alfa3)))^-1);
fsigUalfa = @(alfa1,alfa2,alfa3)(1 - fM12(alfa1,alfa2)*fM23(alfa2,alfa3)/fM123(alfa1,alfa2,alfa3));
fsigUalfaTil = @(alfaTil)(fsigUalfa(falfa1Aj(alfaTil),falfa2Aj(alfaTil),falfa3Aj(alfaTil)));

sigU1 = sigTilu(2,1);
sigU2 = sigTilu(3,1);

%(1) Ajuste 1 (apenas rigidezes)
%Obs.: lmb2 = sig2 (segunda frequência de torção no modelo EF)
a1 = 1; a2 = 1; a3 = 1;

sigAlfa = fsigTalfa(a1,a2,a3);
lmb2 = sigT2;
lmb1b = lmb2*sqrt(2)/(1 + sqrt(sigAlfa));
lmb1 = (1 - sqrt(sigAlfa))*lmb1b/sqrt(2);
K2Aj1 = (1/2)*fI12(a1,a2)*(lmb1 + lmb2);
K1Aj1 = fI123(a1,a2,a3)*(lmb1*lmb2)/K2Aj1;

sigAlfa = fsigUalfa(a1,a2,a3);
lmb2 = sigU2;
lmb1b = lmb2*sqrt(2)/(1 + sqrt(sigAlfa));
lmb1 = (1 - sqrt(sigAlfa))*lmb1b/sqrt(2);
k2Aj1 = (1/2)*fM12(a1,a2)*(lmb1 + lmb2);
k1Aj1 = fM123(a1,a2,a3)*(lmb1*lmb2)/k2Aj1;

%(2) Ajuste 2 (apenas rigidezes)
%Obs.: lmb1b ótimo
a1 = 1; a2 = 1; a3 = 1;

sigAlfa = fsigTalfa(a1,a2,a3);
lmb1b = sqrt(2)*((sigT1*(1 - sqrt(sigAlfa)) + sigT2*(1 + sqrt(sigAlfa)))/...
    ((1 - sqrt(sigAlfa))^2 + (1 + sqrt(sigAlfa))^2));
lmb1 = (1 - sqrt(sigAlfa))*lmb1b/sqrt(2);
lmb2 = (1 + sqrt(sigAlfa))*lmb1b/sqrt(2);
K2Aj2 = (1/2)*fI12(a1,a2)*(lmb1 + lmb2);
K1Aj2 = fI123(a1,a2,a3)*(lmb1*lmb2)/K2Aj2;

sigAlfa = fsigUalfa(a1,a2,a3);
lmb1b = sqrt(2)*((sigU1*(1 - sqrt(sigAlfa)) + sigU2*(1 + sqrt(sigAlfa)))/...
    ((1 - sqrt(sigAlfa))^2 + (1 + sqrt(sigAlfa))^2));
lmb1 = (1 - sqrt(sigAlfa))*lmb1b/sqrt(2);
lmb2 = (1 + sqrt(sigAlfa))*lmb1b/sqrt(2);
k2Aj2 = (1/2)*fM12(a1,a2)*(lmb1 + lmb2);
k1Aj2 = fM123(a1,a2,a3)*(lmb1*lmb2)/k2Aj2;

%(3) Ajuste 3 (apenas rigidezes)
%Obs.: lmb1 = sig1 (primeira frequência de torção no modelo EF)
a1 = 1; a2 = 1; a3 = 1;

sigAlfa = fsigTalfa(a1,a2,a3);
lmb1 = sigT1;
lmb1b = lmb1*sqrt(2)/(1 - sqrt(sigAlfa));
lmb2 = (1 + sqrt(sigAlfa))*lmb1b/sqrt(2);
K2Aj3 = (1/2)*fI12(a1,a2)*(lmb1 + lmb2);
K1Aj3 = fI123(a1,a2,a3)*(lmb1*lmb2)/K2Aj3;

sigAlfa = fsigUalfa(a1,a2,a3);
lmb1 = sigU1;
lmb1b = lmb1*sqrt(2)/(1 - sqrt(sigAlfa));
lmb2 = (1 + sqrt(sigAlfa))*lmb1b/sqrt(2);
k2Aj3 = (1/2)*fM12(a1,a2)*(lmb1 + lmb2);
k1Aj3 = fM123(a1,a2,a3)*(lmb1*lmb2)/k2Aj3;

%(4) Ajuste completo (inércias e rigidezes):
%(4.1) Ajuste sem fixação da inércia total:
falfaAj = @(alfa1,alfa2,alfa3)(1 - 2*fI123(alfa1,alfa2,alfa3)/...
    (fI12(alfa1,alfa2)*fI23(alfa2,alfa3)));
falfaAjaTil = @(alfaTil)(falfaAj(falfa1Aj(alfaTil),falfa2Aj(alfaTil),falfa3Aj(alfaTil)));

%Dinâmica de torção:
lmb1 = sigT1;
lmb2 = sigT2;
sigTalfa0 = ((lmb2 - lmb1)/(lmb2 + lmb1))^2;

fsigCusto = @(alfaTil)((sigTalfa0 - fsigTalfaTil(alfaTil))^2);
x0 = ones(3,1);
xOtm = fminsearch(fsigCusto,x0,optimset('TolX',1e-10));
%{
teste = sigTalfa0 - fsigTalfaTil(xOtm);
%}
a1 = falfa1Aj(xOtm);
a2 = falfa2Aj(xOtm);
a3 = falfa3Aj(xOtm);
I1aj = a1*Ip1;
I2aj = a2*(Ib2 + dm*l2^2);
I3aj = a3*I3;
K2Aj = (1/2)*fI12(a1,a2)*(lmb1 + lmb2);
K1Aj = fI123(a1,a2,a3)*(lmb1*lmb2)/K2Aj;
%{
MtAj = diag([I3aj, I2aj, I1aj]);
KtAj = [ K1Aj,     -K1Aj,     0; 
        -K1Aj, K1Aj+K2Aj, -K2Aj; 
            0,     -K2Aj,  K2Aj];
teste = eig(MtAj\KtAj);
teste = sort(teste);
%}
%(4.2) Ajuste com fixação da inércia total:
Id2 = eye(2);
IpTotal = Ip1 + (Ib2 + dm*l2^2) + I3;
fI3aj5 = @(alfa1,alfa2)(IpTotal - fI1(alfa1) - fI2(alfa2));
fI23aj5 = @(alfa1,alfa2)(((fI2(alfa2))^-1 + fI3aj5(alfa1,alfa2)^-1)^-1);
fI123aj5 = @(alfa1,alfa2)(...
    (1/(fI3aj5(alfa1,alfa2)*fI1(alfa1)) + 1/((fI2(alfa2))*fI1(alfa1)) + 1/((fI2(alfa2))*fI3aj5(alfa1,alfa2)))^-1);
fsigTalfaAj5 = @(alfa1,alfa2)(1 - fI12(alfa1,alfa2)*fI23aj5(alfa1,alfa2)/fI123aj5(alfa1,alfa2));
falfa1Aj5 = @(alfaTil)(Id2(1,:)*alfaTil);
falfa2Aj5 = @(alfaTil)(Id2(2,:)*alfaTil);
fsigTalfaTilAj5 = @(alfaTil)(fsigTalfaAj5(falfa1Aj5(alfaTil),falfa2Aj5(alfaTil)));
fsigCusto = @(alfaTil)((sigTalfa0 - fsigTalfaTilAj5(alfaTil))^2);
x0 = ones(2,1);
xOtm = fminsearch(fsigCusto,x0,optimset('TolX',1e-10));
%{
teste = sigTalfa0 - fsigTalfaTilAj5(xOtm);
%}
a1 = falfa1Aj5(xOtm);
a2 = falfa2Aj5(xOtm);
I1aj5 = fI1(a1);
I2aj5 = fI2(a2);
I3aj5 = fI3aj5(a1,a2);
K2Aj5 = (1/2)*fI12(a1,a2)*(lmb1 + lmb2);
K1Aj5 = fI123aj5(a1,a2)*(lmb1*lmb2)/K2Aj5;
%{
MtAj5 = diag([I3aj5, I2aj5, I1aj5]);
KtAj5 = [ K1Aj5,      -K1Aj5,      0; 
        -K1Aj5, K1Aj5+K2Aj5, -K2Aj5; 
             0,      -K2Aj5,  K2Aj5];
teste = eig(MtAj5\KtAj5);
teste = sort(teste);
%}

%%
%Dinâmica longitudinal:
lmb1 = sigU1;
lmb2 = sigU2;
sigUalfa0 = ((lmb2 - lmb1)/(lmb2 + lmb1))^2;

fsigCusto = @(alfaTil)((sigUalfa0 - fsigUalfaTil(alfaTil))^2);
x0 = ones(3,1);
xOtm = fminsearch(fsigCusto,x0,optimset('TolX',1e-10));
%{
teste = sigUalfa0 - fsigUalfaTil(xOtm);
%}
a1 = falfa1Aj(xOtm);
a2 = falfa2Aj(xOtm);
a3 = falfa3Aj(xOtm);
M1aj = a1*M1;
M2aj = a2*(Mb2 + dm);
M3aj = a3*M3;
k2Aj = (1/2)*fM12(a1,a2)*(lmb1 + lmb2);
k1Aj = fM123(a1,a2,a3)*(lmb1*lmb2)/k2Aj;
%{
MuAj = diag([M3aj, M2aj, M1aj]);
KuAj = [ k1Aj,     -k1Aj,     0; 
        -k1Aj, k1Aj+k2Aj, -k2Aj; 
            0,     -k2Aj,  k2Aj];
teste = eig(MuAj\KuAj);
teste = sort(teste);
%}

%Ajuste com fixação da massa total:
Id2 = eye(2);
Mtotal = M1 + (Mb2 + dm) + M3;
fM3aj5 = @(alfa1,alfa2)(Mtotal - fM1(alfa1) - fM2(alfa2));
fM23aj5 = @(alfa1,alfa2)(((fM2(alfa2))^-1 + fM3aj5(alfa1,alfa2)^-1)^-1);
fM123aj5 = @(alfa1,alfa2)(...
    (1/(fM3aj5(alfa1,alfa2)*fM1(alfa1)) + 1/((fM2(alfa2))*fM1(alfa1)) + 1/((fM2(alfa2))*fM3aj5(alfa1,alfa2)))^-1);
fsigUalfaAj5 = @(alfa1,alfa2)(1 - fM12(alfa1,alfa2)*fM23aj5(alfa1,alfa2)/fM123aj5(alfa1,alfa2));
fsigUalfaTilAj5 = @(alfaTil)(fsigUalfaAj5(falfa1Aj5(alfaTil),falfa2Aj5(alfaTil)));
fsigCusto = @(alfaTil)((sigUalfa0 - fsigUalfaTilAj5(alfaTil))^2);
x0 = ones(2,1);
xOtm = fminsearch(fsigCusto,x0,optimset('TolX',1e-10));
%{
teste = sigUalfa0 - fsigUalfaTilAj5(xOtm);
%}
a1 = falfa1Aj5(xOtm);
a2 = falfa2Aj5(xOtm);
M1aj5 = fM1(a1);
M2aj5 = fM2(a2);
M3aj5 = fM3aj5(a1,a2);
k2Aj5 = (1/2)*fM12(a1,a2)*(lmb1 + lmb2);
k1Aj5 = fM123aj5(a1,a2)*(lmb1*lmb2)/k2Aj5;
%{
MuAj5 = diag([M3aj5, M2aj5, M1aj5]);
KuAj5 = [ k1Aj5,      -k1Aj5,      0; 
        -k1Aj5, k1Aj5+k2Aj5, -k2Aj5; 
             0,      -k2Aj5,  k2Aj5];
teste = eig(MuAj5\KuAj5);
teste = sort(teste);
%}

%%
%Rigidezes aproximadas e inércias originais para o metamodelo:
K1Aprox = G*Jp/L1;
K2Aprox = G*Jp/(L - L1);
k1Aprox = E*A/L1;
k2Aprox = E*A/(L - L1);

%%
%{
%Metamodelo:
betaUT = 0;
Mmm = diag([M3, I3, Mb2 + dm, Ib2 + dm*l2^2, M1, I1]);
Kmm = [k1,   0,   -k1,     0,   0,   0;
        0,  K1,     0,   -K1,   0,   0;
      -k1,   0, k1+k2,     0, -k2,   0; 
        0, -K1,     0, K1+K2,   0, -K2; 
        0,   0,   -k2,     0,  k2,   0; 
        0,   0,     0,   -K2,   0,  K2];
CeMm = betaUT*Kmm;
CcMm = diag([b1, c1, 0, 0, b3, c3]);
Cmm = CeMm + CcMm;
Id6 = eye(6);
Bmm = Id6(:,[2,1]);

%%
%Equações de estado:
Amm = [zeros(6), Id6; -Mmm\Kmm, -Mmm\Cmm];
BbMm = [zeros(6,2); Mmm\Bmm];
C1mm = [zeros(1,6), Id6(6,:)];
C2mm = [Id6(5,:), zeros(1,6)];

%%
%Forma normal do metamodelo:
%y1 = teta3p:
beta0 = C1mm;
beta1 = beta0*Amm;
beta0Bb = beta0*BbMm;
beta2 = beta1*Amm;
beta1Bb = beta1*BbMm;
beta3 = beta2*Amm;
beta2Bb = beta2*BbMm;
beta4 = beta3*Amm;
beta3Bb = beta3*BbMm;
beta5 = beta4*Amm;
beta4Bb = beta4*BbMm;
%y2 = Nr(ub3):
alfa0 = kNM1*C2mm;
alfa1 = alfa0*Amm;
alfa0Bb = alfa0*BbMm;
alfa2 = alfa1*Amm;
alfa1Bb = alfa1*BbMm;
alfa3 = alfa2*Amm;
alfa2Bb = alfa2*BbMm;
alfa4 = alfa3*Amm;
alfa3Bb = alfa3*BbMm;
alfa5 = alfa4*Amm;
alfa4Bb = alfa4*BbMm;
alfa6 = alfa5*Amm;
alfa5Bb = alfa5*BbMm;

%%
%Dinâmica interna
r = 11;
n = 6;
nMi = 2*n - r;

Id12 = eye(12);
AmiAstTotal = Id12([1:6,9:12],:);
BbMmAst = [zeros(6,2); Bmm];
%{
%teste:
AmiTotalBbAst = AmiAstTotal*BbMmAst;
%}

sigDif0 = [beta0; beta1; beta2; beta3; beta4; alfa0; alfa1; alfa2; alfa3; alfa4; alfa5];
[ AmiOp, xOp, x , CUSTO, Rsig, FMI, CSIG, REALMax ] = fAGAst( sigDif0, AmiAstTotal, Amm, Mmm );
%}
%Obs.:
%(1) Essa abordagem não deu certo
%(2) Próxima tentativa: trabalhar com dinâmicas de forma independente

%%
%Inercias e rigidezes utilizadas no metamodelo:
switch IndAjustMod
    case 'sa'
        %Sem ajuste:
        K1mm = K1Aprox;
        K2mm = K2Aprox;
        k1mm = k1Aprox;
        k2mm = k2Aprox;
        I1mm = Ip1;
        I2mm = Ib2 + dm*l2^2;
        I3mm = I3;
        M1mm = M1;
        M2mm = Mb2 + dm;
        M3mm = M3;
        Mu = diag([M3mm, M2mm, M1mm]);
        Ku = [ k1mm,     -k1mm,     0;
              -k1mm, k1mm+k2mm, -k2mm; 
                  0,     -k2mm,  k2mm];
        Mt = diag([I3mm, I2mm, I1mm]);
        Kt = [ K1mm,     -K1mm,     0; 
              -K1mm, K1mm+K2mm, -K2mm; 
                  0,     -K2mm,  K2mm];
    case 1
        %Ajuste 1:
        K1mm = K1Aj1;
        K2mm = K2Aj1;
        k1mm = k1Aj1;
        k2mm = k2Aj1;
        I1mm = Ip1;
        I2mm = Ib2 + dm*l2^2;
        I3mm = I3;
        M1mm = M1;
        M2mm = Mb2 + dm;
        M3mm = M3;
        Mu = diag([M3mm, M2mm, M1mm]);
        Ku = [ k1mm,     -k1mm,     0;
              -k1mm, k1mm+k2mm, -k2mm; 
                  0,     -k2mm,  k2mm];
        Mt = diag([I3mm, I2mm, I1mm]);
        Kt = [ K1mm,     -K1mm,     0; 
              -K1mm, K1mm+K2mm, -K2mm; 
                  0,     -K2mm,  K2mm];
    case 2
        %Ajuste 2:
        K1mm = K1Aj2;
        K2mm = K2Aj2;
        k1mm = k1Aj2;
        k2mm = k2Aj2;
        I1mm = Ip1;
        I2mm = Ib2 + dm*l2^2;
        I3mm = I3;
        M1mm = M1;
        M2mm = Mb2 + dm;
        M3mm = M3;
        Mu = diag([M3mm, M2mm, M1mm]);
        Ku = [ k1mm,     -k1mm,     0;
              -k1mm, k1mm+k2mm, -k2mm; 
                  0,     -k2mm,  k2mm];
        Mt = diag([I3mm, I2mm, I1mm]);
        Kt = [ K1mm,     -K1mm,     0; 
              -K1mm, K1mm+K2mm, -K2mm; 
                  0,     -K2mm,  K2mm];
    case 3
        %Ajuste 3:
        K1mm = K1Aj3;
        K2mm = K2Aj3;
        k1mm = k1Aj3;
        k2mm = k2Aj3;
        I1mm = Ip1;
        I2mm = Ib2 + dm*l2^2;
        I3mm = I3;
        M1mm = M1;
        M2mm = Mb2 + dm;
        M3mm = M3;
        Mu = diag([M3mm, M2mm, M1mm]);
        Ku = [ k1mm,     -k1mm,     0;
              -k1mm, k1mm+k2mm, -k2mm; 
                  0,     -k2mm,  k2mm];
        Mt = diag([I3mm, I2mm, I1mm]);
        Kt = [ K1mm,     -K1mm,     0; 
              -K1mm, K1mm+K2mm, -K2mm; 
                  0,     -K2mm,  K2mm];
    case 4
        %Ajuste 4:
        K1mm = K1Aj;
        K2mm = K2Aj;
        k1mm = k1Aj;
        k2mm = k2Aj;
        I1mm = I1aj;
        I2mm = I2aj;
        I3mm = I3aj;
        M1mm = M1aj;
        M2mm = M2aj;
        M3mm = M3aj;
        Mu = diag([M3mm, M2mm, M1mm]);
        Ku = [ k1mm,     -k1mm,     0;
              -k1mm, k1mm+k2mm, -k2mm; 
                  0,     -k2mm,  k2mm];
        Mt = diag([I3mm, I2mm, I1mm]);
        Kt = [ K1mm,     -K1mm,     0; 
              -K1mm, K1mm+K2mm, -K2mm; 
                  0,     -K2mm,  K2mm];
    case 5
        %Ajuste 5:
        K1mm = K1Aj5;
        K2mm = K2Aj5;
        k1mm = k1Aj5;
        k2mm = k2Aj5;
        I1mm = I1aj5;
        I2mm = I2aj5;
        I3mm = I3aj5;
        M1mm = M1aj5;
        M2mm = M2aj5;
        M3mm = M3aj5;
        Mu = diag([M3mm, M2mm, M1mm]);
        Ku = [ k1mm,     -k1mm,     0;
              -k1mm, k1mm+k2mm, -k2mm; 
                  0,     -k2mm,  k2mm];
        Mt = diag([I3mm, I2mm, I1mm]);
        Kt = [ K1mm,     -K1mm,     0; 
              -K1mm, K1mm+K2mm, -K2mm; 
                  0,     -K2mm,  K2mm];
    otherwise
        disp('other value')
end
argumentos.K1mm = K1mm;
argumentos.K2mm = K2mm;
argumentos.k1mm = k1mm;
argumentos.k2mm = k2mm;
argumentos.I1mm = I1mm;
argumentos.I2mm = I2mm;
argumentos.I3mm = I3mm;
argumentos.M1mm = M1mm;
argumentos.M2mm = M2mm;
argumentos.M3mm = M3mm;

%%
betaU = 0;
CeU = betaU*Ku;
CcU = diag([b1, 0, b3]);
Cu = CeU + CcU;
Bu = Id3(:,1);
W = Mu*ones(3,1)*g;
%%
%Ponto de equilíbrio uTilEq para u1=0:
SigmaR = Id3(:,[2,3]);
Kur = SigmaR.'*Ku*SigmaR;
Wr = SigmaR.'*W;
uTilEqr = Kur\Wr;
uTilEq = SigmaR*uTilEqr;
W0 = W - Ku*uTilEq;
%%
%Força de reação à restrição u1=0: Fu0
%Obs.: minimizar CUSTO(Fu0) = (Bu*Fu0 - W0).'*(Bu*Fu0 - W0)
aF = Bu.'*Bu;
bF = 2*Bu.'*W0;
cF = W0.'*W0;
DeltaF = bF^2 - 4*aF*cF;
Fu0min = -bF/(2*aF);
CUSTOmin = -DeltaF/(4*aF);
Fu0 = Fu0min;
%%
%Força de contato:
u3Eq = uTilEq(3,1);
uC = u3Eq;
nNu3 = (1/pi)*(1 - uC/delta0)^(-1);
switch opcaoSignAt
    case 1
        fNu3 = @(u3)(kNM1*(u3 - uC).*double(u3>=uC));
        fdNu3du3 = @(u3)(kNM1*double(u3>=uC));
    case 2
        fNu3 = @(u3)(kNM1*(u3 - uC).*(tanh(nNu3*pi*(u3/delta0 - 1)) + 1));
        fdNu3du3 = @(u3)(kNM1.*(tanh(nNu3*pi*(u3/delta0 - 1)) + 1) +...
            kNM1*(u3 - uC)*(-(pi*nNu3*(tanh(pi*nNu3*(u3/delta0 - 1))^2 - 1))/delta0));
    otherwise
        disp('other value')
end
repositorio.fNu3 = fNu3;
repositorio.fdNu3du3 = fdNu3du3;
fNub3 = @(ub3)(fNu3(ub3 + u3Eq));
fub3 = @(ubTil)(Id3(3,:)*ubTil);
fNuTil = @(ubTil)(Id3(:,3)*fNub3(fub3(ubTil)));
dub3dubTilT = Id3(3,:);
fdNub3dubTilT = @(ubTil)(fdNu3du3(fub3(ubTil) + u3Eq)*dub3dubTilT);
fdNudubTilT = @(ubTil)(Id3(:,3)*fdNub3dubTilT(ubTil));

%%
%Equações de estado: 
%Obs.: 
%(1) Mu*ubppTilu + Cu*ubpTilu + Ku*ubTil + fNuTil(ubTil) = Bu*Fu
%(2) Zp = Au*Z + FuNr(Z) + BbU*Fu
Kbu = Ku + diag([0, 0, kNM1]);
Abu = [zeros(3), Id3; -Mu\Kbu, -Mu\Cu];
%Au = [zeros(3), Id3; -Mu\Ku, -Mu\Cu];
Au = Abu;
BbU = [zeros(3,1); Mu\Bu];
Id6 = eye(6);
AuTil = Id6(:,1:3);
AvTil = Id6(:,4:6);
fubTil = @(Z)(AuTil.'*Z);
fNz = @(Z)(fNuTil(fubTil(Z)));
FuNrz = @(Z)([zeros(3,1); fNz(Z)]);
C2 = [Id3(3,:), zeros(1,3)];
%%
%Forma normal do metamodelo:
%y2 = Nr:
alfam1 = kNM1*C2*(Abu\Id6);
alfam1Bb = alfam1*BbU;
alfa0 = kNM1*C2;
alfa1 = alfa0*Au;
alfa0Bb = alfa0*BbU;
alfa2 = alfa1*Au;
alfa1Bb = alfa1*BbU;
alfa3 = alfa2*Au;
alfa2Bb = alfa2*BbU;
alfa4 = alfa3*Au;
alfa3Bb = alfa3*BbU;
alfa5 = alfa4*Au;
alfa4Bb = alfa4*BbU;
alfa6 = alfa5*Au;
alfa5Bb = alfa5*BbU;

%Versão descontínua adaptada:
%Obs.: 
%(1) Mu*ubppTil + Cu*ubpTil + fKbu(ubTil)*ubTil = Bu*Fu
fkNM1 = @(u3)(kNM1*double(u3>=uC));
fkNM1u = @(ub3)(fkNM1(ub3 + u3Eq));
fKbu = @(ubTil)(Ku + diag([0, 0, fkNM1u(fub3(ubTil))]));
fKbuz = @(Z)(fKbu(fubTil(Z)));
fAbu = @(Z)([zeros(3), Id3; -Mu\fKbuz(Z), -Mu\Cu]);
falfa1 = @(Z)(alfa0*fAbu(Z));
falfa2 = @(Z)(falfa1(Z)*fAbu(Z));
falfa3 = @(Z)(falfa2(Z)*fAbu(Z));
falfa4 = @(Z)(falfa3(Z)*fAbu(Z));
falfa5 = @(Z)(falfa4(Z)*fAbu(Z));
falfa6 = @(Z)(falfa5(Z)*fAbu(Z));
falfa5Bb = @(Z)(falfa5(Z)*BbU);

%%
%y1 = tetap3
betaT = 0;
CeT = betaT*Kt;
CcT = diag([c1, 0, c3]);
Ct = CeT + CcT;
Bt = Id3(:,1);

%%
%Vetores de atrito
%Obs.: uTil = [teta1, teta2, teta3].'
fteta3T = @(upTil)(Id3(3,:)*upTil);
dteta3duTilT = Id3(3,:);
%Interação de atrito com modelo de atrito contínuo na origem:
fTatc = @(tetap3,y2)(miAtr(tetap3).*y2*R3);
fTatcuTil = @(upTil,y2)(miAtr(fteta3T(upTil))*y2*R3);
fdTatcduTilT = @(upTil,y2)(dmiAtrdtp(fteta3T(upTil))*dteta3duTilT*y2*R3);
FatTc = @(upTil,y2)(...
    [0;
     0;
     fTatcuTil(upTil,y2)]);
dFatTcduTilT = @(upTil,y2)(...
    [zeros(1,3);
     zeros(1,3);
     fdTatcduTilT(upTil,y2)]);
%Interação de atrito com modelo de atrito descontínuo na origem:
fTatDc = @(upTil,y2)(miDc(fteta3T(upTil))*y2*R3);
fdTatDcduTilT = @(upTil,y2)(dmidtXpDc(fteta3T(upTil))*dteta3duTilT*y2*R3);
FatTdc = @(upTil,y2)(...
    [0;
     0;
     fTatDc(upTil,y2)]);
dFatTdcduTilT = @(upTil,y2)(...
    [zeros(1,3);
     zeros(1,3);
     fdTatDcduTilT(upTil,y2)]);
%%
%Equações de estado:
At = [zeros(3), Id3; -Mt\Kt, -Mt\Ct];
argumentos.At = At;
BbT = [zeros(3,1); Mt\Bt];
C1 = [zeros(1,3), Id3(3,:)];
%%
%Forma normal do metamodelo:
%y1 = teta3p:
beta0 = C1;
beta1 = beta0*At;
beta0Bb = beta0*BbT;
beta2 = beta1*At;
beta1Bb = beta1*BbT;
beta3 = beta2*At;
beta2Bb = beta2*BbT;
beta4 = beta3*At;
beta3Bb = beta3*BbT;
beta5 = beta4*At;
beta4Bb = beta4*BbT;
%%
%Dinâmica interna na dinâmica de torção:
%Obs.: 
r = 5;
n = 3;
nMi = 2*n - r; %nM1 = 1
AmiAstTotal = Id6([1:3,5:6],:);
BbTast = [zeros(3,1); Bt];
%{
%Teste:
AmiTotalBbTast = AmiAstTotal*BbTast;
%}
sigDif0t = [beta0; beta1; beta2; beta3; beta4];
[ AmiT, xOp, x , CUSTO, Rsig, FMI, CSIG, REALMax ] = fAGAst( sigDif0t, AmiAstTotal, At, Mt );
%%
%Difeomorfismos:
sigDifu = [alfa0; alfa1; alfa2; alfa3; alfa4; alfa5];
condSigDifu = cond(sigDifu);
sigDift = [beta0; beta1; beta2; beta3; beta4; AmiT];
condSigDift = cond(sigDift);
fi1DifTx = @(Z,y2)(fFi1DifTx( Z, y2, argumentos, repositorio ));
fi2DifTx = @(Z,y2,yp2)(fFi2DifTx( Z, y2, yp2, argumentos, repositorio ));
fi3DifTx = @(Z,y2,yp2,y2p2)(fFi3DifTx( Z, y2, yp2, y2p2, argumentos, repositorio ));
fi4DifTx = @(Z,y2,yp2,y2p2,y3p2)(fFi4DifTx( Z, y2, yp2, y2p2, y3p2, argumentos, repositorio ));
fi5DifTx = @(Z,y2,yp2,y2p2,y3p2,y4p2)(fFi5DifTx( Z, y2, yp2, y2p2, y3p2, y4p2, argumentos, repositorio ));
FiDifTx = @(Z,y2,yp2,y2p2,y3p2)(...
    [0; 
     fi1DifTx(Z,y2); 
     fi2DifTx(Z,y2,yp2); 
     fi3DifTx(Z,y2,yp2,y2p2); 
     fi4DifTx(Z,y2,yp2,y2p2,y3p2); 
     0]);
dfi1DifTxdZT = @(Z,y2)(...
    fdFi1DifTxdZT( Z, y2, argumentos, repositorio ));
dfi2DifTxdZT = @(Z,y2,yp2)(...
    fdFi2DifTxdZT( Z, y2, yp2, argumentos, repositorio ));
dfi3DifTxdZT = @(Z,y2,yp2,y2p2)(...
    fdFi3DifTxdZT( Z, y2, yp2, y2p2, argumentos, repositorio ));
dfi4DifTxdZT = @(Z,y2,yp2,y2p2,y3p2)(...
    fdFi4DifTxdZT( Z, y2, yp2, y2p2, y3p2, argumentos, repositorio ));
dfi5DifTxdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(...
    fdFi5DifTxdZT( Z, y2, yp2, y2p2, y3p2, y4p2, argumentos, repositorio ));
dFiDifTxdZT = @(Z,y2,yp2,y2p2,y3p2)(...
    [zeros(1,6); 
     dfi1DifTxdZT(Z,y2); 
     dfi2DifTxdZT(Z,y2,yp2); 
     dfi3DifTxdZT(Z,y2,yp2,y2p2); 
     dfi4DifTxdZT(Z,y2,yp2,y2p2,y3p2); 
     zeros(1,6)]);
QsiDifTx = @(Z,y2,yp2,y2p2,y3p2)(sigDift*Z - FiDifTx(Z,y2,yp2,y2p2,y3p2));
repositorio.QsiDifTx = QsiDifTx;
dQsiDifTxdZT = @(Z,y2,yp2,y2p2,y3p2)(sigDift - dFiDifTxdZT(Z,y2,yp2,y2p2,y3p2));
repositorio.dQsiDifTxdZT = dQsiDifTxdZT;
iQsiDifTx = @(X,y2,yp2,y2p2,y3p2,Z0)(fiQsiDifTx( X, y2, yp2, y2p2, y3p2, Z0, argumentos, repositorio ));
repositorio.iQsiDifTx = iQsiDifTx;
diQsiDifTxdXT = @(X,y2,yp2,y2p2,y3p2,Z0)(...
    dQsiDifTxdZT(iQsiDifTx(X,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2)\Id6);

%Difeomorfismos, considerando contato e não contato de atrito:
fsigDifu = @(Z)([alfa0; falfa1(Z); falfa2(Z); falfa3(Z); falfa4(Z); falfa5(Z)]);
QsiDifU = @(Z)(fsigDifu(Z)*Z);
repositorio.QsiDifU = QsiDifU;
dQsiDifUdZT = @(Z)(fsigDifu(Z));
repositorio.dQsiDifUdZT = dQsiDifUdZT;
%{
iQsiDifU = @(X,Z0)(fiQsiDifU( X, Z0, repositorio ));
diQsiDifUdXT = @(X,Z0)(...
    dQsiDifUdZT(iQsiDifU(X,Z0))\Id6);
%}
iQsiDifU = @(X,Z)(fsigDifu(Z)\X);
diQsiDifUdXT = @(X,Z)(...
    dQsiDifUdZT(iQsiDifU(X,Z))\Id6);
repositorio.iQsiDifU = iQsiDifU;
repositorio.diQsiDifUdXT = diQsiDifUdXT;

%%
%Trajetória desejada:
%Força normal:
%Nref = Nref + kNM1*uNM1aj;
wRef = argumentos.wRef;
gama_TD = 0.1; %(N/(rad*s^-1)) Parâmetro de ajuste para cada ponto de operação
Psi = @(beta)(beta*exp(wRef/beta) - wRef - Nref/gama_TD - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
beta0TD = 10^6; %rad/s
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
while fiBeta(beta0TD) <= beta0TD/2
    beta0TD = beta0TD/2;
end
[ beta_TD ] = f_Newton_Raphson( beta0TD, Psi, dPsidbeta, argumentos );

Nr0b = @(w)(gama_TD*beta_TD*exp(w/beta_TD) - gama_TD*(w + beta_TD));
repositorio.Nr0b = Nr0b;
dNr0bdw = @(w)(gama_TD*(exp(w/beta_TD) - 1));
repositorio.dNr0bdw = dNr0bdw;
d2Nr0bdw2 = @(w)((gama_TD*exp(w/beta_TD))/beta_TD);
repositorio.d2Nr0bdw2 = d2Nr0bdw2;
d3Nr0bdw3 = @(w)((gama_TD*exp(w/beta_TD))/beta_TD^2);
repositorio.d3Nr0bdw3 = d3Nr0bdw3;
d4Nr0bdw4 = @(w)((gama_TD*exp(w/beta_TD))/beta_TD^3);
repositorio.d4Nr0bdw4 = d4Nr0bdw4;
d5Nr0bdw5 = @(w)((gama_TD*exp(w/beta_TD))/beta_TD^4);
repositorio.d5Nr0bdw5 = d5Nr0bdw5;
d6Nr0bdw6 = @(w)((gama_TD*exp(w/beta_TD))/beta_TD^5);
repositorio.d6Nr0bdw6 = d6Nr0bdw6;
d7Nr0bdw7 = @(w)((gama_TD*exp(w/beta_TD))/beta_TD^6);
repositorio.d7Nr0bdw7 = d7Nr0bdw7;

%Força normal estimada:
%Obs.: 
%(1) para compensação de movimento longitudinal devido ao mov lateral
%(2) aplicar na restrição do modelo EF
NrefAj = Nref + kNM1*uNM1aj;
wRef = argumentos.wRef;
gama_TD = 0.1; %(N/(rad*s^-1)) Parâmetro de ajuste para cada ponto de operação
Psi = @(beta)(beta*exp(wRef/beta) - wRef - NrefAj/gama_TD - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
beta0TD = 10^6; %rad/s
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
while fiBeta(beta0TD) <= beta0TD/2
    beta0TD = beta0TD/2;
end
[ beta_TDaj ] = f_Newton_Raphson( beta0TD, Psi, dPsidbeta, argumentos );
Nr0bAj = @(w)(gama_TD*beta_TDaj*exp(w/beta_TDaj) - gama_TD*(w + beta_TDaj));
dNr0bAjdw = @(w)(gama_TD*(exp(w/beta_TDaj) - 1));
d2Nr0bAjdw2 = @(w)((gama_TD*exp(w/beta_TDaj))/beta_TDaj);

%Força normal estimada com correção de rotação:
%Obs.: 
%(1) para compensação de movimento longitudinal devido ao mov lateral
%(2) aplicar na restrição do metamodelo
nubMax = 6;
dubNM1max = nubMax*uNM1aj0;
NrefAjMM = Nref + kNM1*(uNM1aj - dubNM1max);
wRef = argumentos.wRef;
gama_TD = 0.1; %(N/(rad*s^-1)) Parâmetro de ajuste para cada ponto de operação
Psi = @(beta)(beta*exp(wRef/beta) - wRef - NrefAjMM/gama_TD - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
beta0TD = 10^6; %rad/s
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
while fiBeta(beta0TD) <= beta0TD/2
    beta0TD = beta0TD/2;
end
[ beta_TDaj ] = f_Newton_Raphson( beta0TD, Psi, dPsidbeta, argumentos );
Nr0bAjMM = @(w)(gama_TD*beta_TDaj*exp(w/beta_TDaj) - gama_TD*(w + beta_TDaj));
dNr0bAjMMdw = @(w)(gama_TD*(exp(w/beta_TDaj) - 1));
d2Nr0bAjMMdw2 = @(w)((gama_TD*exp(w/beta_TDaj))/beta_TDaj);

%(1) Trajetória desejada: 
%(1.a) y1d = tetaXp3d:
t09 = argumentos.t0_9wRef;
alfa = argumentos.alfa_trajetoria_desejada;
fi = @(nu)(tanh(alfa*nu) - (alfa*nu)/((alfa^2*nu^2)/3 + 1));
Fi = @(nu)(nu - (log(tanh(alfa*nu) + 1) + (3*log(alfa^2*nu^2 + 3))/2)/alfa);
y1d = @(t)(wRef*fi(t/t09));
Y1d = @(t)(wRef*t09*Fi(t/t09));
repositorio.y1d = y1d;
dfidnu = @(nu)((6*alfa^3*nu^2)/(alfa^2*nu^2 + 3)^2 -...
    (3*alfa)/(alfa^2*nu^2 + 3) -...
    alfa*(tanh(alfa*nu)^2 - 1));
y1pd = @(t)((wRef/t09)*dfidnu(t/t09));
repositorio.y1pd = y1pd;
d2fidnu2 = @(nu)(2*alfa^2*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1) +...
    (18*alfa^3*nu)/(alfa^2*nu^2 + 3)^2 -...
    (24*alfa^5*nu^3)/(alfa^2*nu^2 + 3)^3);
y1ppd = @(t)((wRef/t09^2)*d2fidnu2(t/t09));
repositorio.y1ppd = y1ppd;
d3fidnu3 = @(nu)((18*alfa^3)/(alfa^2*nu^2 + 3)^2 -...
    2*alfa^3*(tanh(alfa*nu)^2 - 1)^2 -...
    4*alfa^3*tanh(alfa*nu)^2*(tanh(alfa*nu)^2 - 1) -...
    (144*alfa^5*nu^2)/(alfa^2*nu^2 + 3)^3 +...
    (144*alfa^7*nu^4)/(alfa^2*nu^2 + 3)^4);
y1pppd = @(t)((wRef/t09^3)*d3fidnu3(t/t09));
repositorio.y1pppd = y1pppd;
d4fidnu4 = @(nu)(16*alfa^4*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1)^2 -...
    (360*alfa^5*nu)/(alfa^2*nu^2 + 3)^3 +...
    8*alfa^4*tanh(alfa*nu)^3*(tanh(alfa*nu)^2 - 1) +...
    (1440*alfa^7*nu^3)/(alfa^2*nu^2 + 3)^4 -...
    (1152*alfa^9*nu^5)/(alfa^2*nu^2 + 3)^5);
y1ppppd = @(t)((wRef/t09^4)*d4fidnu4(t/t09));
repositorio.y1ppppd = y1ppppd;
d5fidnu5 = @(nu)((6480*alfa^7*nu^2)/(alfa^2*nu^2 + 3)^4 -...
    16*alfa^5*(tanh(alfa*nu)^2 - 1)^3 -...
    16*alfa^5*tanh(alfa*nu)^4*(tanh(alfa*nu)^2 - 1) -...
    (360*alfa^5)/(alfa^2*nu^2 + 3)^3 -...
    (17280*alfa^9*nu^4)/(alfa^2*nu^2 + 3)^5 +...
    (11520*alfa^11*nu^6)/(alfa^2*nu^2 + 3)^6 -...
    88*alfa^5*tanh(alfa*nu)^2*(tanh(alfa*nu)^2 - 1)^2);
y1pppppd = @(t)((wRef/t09^5)*d5fidnu5(t/t09));
repositorio.y1pppppd = y1pppppd;
d6fidnu6 = @(nu)((15120*alfa^7*nu)/(alfa^2*nu^2 + 3)^4 +...
    272*alfa^6*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1)^3 +...
    32*alfa^6*tanh(alfa*nu)^5*(tanh(alfa*nu)^2 - 1) -...
    (120960*alfa^9*nu^3)/(alfa^2*nu^2 + 3)^5 +...
    (241920*alfa^11*nu^5)/(alfa^2*nu^2 + 3)^6 -...
    (138240*alfa^13*nu^7)/(alfa^2*nu^2 + 3)^7 +...
    416*alfa^6*tanh(alfa*nu)^3*(tanh(alfa*nu)^2 - 1)^2);
y1ppppppd = @(t)((wRef/t09^6)*d6fidnu6(t/t09));
repositorio.y1ppppppd = y1ppppppd;
d7fidnu7 = @(nu)((15120*alfa^7)/(alfa^2*nu^2 + 3)^4 -...
    272*alfa^7*(tanh(alfa*nu)^2 - 1)^4 -...
    64*alfa^7*tanh(alfa*nu)^6*(tanh(alfa*nu)^2 - 1) -...
    (483840*alfa^9*nu^2)/(alfa^2*nu^2 + 3)^5 +...
    (2419200*alfa^11*nu^4)/(alfa^2*nu^2 + 3)^6 -...
    (3870720*alfa^13*nu^6)/(alfa^2*nu^2 + 3)^7 +...
    (1935360*alfa^15*nu^8)/(alfa^2*nu^2 + 3)^8 -...
    2880*alfa^7*tanh(alfa*nu)^2*(tanh(alfa*nu)^2 - 1)^3 -...
    1824*alfa^7*tanh(alfa*nu)^4*(tanh(alfa*nu)^2 - 1)^2);
y1pppppppd = @(t)((wRef/t09^7)*d7fidnu7(t/t09));
repositorio.y1pppppppd = y1pppppppd;
%(1.b) y2d = ub1d:
y2d = @(t)(Nr0b(y1d(t)));
repositorio.y2d = y2d;
y2pd = @(t)(dNr0bdw(y1d(t))*y1pd(t));
repositorio.y2pd = y2pd;
y2ppd = @(t)(d2Nr0bdw2(y1d(t))*(y1pd(t))^2 +...
    dNr0bdw(y1d(t))*y1ppd(t));
repositorio.y2ppd = y2ppd;
y2pppd = @(t)(d3Nr0bdw3(y1d(t))*(y1pd(t))^3 +...
    3*d2Nr0bdw2(y1d(t))*(y1pd(t))*y1ppd(t) +...
    dNr0bdw(y1d(t))*y1pppd(t));
repositorio.y2pppd = y2pppd;
y2ppppd = @(t)(d4Nr0bdw4(y1d(t))*(y1pd(t))^4 +...
    6*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1ppd(t) +...
    3*d2Nr0bdw2(y1d(t))*(y1ppd(t))^2 +...
    4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1pppd(t) +...
    dNr0bdw(y1d(t))*y1ppppd(t));
repositorio.y2ppppd = y2ppppd;
y2pppppd = @(t)(d5Nr0bdw5(y1d(t))*(y1pd(t))^5 +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppd(t) +...
    15*d3Nr0bdw3(y1d(t))*y1pd(t)*(y1ppd(t))^2 +...
    10*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppd(t) +...
    10*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppd(t) +...
    dNr0bdw(y1d(t))*y1pppppd(t));
repositorio.y2pppppd = y2pppppd;
y2ppppppd = @(t)(d6Nr0bdw6(y1d(t))*(y1pd(t))^6 +...
    d5Nr0bdw5(y1d(t))*5*(y1pd(t))^4*y1ppd(t) +...
    10*d5Nr0bdw5(y1d(t))*(y1pd(t))^4*y1ppd(t) +...
    10*d4Nr0bdw4(y1d(t))*3*(y1pd(t))^2*(y1ppd(t))^2 +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1pppd(t) +...
    15*d4Nr0bdw4(y1d(t))*(y1pd(t))^2*(y1ppd(t))^2 +...
    15*d3Nr0bdw3(y1d(t))*(y1ppd(t))^3 +...
    15*d3Nr0bdw3(y1d(t))*y1pd(t)*2*y1ppd(t)*y1pppd(t) +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1pppd(t) +...
    10*d3Nr0bdw3(y1d(t))*2*y1pd(t)*y1ppd(t)*y1pppd(t) +...
    10*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1ppppd(t) +...
    10*d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1pppd(t) +...
    10*d2Nr0bdw2(y1d(t))*(y1pppd(t))^2 +...
    10*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1ppppd(t) +...
    4*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1ppppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1ppppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1pppppd(t) +...
    d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1ppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1ppd(t)*y1ppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1pd(t)*y1pppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1pd(t)*y1pppppd(t) +...
    dNr0bdw(y1d(t))*y1ppppppd(t));
repositorio.y2ppppppd = y2ppppppd;
y2pppppppd01 = @(t)(d7Nr0bdw7(y1d(t))*(y1pd(t))^7 + d6Nr0bdw6(y1d(t))*6*(y1pd(t))^5*y1ppd(t));
y2pppppppd02 = @(t)(d6Nr0bdw6(y1d(t))*5*(y1pd(t))^5*y1ppd(t) + d5Nr0bdw5(y1d(t))*5*4*(y1pd(t))^3*y1ppd(t)^2 + d5Nr0bdw5(y1d(t))*5*(y1pd(t))^4*y1pppd(t));
y2pppppppd03 = @(t)(10*d6Nr0bdw6(y1d(t))*(y1pd(t))^5*y1ppd(t) + 10*d5Nr0bdw5(y1d(t))*4*(y1pd(t))^3*y1ppd(t)^2 + 10*d5Nr0bdw5(y1d(t))*(y1pd(t))^4*y1pppd(t));
y2pppppppd04 = @(t)(10*d5Nr0bdw5(y1d(t))*3*(y1pd(t))^3*(y1ppd(t))^2 + 10*d4Nr0bdw4(y1d(t))*3*2*(y1pd(t))*(y1ppd(t))^3 + 10*d4Nr0bdw4(y1d(t))*3*(y1pd(t))^2*2*(y1ppd(t))*y1pppd(t));
y2pppppppd05 = @(t)(10*d5Nr0bdw5(y1d(t))*(y1pd(t))^4*y1pppd(t) + 10*d4Nr0bdw4(y1d(t))*3*(y1pd(t))^2*y1ppd(t)*y1pppd(t) + 10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t));
y2pppppppd06 = @(t)(15*d5Nr0bdw5(y1d(t))*(y1pd(t))^3*(y1ppd(t))^2 + 15*d4Nr0bdw4(y1d(t))*2*(y1pd(t))*(y1ppd(t))^3 + 15*d4Nr0bdw4(y1d(t))*(y1pd(t))^2*2*(y1ppd(t))*y1pppd(t));
y2pppppppd07 = @(t)(15*d4Nr0bdw4(y1d(t))*y1pd(t)*(y1ppd(t))^3 + 15*d3Nr0bdw3(y1d(t))*3*(y1ppd(t))^2*y1pppd(t));
y2pppppppd08 = @(t)(15*d4Nr0bdw4(y1d(t))*y1pd(t)^2*2*y1ppd(t)*y1pppd(t) + 15*d3Nr0bdw3(y1d(t))*2*y1ppd(t)^2*y1pppd(t) + 15*d3Nr0bdw3(y1d(t))*y1pd(t)*2*y1pppd(t)^2 + 15*d3Nr0bdw3(y1d(t))*y1pd(t)*2*y1ppd(t)*y1ppppd(t));
y2pppppppd09 = @(t)(10*d5Nr0bdw5(y1d(t))*(y1pd(t))^4*y1pppd(t) + 10*d4Nr0bdw4(y1d(t))*3*(y1pd(t))^2*y1ppd(t)*y1pppd(t) + 10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t));
y2pppppppd10 = @(t)(10*d4Nr0bdw4(y1d(t))*2*y1pd(t)^2*y1ppd(t)*y1pppd(t) + 10*d3Nr0bdw3(y1d(t))*2*y1ppd(t)^2*y1pppd(t) + 10*d3Nr0bdw3(y1d(t))*2*y1pd(t)*y1pppd(t)^2 + 10*d3Nr0bdw3(y1d(t))*2*y1pd(t)*y1ppd(t)*y1ppppd(t));
y2pppppppd11 = @(t)(10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t) + 10*d3Nr0bdw3(y1d(t))*2*(y1pd(t))*y1ppd(t)*y1ppppd(t) + 10*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t));
y2pppppppd12 = @(t)(10*d4Nr0bdw4(y1d(t))*y1pd(t)^2*y1ppd(t)*y1pppd(t) + 10*d3Nr0bdw3(y1d(t))*y1ppd(t)^2*y1pppd(t) + 10*d3Nr0bdw3(y1d(t))*y1pd(t)*y1pppd(t)^2 + 10*d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1ppppd(t));
y2pppppppd13 = @(t)(10*d3Nr0bdw3(y1d(t))*y1pd(t)*(y1pppd(t))^2 + 10*d2Nr0bdw2(y1d(t))*2*(y1pppd(t))*y1ppppd(t));
y2pppppppd14 = @(t)(10*d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1ppppd(t) + 10*d2Nr0bdw2(y1d(t))*y1pppd(t)*y1ppppd(t) + 10*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t));
y2pppppppd15 = @(t)(4*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t) + 4*d3Nr0bdw3(y1d(t))*2*(y1pd(t))*y1ppd(t)*y1ppppd(t) + 4*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t));
y2pppppppd16 = @(t)(4*d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1ppppd(t) + 4*d2Nr0bdw2(y1d(t))*y1pppd(t)*y1ppppd(t) + 4*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t));
y2pppppppd17 = @(t)(4*d3Nr0bdw3(y1d(t))*y1pd(t)^2*y1pppppd(t) + 4*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) + 4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppppd(t));
y2pppppppd18 = @(t)(d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t) + d3Nr0bdw3(y1d(t))*2*(y1pd(t))*y1ppd(t)*y1ppppd(t) + d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t));
y2pppppppd19 = @(t)(d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1ppppd(t) + d2Nr0bdw2(y1d(t))*y1pppd(t)*y1ppppd(t) + d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t));
y2pppppppd20 = @(t)(d3Nr0bdw3(y1d(t))*y1pd(t)^2*y1pppppd(t) + d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) + d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppppd(t));
y2pppppppd21 = @(t)(d3Nr0bdw3(y1d(t))*y1pd(t)^2*y1pppppd(t) + d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) + d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppppd(t));
y2pppppppd22 = @(t)(d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppppd(t) + dNr0bdw(y1d(t))*y1pppppppd(t));
y2pppppppd = @(t)(y2pppppppd01(t) + y2pppppppd02(t) + y2pppppppd03(t) + y2pppppppd04(t) +...
    y2pppppppd05(t) + y2pppppppd06(t) + y2pppppppd07(t) + y2pppppppd08(t) +...
    y2pppppppd09(t) + y2pppppppd10(t) + y2pppppppd11(t) + y2pppppppd12(t) +...
    y2pppppppd13(t) + y2pppppppd14(t) + y2pppppppd15(t) + y2pppppppd16(t) +...
    y2pppppppd17(t) + y2pppppppd18(t) + y2pppppppd19(t) + y2pppppppd20(t) +...
    y2pppppppd21(t) + y2pppppppd22(t));
repositorio.y2pppppppd = y2pppppppd;

y1dTil = @(t)([y1d(t); y1pd(t); y1ppd(t); y1pppd(t); y1ppppd(t)]);
y2dTil = @(t)([y2d(t); y2pd(t); y2ppd(t); y2pppd(t); y2ppppd(t); y2pppppd(t)]);
y1pdTil = @(t)([y1pd(t); y1ppd(t); y1pppd(t); y1ppppd(t); y1pppppd(t)]);
y2pdTil = @(t)([y2pd(t); y2ppd(t); y2pppd(t); y2ppppd(t); y2pppppd(t); y2ppppppd(t)]);

%%
%FORMA NORMAL, dinâmica externa:
a1fn = @(Z)(beta5*Z);
a1fnAt = @(Z,y2,yp2,y2p2,y3p2,y4p2)(beta5*Z - fi5DifTx(Z,y2,yp2,y2p2,y3p2,y4p2));
b1fn = beta4Bb;
a2fn = @(Z)(alfa6*Z);
b2fn = alfa5Bb;
%{
fa2fn = @(Z)(falfa6(Z)*Z);
fb2fn = @(Z)(falfa5Bb(Z));
%}

%sigmaA1 = 0.1;
%sigmaA1 = 1.5; %Funciona na ausência de torque de atrito (da1fn)
sigmaA1 = 5.0; %E na presença de torque de atrito como perturbação?
da1fn = @(Z)(sigmaA1*a1fn(Z));
da1fnAt = @(Z,y2,yp2,y2p2,y3p2,y4p2)(sigmaA1*a1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2));
sigmaA2 = 0.1;
da2fn = @(Z)(sigmaA2*a2fn(Z));
%{
fda2fn = @(Z)(sigmaA2*fa2fn(Z));
%}
%Impacto da diferença entre b1fnAc e b1fn sobre lei de controle:
%{
b1fnAcExp = K1*K2/(I1*I3*(Ib2 + l2^2*dm*Mb2/(Mb2 + dm)));
b1fnExp = K1*K2/(I1*I3*(Ib2 + l2^2*dm));
alfaM = Mb2/(Mb2+dm);
betaM = dm*l2^2/Ib2;
db1b1m1 = (b1fnAcExp - b1fnExp)/b1fnExp;
teste = db1b1m1 - (1 - alfaM)*betaM/(1 + alfaM*betaM);
%}
%sigmaB1 = 0.1;
sigmaB1 = 0.5; %Funciona na ausência de torque de atrito
db1fn = sigmaB1*b1fn;
sigmaB2 = 0.1;
db2fn = sigmaB2*b2fn;
%{
fdb2fn = @(Z)(sigmaB2*fb2fn(Z));
%}

Flc1 = @(Z)(abs(da1fn(Z)));
Flc1At = @(Z,y2,yp2,y2p2,y3p2,y4p2)(abs(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
Flc2 = @(Z)(abs(da2fn(Z)));
%{
fFlc2 = @(Z)(abs(fda2fn(Z)));
%}
D11fn = abs(db1fn)/b1fn;
D22fn = abs(db2fn)/b2fn;
%{
fD22fn = @(Z)(abs(fdb2fn(Z))/fb2fn(Z));
%}

y1 = @(X)(Id6(1,:)*X);
y1p = @(X)(Id6(2,:)*X);
y1pp = @(X)(Id6(3,:)*X);
y1ppp = @(X)(Id6(4,:)*X);
y1pppp = @(X)(Id6(5,:)*X);
y2 = @(X)(Id6(1,:)*X);
y2p = @(X)(Id6(2,:)*X);
y2pp = @(X)(Id6(3,:)*X);
y2ppp = @(X)(Id6(4,:)*X);
y2pppp = @(X)(Id6(5,:)*X);
y2ppppp = @(X)(Id6(6,:)*X);

%Frequência de corte:
I12 = (I1mm^-1 + I2mm^-1)^-1;
I23 = (I2mm^-1 + I3mm^-1)^-1;
I123 = (1/(I3mm*I1mm) + 1/(I2mm*I1mm) + 1/(I2mm*I3mm))^-1;
b = K2mm/I12 + K1mm/I23;
%c = K1*K2/I123;
delta = (K2mm/I12 - K1mm/I23)^2 + 4*K1mm*K2mm/I2mm^2;
lmb1t = 0;
lmb2t = (1/2)*(b - sqrt(delta));
lmb3t = (1/2)*(b + sqrt(delta));
wT2 = sqrt(lmb2t);
wT3 = sqrt(lmb3t);

M12 = (M1mm^-1 + M2mm^-1)^-1;
M23 = (M2mm^-1 + M3mm^-1)^-1;
M123 = (1/(M3mm*M1mm) + 1/(M2mm*M1mm) + 1/(M2mm*M3mm))^-1;
b = k2mm/M12 + k1mm/M23;
%c = k1*k2/M123;
delta = (k2mm/M12 - k1mm/M23)^2 + 4*k1mm*k2mm/M2mm^2;
lmb1u = 0;
lmb2u = (1/2)*(b - sqrt(delta));
lmb3u = (1/2)*(b + sqrt(delta));
wU2 = sqrt(lmb2u);
wU3 = sqrt(lmb3u);

%%
%Formulação por perturbações singulares:
sAt = size(At);
sAu = size(Au);
nBbT = length(BbT);
nBbU = length(BbU);
%epsPS = lmb2t/lmb2u;
epsPS = (wT2/wU2)*10^-3;
A11 = At; A12 = zeros(sAu);
A21 = zeros(sAt); A22 = Au*epsPS;
B1 = [BbT, zeros(nBbU,1)];
B2 = [zeros(nBbT,1), BbU*epsPS];
A0 = A11 - A12*(A22\A21);
B0 = B1 - A12*(A22\B2);
%Obs.: A12 = 0. Logo: A0 = A11 e B0 = B1

%y2 = Nr:
alfa1PS = alfa0*A22;
alfa0BbPS = alfa0*B2(:,2);
alfa2PS = alfa1PS*A22;
alfa1BbPS = alfa1PS*B2(:,2);
alfa3PS = alfa2PS*A22;
alfa2BbPS = alfa2PS*B2(:,2);
alfa4PS = alfa3PS*A22;
alfa3BbPS = alfa3PS*B2(:,2);
alfa5PS = alfa4PS*A22;
alfa4BbPS = alfa4PS*B2(:,2);
alfa6PS = alfa5PS*A22;
alfa5BbPS = alfa5PS*B2(:,2);
%%
%Difeomorfismos:
sigDifuPS = [alfa0; alfa1PS; alfa2PS; alfa3PS; alfa4PS; alfa5PS];
condSigDifuPS = cond(sigDifuPS);

y2dTilPS = @(tal)([0; 0; 0; 0; 0; 0]);
y2pdTilPS = @(tal)([0; 0; 0; 0; 0; 0]);

%%
%FORMA NORMAL, dinâmica externa (com perturbações singulares):
a2fnPS = @(Zf)(alfa6PS*Zf);
b2fnPS = alfa5BbPS;

sigmaA2 = 0.1;
da2fnPS = @(Zf)(sigmaA2*a2fnPS(Zf));
sigmaB2 = 0.1;
db2fnPS = sigmaB2*b2fnPS;

Flc2PS = @(Zf)(abs(da2fnPS(Zf)));
D22fnPS = abs(db2fnPS)/b2fnPS;

%%
%FORMA NORMAL, dinâmica externa:
avltX2 = Abvl(jPtX2,1);
avltX2M1 = Abvl(jPtX2+1,1);
nlmbTx = 0.5; 
%Obs.: 0 < nlmbTx <= 1
lmbTxRef = (nlmbTx*sqrt(avltX2) + (1 - nlmbTx)*(1/3)*sqrt(avltX2M1))*...
    double(sqrt(avltX2) < (1/3)*sqrt(avltX2M1)) +...
    (1.1*sqrt(avltX2))*...
    double(sqrt(avltX2) >= (1/3)*sqrt(avltX2M1));
avlu2 = Abvl(jPu2,1);
avlu2M1 = Abvl(jPu2+1,1);
nlmbU = 0.5; 
%Obs.: 0 < nlmbU <= 1
lmbUref = (nlmbU*sqrt(avlu2) + (1 - nlmbU)*(1/3)*sqrt(avlu2M1))*...
    double(sqrt(avlu2) < (1/3)*sqrt(avlu2M1)) +...
    (1.1*sqrt(avlu2))*...
    double(sqrt(avlu2) >= (1/3)*sqrt(avlu2M1));

lambda1 = lmbTxRef;
lambda2 = lmbUref;
lambda2PS = lambda2*epsPS;

y4r1 = @(t,X)(y1ppppd(t) -...
    4*lambda1*(y1ppp(X) - y1pppd(t)) -...
    6*lambda1^2*(y1pp(X) - y1ppd(t)) -...
    4*lambda1^3*(y1p(X) - y1pd(t)) -...
    lambda1^4*(y1(X) - y1d(t)));
y5r2 = @(t,X)(y2pppppd(t) -...
    5*lambda2*(y2pppp(X) - y2ppppd(t)) -...
    10*lambda2^2*(y2ppp(X) - y2pppd(t)) -...
    10*lambda2^3*(y2pp(X) - y2ppd(t)) -...
    5*lambda2^4*(y2p(X) - y2pd(t)) -...
    lambda2^5*(y2(X) - y2d(t)));
y5r2PS = @(tal,Xf)(-...
    5*lambda2PS*y2pppp(Xf) -...
    10*lambda2PS^2*y2ppp(Xf) -...
    10*lambda2PS^3*y2pp(Xf) -...
    5*lambda2PS^4*y2p(Xf) -...
    lambda2PS^5*y2(Xf));

%Variáveis de escorregamento:
s1 = @(t,X)(y1pppp(X) - y4r1(t,X));
s2 = @(t,X)(y2ppppp(X) - y5r2(t,X));
s2PS = @(tal,Xf)(y2ppppp(Xf) - y5r2PS(tal,Xf));

y5r1 = @(t,X)(y1pppppd(t) -...
    4*lambda1*(y1pppp(X) - y1ppppd(t)) -...
    6*lambda1^2*(y1ppp(X) - y1pppd(t)) -...
    4*lambda1^3*(y1pp(X) - y1ppd(t)) -...
    lambda1^4*(y1p(X) - y1pd(t)));
y6r2 = @(t,X)(y2ppppppd(t) -...
    5*lambda2*(y2ppppp(X) - y2pppppd(t)) -...
    10*lambda2^2*(y2pppp(X) - y2ppppd(t)) -...
    10*lambda2^3*(y2ppp(X) - y2pppd(t)) -...
    5*lambda2^4*(y2pp(X) - y2ppd(t)) -...
    lambda2^5*(y2p(X) - y2pd(t)));
y6r2PS = @(tal,Xf)(-...
    5*lambda2PS*y2ppppp(Xf) -...
    10*lambda2PS^2*y2pppp(Xf) -...
    10*lambda2PS^3*y2ppp(Xf) -...
    5*lambda2PS^4*y2pp(Xf) -...
    lambda2PS^5*y2p(Xf));
y5r1d = @(t)(y1pppppd(t));
y6r2d = @(t)(y2ppppppd(t));
y6r2dPS = @(tal)(0);

y5pr1d = @(t)(y1ppppppd(t));
y6pr2d = @(t)(y2pppppppd(t));
y6pr2dPS = @(tal)(0);

%Vetores de estado desejados:
nT = 3;
nSigDifT0 = length(y1dTil(0));

AdT = Id6(:,nSigDifT0+1:2*nT);
BdT = Id6(:,1:nSigDifT0);

fXtd = @(t,mid)(AdT*mid + BdT*y1dTil(t)); %Obs.: Xd = [ydTil(t); mid]
fDXtdDt = @(t,mipd)(AdT*mipd + BdT*y1pdTil(t));
fXud = @(t)(y2dTil(t));
fXudPS = @(tal)(y2dTilPS(tal));
fZtd = @(t,mid)(sigDift\fXtd(t,mid));
fZtdAt = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(iQsiDifTx(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0));
fZud = @(t)(sigDifu\fXud(t));
%{
fZudPlus = @(t,Z0)(iQsiDifU(fXud(t),Z0));
%}
fZudPS = @(tal)(sigDifuPS\fXudPS(tal));

UcThat = @(t,Z,X)(y5r1(t,X) - a1fn(Z));
UcThatAt = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(y5r1(t,X) - a1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2));
UcUhat = @(t,Z,X)(y6r2(t,X) - a2fn(Z));
%{
fUcUhat = @(t,Z,X)(y6r2(t,X) - fa2fn(Z));
%}
UcUhatPS = @(tal,Zf,Xf)(y6r2PS(tal,Xf) - a2fnPS(Zf));
UcTdHat = @(t,mid)(y5r1d(t) - a1fn(fZtd(t,mid)));
UcTdHatAt = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    y5r1d(t) - a1fnAt(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2));
UcUdHat = @(t)(y6r2d(t) - a2fn(fZud(t)));
%{
fUcUdHat = @(t,Z0)(y6r2d(t) - a2fn(fZudPlus(t,Z0)));
%}
UcUdHatPS = @(tal)(y6r2dPS(tal) - a2fnPS(fZudPS(tal)));

%%
%Ganhos:
etaLc1 = 1;
etaLc2 = 1;
etaLc2PS = etaLc2*epsPS;
%(1) Expressões diretas:
KlcT = @(t,Z,X)((1 - D11fn)\(Flc1(Z) + etaLc1 + D11fn*abs(UcThat(t,Z,X))));
KlcTat = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
    (1 - D11fn)\(Flc1At(Z,y2,yp2,y2p2,y3p2,y4p2) + etaLc1 +...
                 D11fn*abs(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2))));
KlcTd = @(t,mid)((1 - D11fn)\(Flc1(fZtd(t,mid)) + etaLc1 + D11fn*abs(UcTdHat(t,mid))));
KlcTdAt = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    (1 - D11fn)\(Flc1At(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2) + etaLc1 +...
                 D11fn*abs(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0))));
KlcU = @(t,Z,X)((1 - D22fn)\(Flc2(Z) + etaLc2 + D22fn*abs(UcUhat(t,Z,X))));
%{
fKlcU = @(t,Z,X)((1 - fD22fn(Z))\(Flc2(Z) + etaLc2 + fD22fn(Z)*abs(UcUhat(t,Z,X))));
%}
KlcUd = @(t)((1 - D22fn)\(Flc2(fZud(t)) + etaLc2 + D22fn*abs(UcUdHat(t))));
%{
KlcUdPlus = @(t,Z0)((1 - D22fn)\(Flc2(fZudPlus(t,Z0)) + etaLc2 + D22fn*abs(fUcUdHat(t,Z0))));
%}
KlcUps = @(tal,Zf,Xf)((1 - D22fnPS)\(Flc2PS(Zf) + etaLc2PS + D22fnPS*abs(UcUhatPS(tal,Zf,Xf))));
KlcUdPS = @(tal)((1 - D22fnPS)\(Flc2PS(fZudPS(tal)) + etaLc2PS + D22fnPS*abs(UcUdHatPS(tal))));

KfiT = lambda1/(1 - D11fn);
KfiU = lambda2/(1 - D22fn);
KfiUps = lambda2PS/(1 - D22fnPS);
KbLcT = @(t,Z,X,mid,fiBd)(KlcT(t,Z,X) - KlcTd(t,mid) + KfiT*fiBd);
KbLcTat = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    KlcTat(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) - KlcTdAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0) + KfiT*fiBd);
KbLcU = @(t,Z,X,fiBd)(KlcU(t,Z,X) - KlcUd(t) + KfiU*fiBd);
%{
KbLcUplus = @(t,Z,X,fiBd)(KlcU(t,Z,X) - KlcUdPlus(t,Z0) + KfiU*fiBd);
%}
KbLcUps = @(tal,Zf,Xf,fiBd)(KlcUps(tal,Zf,Xf) - KlcUdPS(tal) + KfiUps*fiBd);

%(2) Derivadas:
da1fndZT = beta5;
da1fnAtdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(beta5 - dfi5DifTxdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
da2fndZT = alfa6;
%{
fda2fndZT = @(Z)(falfa6(Z));
%}
da2fndZTps = alfa6PS;
dda1fndZT = sigmaA1*da1fndZT;
dda1fnAtdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(sigmaA1*da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
dda2fndZT = sigmaA2*da2fndZT;
%{
fdda2fndZT = @(Z)(sigmaA2*fda2fndZT(Z));
%}
dda2fndZTps = sigmaA2*da2fndZTps;
dFlcTdZT = @(Z)(sign(da1fn(Z))*dda1fndZT);
dFlcTatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(...
    sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2))*dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
dFlcUdZT = @(Z)(sign(da2fn(Z))*dda2fndZT);
%{
fdFlcUdZT = @(Z)(sign(da2fn(Z))*fdda2fndZT(Z));
%}
dFlcUdZTps = @(Zf)(sign(da2fnPS(Zf))*dda2fndZTps);
dZtdXtT = sigDift\Id6;
fdZtdXtT = @(X,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(X,y2,yp2,y2p2,y3p2,Z0));
dZudXuT = sigDifu\Id6;
%{
fdZudXuT = @(X,Z0)(diQsiDifUdXT(X,Z0));
%}
dZudXuTps = sigDifuPS\Id6;
dXtddt = @(t)(BdT*y1pdTil(t));
dXuddt = @(t)(y2pdTil(t));
dXuddtPS = @(tal)(y2pdTilPS(tal));
dZtddt = @(t)(dZtdXtT*dXtddt(t));
fdZtddt = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(...
    fdZtdXtT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXtddt(t));
dZuddt = @(t)(dZudXuT*dXuddt(t));
%{
fdZuddt = @(t,Z0)(fdZudXuT(fXud(t),Z0)*dXuddt(t));
%}
dZuddtPS = @(tal)(dZudXuTps*dXuddtPS(tal));

dFlcTddt = @(t,mid)(dFlcTdZT(fZtd(t,mid))*dZtddt(t));%
dFlcTddtAt = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    dFlcTatdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*fdZtddt(t,mid,y2,yp2,y2p2,y3p2,Z0));%
dFlcUddt = @(t)(dFlcUdZT(fZud(t))*dZuddt(t));%
%{
fdFlcUddt = @(t,Z0)(fdFlcUdZT(fZudPlus(t,Z0))*fdZuddt(t,Z0));%
%}
dFlcUddtPS = @(tal)(dFlcUdZTps(fZudPS(tal))*dZuddtPS(tal));
dUcTdHatdt = @(t)(y5pr1d(t) - da1fndZT*dZtddt(t));%
fdUcTdHatdt = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    y5pr1d(t) - da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*fdZtddt(t,mid,y2,yp2,y2p2,y3p2,Z0));%
dUcUdHatdt = @(t)(y6pr2d(t) - da2fndZT*dZuddt(t));%
%{
fdUcUdHatdt = @(t,Z0)(y6pr2d(t) - fda2fndZT(fZudPlus(t,Z0))*fdZuddt(t,Z0));%
%}
dUcUdHatdtPS = @(tal)(y6pr2dPS(tal) - da2fndZTps*dZuddtPS(tal));
dabsUcTdHatdt = @(t,mid)(sign(UcTdHat(t,mid))*dUcTdHatdt(t));%
fdabsUcTdHatdt = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0))*...
    fdUcTdHatdt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));%
dabsUcUdHatdt = @(t)(sign(UcUdHat(t))*dUcUdHatdt(t));%
%{
fdabsUcUdHatdt = @(t,Z0)(sign(fUcUdHat(t,Z0))*fdUcUdHatdt(t,Z0));%
%}
dabsUcUdHatdtPS = @(tal)(sign(UcUdHatPS(tal))*dUcUdHatdtPS(tal));
dKlcTddt = @(t,mid)((1 - D11fn)\(dFlcTddt(t,mid) + D11fn*dabsUcTdHatdt(t,mid)));%
fdKlcTddt = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(dFlcTddtAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0) +...
    D11fn*fdabsUcTdHatdt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));%
dKlcUddt = @(t)((1 - D22fn)\(dFlcUddt(t) + D22fn*dabsUcUdHatdt(t)));%
%{
fdKlcUddt = @(t,Z0)((1 - D22fn)\(fdFlcUddt(t,Z0) + D22fn*fdabsUcUdHatdt(t,Z0)));%
%}
dKlcUddtPS = @(tal)((1 - D22fnPS)\(dFlcUddtPS(tal) + D22fnPS*dabsUcUdHatdtPS(tal)));
dXtddmid = AdT;
dZtddmidT = dZtdXtT*dXtddmid;%
fdZtddmidT = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(fdZtdXtT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXtddmid);%
dFlcTdmidT = @(t,mid)(dFlcTdZT(fZtd(t,mid))*dZtddmidT);%
fdFlcTdmidT = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    dFlcTatdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
    fdZtddmidT(t,mid,y2,yp2,y2p2,y3p2,Z0));%
dUcTdHatdmidT = -da1fndZT*dZtddmidT;
fdUcTdHatdmidT = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
    fdZtddmidT(t,mid,y2,yp2,y2p2,y3p2,Z0));
%Obs.:
%UcTdHat = @(t,mid)(y5r1d(t) - a1fn(fZtd(t,mid)));
dabsUcTdHatdmidT = @(t,mid)(sign(UcTdHat(t,mid))*dUcTdHatdmidT);%
fdabsUcTdHatdmidT = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0))*fdUcTdHatdmidT(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));%
dKlcTddmidT = @(t,mid)((1 - D11fn)\(dFlcTdmidT(t,mid) + D11fn*dabsUcTdHatdmidT(t,mid)));%
fdKlcTddmidT = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(fdFlcTdmidT(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0) +...
    D11fn*fdabsUcTdHatdmidT(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));%
DKlcTdDt = @(t,mid,mipd)(dKlcTddt(t,mid) + dKlcTddmidT(t,mid)*mipd);%
fDKlcTdDt = @(t,mid,mipd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    fdKlcTddt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0) +...
    fdKlcTddmidT(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*mipd);%
DKlcUdDt = @(t)(dKlcUddt(t));%
%{
fDKlcUdDt = @(t,Z0)(fdKlcUddt(t,Z0));%
%}
DKlcUdDtPS = @(tal)(dKlcUddtPS(tal));

%%
%Dinâmica interna: 
%Obs.: mip = GamaMi*mi + SigmaMi*yTil

%(1) Expressões diretas:
%(1.a) Sem atrito:
GamaMi = AmiT*At*(sigDift\AdT);
SigmaMi = AmiT*At*(sigDift\BdT);

%(2) Derivadas:
fmipd = @(t,mid)(GamaMi*mid + SigmaMi*y1dTil(t));
fmipdAt = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(...
    AmiT*At*iQsiDifTx(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0));
fdmipdAtdmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(...
    AmiT*At*diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXtddmid);
fmippd = @(t,mipd)(GamaMi*mipd + SigmaMi*y1pdTil(t));
fmippdAt = @(t,mid,mipd,y2,yp2,y2p2,y3p2,Z0)(...
    AmiT*At*diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*fDXtdDt(t,mipd));

%%
%Camada limite dinâmica:
%Obs.: DfiBdT(fipBdT)*fipBdT + KfiT*fiBdT - KlcTd(t,mid) = 0
%      DfiBdU(fipBdU)*fipBdU + KfiU*fiBdU - KlcUd(t) = 0
sat = @(nu)(tanh(pi*nu));
dsatdnu = @(nu)(-pi*(tanh(pi*nu)^2 - 1));

%(1) Expressões diretas:
DfiBdT = @(fipBd)(double(fipBd>=0)/(1+D11fn) + double(fipBd<0)/(1-D11fn));
DfiBdU = @(fipBd)(double(fipBd>=0)/(1+D22fn) + double(fipBd<0)/(1-D22fn));
DfiBdUps = @(fipBd)(double(fipBd>=0)/(1+D22fnPS) + double(fipBd<0)/(1-D22fnPS));
GamaFiBdT = @(fipBd)(DfiBdT(fipBd)*fipBd);
GamaFiBdU = @(fipBd)(DfiBdU(fipBd)*fipBd);
GamaFiBdUps = @(fipBd)(DfiBdUps(fipBd)*fipBd);

%(2) Derivadas:
dGamaFiTdfipBdT = @(fipBd)(DfiBdT(fipBd));
dGamaFiUdfipBdU = @(fipBd)(DfiBdU(fipBd));
dGamaFiUdfipBdUps = @(fipBd)(DfiBdUps(fipBd));
ffippBdT = @(t,mid,mipd,fipBd)(DfiBdT(fipBd)\(-KfiT*fipBd + DKlcTdDt(t,mid,mipd)));
ffippBdTat = @(t,mid,mipd,fipBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
    DfiBdT(fipBd)\(-KfiT*fipBd + fDKlcTdDt(t,mid,mipd,y2,yp2,y2p2,y3p2,y4p2,Z0)));
ffippBdU = @(t,fipBd)(DfiBdU(fipBd)\(-KfiU*fipBd + DKlcUdDt(t)));
ffippBdUps = @(t,fipBd)(DfiBdUps(fipBd)\(-KfiUps*fipBd + DKlcUdDtPS(t)));

%%
fZ = @(uTil,upTil)(AuTil*uTil + AvTil*upTil);
fXt = @(uTil,upTil)(sigDift*fZ(uTil,upTil));
fXtAt = @(uTil,upTil,y2,yp2,y2p2,y3p2)(QsiDifTx(fZ(uTil,upTil),y2,yp2,y2p2,y3p2));
fXu = @(uTil,upTil)(sigDifu*fZ(uTil,upTil));
fXuPlus = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
fXuPS = @(uTil,upTil)(sigDifuPS*fZ(uTil,upTil));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modelo com acoplamento entre as dinâmicas lateral e longitudinal
%Obs.: uTilAc = [u1, teta1, u2, vL, wL, teta2, u3, teta3].'
Id8 = eye(8);
argumentos.Id8 = Id8;
fu1 = @(uTil)(Id8(1,:)*uTil);
repositorio.fu1 = fu1;
fteta1 = @(uTil)(Id8(2,:)*uTil);
repositorio.fteta1 = fteta1;
fu2 = @(uTil)(Id8(3,:)*uTil);
repositorio.fu2 = fu2;
fvL = @(uTil)(Id8(4,:)*uTil);
repositorio.fvL = fvL;
fwL = @(uTil)(Id8(5,:)*uTil);
repositorio.fwL = fwL;
fteta2 = @(uTil)(Id8(6,:)*uTil);
repositorio.fteta2 = fteta2;
fu3 = @(uTil)(Id8(7,:)*uTil);
repositorio.fu3 = fu3;
fteta3 = @(uTil)(Id8(8,:)*uTil);
repositorio.fteta3 = fteta3;
%%
%Elementos provenientes da energia cinética e suas derivadas:
M1ac = [M3mm,    0,    0,                  0,                 0,                  0,    0,    0; 
           0, I3mm,    0,                  0,                 0,                  0,    0,    0;
           0,    0, M2mm,                  0,                 0,                  0,    0,    0;
           0,    0,    0,               M2mm,                 0, -dm*l2*sin(alfaAc),    0,    0;
           0,    0,    0,                  0,              M2mm,  dm*l2*cos(alfaAc),    0,    0;
           0,    0,    0, -dm*l2*sin(alfaAc), dm*l2*cos(alfaAc),               I2mm,    0,    0;
           0,    0,    0,                  0,                 0,                  0, M1mm,    0;
           0,    0,    0,                  0,                 0,                  0,    0, I1mm];
fM2ac = @(uTil)(...
    [0, 0, 0,                                                     0,                                                    0,                                                     0, 0, 0; 
     0, 0, 0,                                                     0,                                                    0,                                                     0, 0, 0;
     0, 0, 0,                                                     0,                                                    0,                                                     0, 0, 0;
     0, 0, 0,                                                     0,                                                    0, -dm*l2*sin(fteta2(uTil) + alfaAc) + dm*l2*sin(alfaAc), 0, 0;
     0, 0, 0,                                                     0,                                                    0,  dm*l2*cos(fteta2(uTil) + alfaAc) - dm*l2*cos(alfaAc), 0, 0;
     0, 0, 0, -dm*l2*sin(fteta2(uTil) + alfaAc) + dm*l2*sin(alfaAc), dm*l2*cos(fteta2(uTil) + alfaAc) - dm*l2*cos(alfaAc),                                                     0, 0, 0;
     0, 0, 0,                                                     0,                                                    0,                                                     0, 0, 0;
     0, 0, 0,                                                     0,                                                    0,                                                     0, 0, 0]);
repositorio.fM2ac = fM2ac;
fGac = @(uTil,upTil)(...
    [0; 
     0; 
     0; 
     -dm*l2*fteta2(upTil)^2*cos(fteta2(uTil) + alfaAc); 
     -dm*l2*fteta2(upTil)^2*sin(fteta2(uTil) + alfaAc); 
     0; 
     0; 
     0]);
fM2upp = @(uTil,uppTil)(...
    [0;
     0;
     0;
     fteta2(uppTil)*(dm*l2*sin(alfaAc) - dm*l2*sin(alfaAc + fteta2(uTil)));
     fteta2(uppTil)*(dm*l2*cos(alfaAc + fteta2(uTil)) - dm*l2*cos(alfaAc));
     fvL(uppTil)*(dm*l2*sin(alfaAc) - dm*l2*sin(alfaAc + fteta2(uTil))) + fwL(uppTil)*(dm*l2*cos(alfaAc + fteta2(uTil)) - dm*l2*cos(alfaAc));
     0;
     0]);
%%
%Elementos provenientes da energia potencial e suas derivadas:
Kac = [k1mm,     0,     -k1mm, 0, 0,         0,     0,     0;
          0,  K1mm,         0, 0, 0,     -K1mm,     0,     0;
      -k1mm,     0, k1mm+k2mm, 0, 0,         0, -k2mm,     0; 
          0,     0,         0, 0, 0,         0,     0,     0; 
          0,     0,         0, 0, 0,         0,     0,     0; 
          0, -K1mm,         0, 0, 0, K1mm+K2mm,     0, -K2mm; 
          0,     0,     -k2mm, 0, 0,         0,  k2mm,     0; 
          0,     0,         0, 0, 0,     -K2mm,     0,  K2mm];
fFu = @(uTil)(fPCFu( uTil, argumentos, repositorio ));
Wac = [M3mm*g; 0; M2mm*g; 0; 0; 0; M1mm*g; 0];
%%
%Matriz de ganho e amortecimento:
Bac = Id8(:,[2,1]);
betaAc = 0;
CeAc = betaAc*Kac;
CcAc = diag([b1, c1, 0, 0, 0, 0, b3, c3]);
Cac = CeAc + CcAc;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vetores de atrito e rubbing
%Obs.: uTilAc = [u1, teta1, u2, vL, wL, teta2, u3, teta3].'
fN = @(uTil)(fNu3(fu3(uTil)));
%Interação de contato, sem atrito:
Fnr = @(uTil)(...
    [0; 
     0; 
     0; 
     0; 
     0; 
     0;
     fN(uTil);
     0]);
%Interação de atrito com modelo de atrito contínuo na origem:
fTatcAc = @(uTil,upTil)(miAtr(fteta3(upTil))*fNu3(fu3(uTil))*R3);
Fatc = @(uTil,upTil)(...
    [0; 
     0; 
     0; 
     0; 
     0; 
     0;
     fN(uTil);
     fTatcAc(uTil,upTil)]);
%Interação de atrito com modelo de atrito descontínuo na origem:
fTatDc = @(uTil,upTil)(miDc(fteta3(upTil))*fNu3(fu3(uTil))*R3);
FatDc = @(uTil,upTil)(...
    [0; 
     0; 
     0; 
     0; 
     0; 
     0;
     fN(uTil);
     fTatDc(uTil,upTil)]);
%Interação de rubbing com modelo de atrito contínuo na origem:
fFrubcAc = @(uTil,upTil)(...
    [0; 
     0; 
     0; 
     fFyc(fvL(uTil),fwL(uTil),fteta2(upTil),fvL(upTil),fwL(upTil)); 
     fFzc(fvL(uTil),fwL(uTil),fteta2(upTil),fvL(upTil),fwL(upTil)); 
     fT2c(fvL(uTil),fwL(uTil),fteta2(upTil),fvL(upTil),fwL(upTil)); 
     0; 
     0]);
%Interação de rubbing com modelo de atrito descontínuo na origem:
fFrubDcAc = @(uTil,upTil)(...
    [0; 
     0; 
     0; 
     fFydc(fvL(uTil),fwL(uTil),fteta2(upTil),fvL(upTil),fwL(upTil)); 
     fFzdc(fvL(uTil),fwL(uTil),fteta2(upTil),fvL(upTil),fwL(upTil)); 
     fT2dc(fvL(uTil),fwL(uTil),fteta2(upTil),fvL(upTil),fwL(upTil)); 
     0; 
     0]);
%%
%Derivadas:
%Matriz de inércia dependente de uTil:
fdM2AcuppduTilT = @(uTil,uppTil)(fPCdM2uppduTilT( uTil, uppTil, argumentos, repositorio ));
%Matriz giroscópica:
fdGacduTilT = @(uTil,upTil)(fPCdGacduTilT( uTil, upTil, argumentos, repositorio ));
fdGacdupTilT = @(uTil,upTil)(fPCdGacdupTilT( uTil, upTil, argumentos, repositorio ));
%Termo não linear da energia potencial:
fdFuduTilT = @(uTil)(fPCdFuduT( uTil, argumentos, repositorio ));
%Contato sem atrito:
fdFnrduTilT = @(uTil)(fPCdFnrduTilT( uTil, argumentos, repositorio ));
%Atrito contínuo e descontínuo na origem:
fdFatcduTilT = @(uTil,upTil)(fPCdFatcduTilT( uTil, upTil, argumentos, repositorio ));
fdFatDcduTilT = @(uTil,upTil)(fPCdFatDcduTilT( uTil, upTil, argumentos, repositorio ));
fdFatcdupTilT = @(uTil,upTil)(fPCdFatcdupTilT( uTil, upTil, argumentos, repositorio ));
fdFatDcdupTilT = @(uTil,upTil)(fPCdFatDcdupTilT( uTil, upTil, argumentos, repositorio ));
%Rubbing com modelos de atrito contínuo e descontínuo na origem:
fdFrubcduTilT = @(uTil,upTil)(fPCdFRubcduT( uTil, upTil, argumentos, repositorio ));
fdFrubDcduTilT = @(uTil,upTil)(fPCdFRubDcduT( uTil, upTil, argumentos, repositorio ));
fdFrubcdupTilT = @(uTil,upTil)(fPCdFRubcdupT( uTil, upTil, argumentos, repositorio ));
fdFrubDcdupTilT = @(uTil,upTil)(fPCdFRubDcdupT( uTil, upTil, argumentos, repositorio ));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ponto de equilíbrio uTilEq para u1=0 (Sistema com acoplamento):
SigmaRac = Id8(:,[3,7]);
Kacr = SigmaRac.'*Kac*SigmaRac;
FuAcr = @(uTilr)(SigmaRac.'*fFu(SigmaRac*uTilr));
Wacr = SigmaRac.'*Wac;
fdFuAcrduTilrT = @(uTilr)(SigmaRac.'*fdFuduTilT(SigmaRac*uTilr)*SigmaRac);
Eta0 = SigmaRac.'*zeros(8,1);
Psi = @(eta)(Kacr*eta + FuAcr(eta) - Wacr);
dPsidEtaT = @(eta)(Kacr + fdFuAcrduTilrT(eta));
[ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
uTilEqrAc = EtaConv;
uTilEqAc = SigmaRac*uTilEqrAc;
Wac0 = Wac - Kac*uTilEqAc;
%%
%Força de reação à restrição u1=0: Fu0
%Obs.: minimizar CUSTO(Fu0) = (Bu*Fu0 - W0).'*(Bu*Fu0 - W0)
aF = Bac(:,2).'*Bac(:,2);
bF = 2*Bac(:,2).'*Wac0;
cF = Wac0.'*Wac0;
DeltaF = bF^2 - 4*aF*cF;
FuAc0min = -bF/(2*aF);
CUSTOminAc = -DeltaF/(4*aF);
FuAc0 = FuAc0min;
%%
%Força de contato:
u3Eq = uTilEqAc(7,1);
uCac = u3Eq;
fNu3Ac = @(u3)(kNM1*(u3 - uCac).*double(u3>=uCac));
repositorio.fNu3Ac = fNu3Ac;
fdNu3Acdu3 = @(u3)(kNM1*double(u3>=uCac));
repositorio.fdNu3Acdu3 = fdNu3Acdu3;
fNbu3Ac = @(ub3)(fNu3Ac(ub3+u3Eq));
%%
%Equações da dinâmica em torno do ponto de equilíbrio:
fMb2ac = @(ubTil)(fM2ac(ubTil + uTilEqAc));
fGbac = @(ubTil,ubpTil)(fGac(ubTil + uTilEqAc,ubpTil));
fMb2upp = @(ubTil,ubppTil)(fM2upp(ubTil + uTilEqAc,ubppTil));
fFbu = @(ubTil)(Kac*uTilEqAc + fFu(ubTil + uTilEqAc) - Wac + Wac0);
FbNr = @(ubTil)(Fnr(ubTil + uTilEqAc));
FbAtc = @(ubTil,ubpTil)(Fatc(ubTil + uTilEqAc,ubpTil));
FbAtDc = @(ubTil,ubpTil)(FatDc(ubTil + uTilEqAc,ubpTil));
fFbRubcAc = @(ubTil,ubpTil)(fFrubcAc(ubTil + uTilEqAc,ubpTil));
fFbRubDcAc = @(ubTil,ubpTil)(fFrubDcAc(ubTil + uTilEqAc,ubpTil));
%Derivadas:
fdMb2AcuppdubT = @(ubTil,ubppTil)(fdM2AcuppduTilT(ubTil + uTilEqAc,ubppTil));
fdMb2AcuppdubppT = @(ubTil)(fMb2ac(ubTil));
fdGbacdubT = @(ubTil,ubpTil)(fdGacduTilT(ubTil + uTilEqAc,ubpTil));
fdGbacdubpT = @(ubTil,ubpTil)(fdGacdupTilT(ubTil + uTilEqAc,ubpTil));
fdFbudubT = @(ubTil)(fdFuduTilT(ubTil + uTilEqAc));
fdFbNrdubT = @(ubTil)(fdFnrduTilT(ubTil + uTilEqAc));
fdFbatcdubT = @(ubTil,ubpTil)(fdFatcduTilT(ubTil + uTilEqAc,ubpTil));
fdFbatcdubpT = @(ubTil,ubpTil)(fdFatcdupTilT(ubTil + uTilEqAc,ubpTil));
fdFbatDcdubT = @(ubTil,ubpTil)(fdFatDcduTilT(ubTil + uTilEqAc,ubpTil));
fdFbatDcdubpT = @(ubTil,ubpTil)(fdFatDcdupTilT(ubTil + uTilEqAc,ubpTil));
fdFbrubcdubT = @(ubTil,ubpTil)(fdFrubcduTilT(ubTil + uTilEqAc,ubpTil));
fdFbrubcdubpT = @(ubTil,ubpTil)(fdFrubcdupTilT(ubTil + uTilEqAc,ubpTil));
fdFbrubDcdubT = @(ubTil,ubpTil)(fdFrubDcduTilT(ubTil + uTilEqAc,ubpTil));
fdFbrubDcdubpT = @(ubTil,ubpTil)(fdFrubDcdupTilT(ubTil + uTilEqAc,ubpTil));

%%
%Parcela não linear:
fTilSmAt = @(ubTil,ubpTil,ubppTil)(...
    fGbac(ubTil,ubpTil) + fMb2upp(ubTil,ubppTil) +...
    fFbu(ubTil) +...
    FbNr(ubTil) + fFbRubcAc(ubTil,ubpTil));
fTilc = @(ubTil,ubpTil,ubppTil)(...
    fGbac(ubTil,ubpTil) + fMb2upp(ubTil,ubppTil) +...
    fFbu(ubTil) +...
    FbAtc(ubTil,ubpTil) + fFbRubcAc(ubTil,ubpTil));
fTildc = @(ubTil,ubpTil,ubppTil)(...
    fGbac(ubTil,ubpTil) + fMb2upp(ubTil,ubppTil) +...
    fFbu(ubTil) +...
    FbAtDc(ubTil,ubpTil) + fFbRubDcAc(ubTil,ubpTil));
%Correção de imprecisões, por contato do cálculo de uTilEq:
fTil0 = fTilSmAt(zeros(8,1),zeros(8,1),zeros(8,1));
fTilSmAt = @(ubTil,ubpTil,ubppTil)(fTilSmAt(ubTil,ubpTil,ubppTil) - fTil0);
fTilc = @(ubTil,ubpTil,ubppTil)(fTilc(ubTil,ubpTil,ubppTil) - fTil0);
fTildc = @(ubTil,ubpTil,ubppTil)(fTildc(ubTil,ubpTil,ubppTil) - fTil0);
%Derivadas:
dfTilSmAtdubT = @(ubTil,ubpTil,ubppTil)(...
    fdMb2AcuppdubT(ubTil,ubppTil) + fdGbacdubT(ubTil,ubpTil) +...
    fdFbudubT(ubTil) +...
    fdFbNrdubT(ubTil) + fdFbrubcdubT(ubTil,ubpTil));
dfTilcdubT = @(ubTil,ubpTil,ubppTil)(...
    fdMb2AcuppdubT(ubTil,ubppTil) + fdGbacdubT(ubTil,ubpTil) +...
    fdFbudubT(ubTil) +...
    fdFbatcdubT(ubTil,ubpTil) + fdFbrubcdubT(ubTil,ubpTil));
dfTildcdubT = @(ubTil,ubpTil,ubppTil)(...
    fdMb2AcuppdubT(ubTil,ubppTil) + fdGbacdubT(ubTil,ubpTil) +...
    fdFbudubT(ubTil) +...
    fdFbatDcdubT(ubTil,ubpTil) + fdFbrubDcdubT(ubTil,ubpTil));
dfTilSmAtdubpT = @(ubTil,ubpTil)(...
    fdGbacdubpT(ubTil,ubpTil) +...
    fdFbrubcdubpT(ubTil,ubpTil));
dfTilcdubpT = @(ubTil,ubpTil)(...
    fdGbacdubpT(ubTil,ubpTil) +...
    fdFbatcdubpT(ubTil,ubpTil) + fdFbrubcdubpT(ubTil,ubpTil));
dfTildcdubpT = @(ubTil,ubpTil)(...
    fdGbacdubpT(ubTil,ubpTil) +...
    fdFbatDcdubpT(ubTil,ubpTil) + fdFbrubDcdubpT(ubTil,ubpTil));
dfTildubppT = @(ubTil)(fdMb2AcuppdubppT(ubTil));

%Restrição em ub1:
GamaRub1Ac = Id8(1,:);

%%
%Modelo em Elementos Finitos:
%(1) Restrições via multiplicadores de Lagrange
%(2) Restriçao variável no tempo em ub1
GamaRub1 = IdNg(6*1-5,:);
GamaRtX1 = IdNg(6*1,:);
GamaRt = [GamaRub1; GamaR];
GamaRtX = [GamaRub1; GamaRtX1; GamaR];
%Dinânica longitudinal (PC):
hlU = @(t)(y2d(t)/kNM1);
%Dinâmica completa (EF):
y2dAj = @(t)(Nr0bAj(y1d(t)));
hU = @(t)(y2dAj(t)/kNM1);
fhRestr = @(t)([hU(t); hRestr]);
rEstrT = rEstr + 1;

y2dAjMM = @(t)(Nr0bAjMM(y1d(t)));
hUmm = @(t)(y2dAjMM(t)/kNM1);
%Obs.: hU <-> (hU - dhU) no metamodelo

nUy1 = 0.01;
hUmmAst = @(t,yb1)(hU(t) - nUy1*(lambda1^4/kNM1)*(yb1 - y1d(t)));

%%
%Compensação do movimento longitudinal, devido ao deslocamento lateral:
%{
fv2t = @(yL,zL,ypL,zpL)([ypL, zpL]*fuTilHatfiL(yL,zL));
M2 = Mb2 + dm;
L2 = L - L1;
fuLEst = @(yL,zL,ypL,zpL)(fv2t(yL,zL,ypL,zpL)*sqrt(M2*((M2/L1 + M1/L2)*g)^(-1)));
fduNM1 = @(yL,zL,ypL,zpL)(L - (sqrt(L1^2 - fuLEst(yL,zL,ypL,zpL)^2) + sqrt(L2^2 - fuLEst(yL,zL,ypL,zpL)^2)));
duNM10 = uNM1aj0;
%}

%%
%Tentativas de metamodelo a partir do modelo em EF:
%{
SigmaReFmm = [IdNg(:,6*1-5),IdNg(:,6*1),IdNg(:,(6*2-5):(6*(N1+1)-5)),IdNg(:,6*(N1+1):(6*N)),IdNg(:,6*(N+1)-5),IdNg(:,6*(N+1))];
M1rMM = SigmaReFmm.'*M1gEF*SigmaReFmm;
KbrMM = SigmaReFmm.'*KgEFb*SigmaReFmm;
argumentos.SigmaReF = SigmaReFmm;
[ PtXmm, sigTiltXmm, PuMM, sigTiluMM, PvwMM, sigTilvwMM,...
    AbvtMM, AbvlMM, infoMM ] =...
    fModMMod( M1rMM, KbrMM, argumentos );
argumentos.SigmaReF = SigmaReF;
%}
%Metamodelo apenas com dinâmica de torção (linear):
%{
MtXref = PtX.'*(M1restr*PtX);
MtX = Id3;
CtX = PtX.'*(Crestr*PtX);
KtXref = PtX.'*(KbRestr*PtX); 
KtX = diag(sigTiltX);
BtX = PtX.'*Br(:,1);

AktX = [zeros(3), Id3; -KtX, -CtX];
BktX = [zeros(3,1); BtX];
PtXg = SigmaReF*PtX;
C1tX = [zeros(1,3), PtXg(6*(N+1),:)];
%}
%{
MtXref = PtXmm.'*(M1rMM*PtXmm);
MtX = Id3;
CrMM = SigmaReFmm.'*CgEF*SigmaReFmm;
CtX = PtXmm.'*(CrMM*PtXmm);
KtXref = PtXmm.'*(KbrMM*PtXmm); 
KtX = diag(sigTiltXmm);
Br = SigmaReFmm.'*Bef;
BtX = PtXmm.'*Br(:,1);

AktX = [zeros(3), Id3; -KtX, -CtX];
BktX = [zeros(3,1); BtX];
PtXg = SigmaReFmm*PtXmm;
C1tX = [zeros(1,3), PtXg(6*(N+1),:)];
%}
%
sigTilk = Abvl(1:jPu2,1);
Pk = Abvt(:,1:jPu2);
sPk = size(Pk);
cPk = sPk(1,2);
MtXref = Pk.'*(M1restr*Pk);
MtX = eye(cPk);
CtX = Pk.'*(Crestr*Pk);
KtXref = Pk.'*(KbRestr*Pk); 
KtX = diag(sigTilk);
BtX = Pk.'*Br(:,1);
Bl = Pk.'*Br(:,2);

%Contato com pinos de atrito (mu = 0):
fubNM1 = @(ubEF)(IdNg(Ng-5,:)*ubEF);
fKl = @(ubEF)(kNM1*(IdNg(:,Ng-5)*IdNg(:,Ng-5).')*double(fubNM1(ubEF)>=0));
fKlr = @(ubEFr)(SigmaReF.'*fKl(SigmaReF*ubEFr)*SigmaReF);
fKlrk = @(ubrk)(Pk.'*(fKlr(Pk*ubrk)*Pk));
fKtX = @(ubrk)(KtX + fKlrk(ubrk));

%Restriçao extrínseca:
GamaRrk = IdNg(1,:)*SigmaReF*Pk;
hRrk = @(t)(hU(t));

fAktX = @(ubrk)([zeros(cPk), eye(cPk); -fKtX(ubrk), -CtX]);
BktX = [zeros(cPk,1); BtX];
Bku = [zeros(cPk,1); Bl];
PkGlb = SigmaReF*Pk;
Cu1 = [PkGlb(6*1-5,:), zeros(1,cPk)];
CtpX1 = [zeros(1,cPk), PkGlb(6*1,:)];
CuN1M1 = [PkGlb(6*(N1+1)-5,:), zeros(1,cPk)];
CvN1M1 = [PkGlb(6*(N1+1)-4,:), zeros(1,cPk)];
CwN1M1 = [PkGlb(6*(N1+1)-2,:), zeros(1,cPk)];
CtpXN1M1 = [zeros(1,cPk), PkGlb(6*(N1+1),:)];
CuNM1 = [PkGlb(6*(N+1)-5,:), zeros(1,cPk)];
CtpXNM1 = [zeros(1,cPk), PkGlb(6*(N+1),:)];
%}
%{
jPu2 = infoMM.jPu2;
sigTilk = AbvlMM(1:jPu2,1);
Pk = AbvtMM(:,1:jPu2);
sPk = size(Pk);
cPk = sPk(1,2);
MtXref = Pk.'*(M1rMM*Pk);
MtX = eye(cPk);
CrMM = SigmaReFmm.'*CgEF*SigmaReFmm;
CtX = Pk.'*(CrMM*Pk);
KtXref = Pk.'*(KbrMM*Pk); 
KtX = diag(sigTilk);
Br = SigmaReFmm.'*Bef;
BtX = Pk.'*Br(:,1);
Bl = Pk.'*Br(:,2);

AktX = [zeros(cPk), eye(cPk); -KtX, -CtX];
BktX = [zeros(cPk,1); BtX];
Bku = [zeros(cPk,1); Bl];
PtXg = SigmaReFmm*Pk;
C1tX = [zeros(1,cPk), PtXg(6*(N+1),:)];
C2u = [PtXg(6*(N+1)-5,:), zeros(1,cPk)];
%}
%Forma normal do metamodelo:
%{
Kl = kNM1*(IdNg(:,Ng-5)*IdNg(:,Ng-5).');
Klr = SigmaReF.'*Kl*SigmaReF;
Klrk = Pk.'*(Klr*Pk);
KtXc = KtX + Klrk;
AktX = [zeros(cPk), eye(cPk); -KtXc, -CtX];

beta0tX = CtpXNM1;
beta1tX = beta0tX*AktX;
beta0BktX = beta0tX*BktX;
beta0Bku = beta0tX*Bku;
a11 = @(ZktX)(beta1tX*ZktX);
b111 = beta0BktX;
b121 = beta0Bku;

beta2tX = beta1tX*AktX;
beta1BktX = beta1tX*BktX;
beta1Bku = beta1tX*Bku;
a12 = @(ZktX)(beta2tX*ZktX);
b112 = beta1BktX;
b122 = beta1Bku;

beta3tX = beta2tX*AktX;
beta2BktX = beta2tX*BktX;
beta2Bku = beta2tX*Bku;
a13 = @(ZktX)(beta3tX*ZktX);
b113 = beta2BktX;
b123 = beta2Bku;

beta4tX = beta3tX*AktX;
beta3BktX = beta3tX*BktX;
beta3Bku = beta3tX*Bku;
a14 = @(ZktX)(beta4tX*ZktX);
b114 = beta3BktX;
b124 = beta3Bku;

beta5tX = beta4tX*AktX;
beta4BktX = beta4tX*BktX;
beta4Bku = beta4tX*Bku;
a15 = @(ZktX)(beta5tX*ZktX);
b115 = beta4BktX;
b125 = beta4Bku;

beta6tX = beta5tX*AktX;
beta5BktX = beta5tX*BktX;
beta5Bku = beta5tX*Bku;
a16 = @(ZktX)(beta6tX*ZktX);
b116 = beta5BktX;
b126 = beta5Bku;

beta7tX = beta6tX*AktX;
beta6BktX = beta6tX*BktX;
beta6Bku = beta6tX*Bku;
a17 = @(ZktX)(beta7tX*ZktX);
b117 = beta6BktX;
b127 = beta6Bku;

adTx = [b111, b121, b112, b122, b113, b123, b114, b124, b115, b125,...
    b116, b126, b117, b127];

alfa0u = CuNM1;
alfa1u = alfa0u*AktX;
alfa0BktX = alfa0u*BktX;
alfa0Bku = alfa0u*Bku;
a21 = @(Zku)(alfa1u*Zku);
b211 = alfa0BktX;
b221 = alfa0Bku;

alfa2u = alfa1u*AktX;
alfa1BktX = alfa1u*BktX;
alfa1Bku = alfa1u*Bku;
a22 = @(Zku)(alfa2u*Zku);
b212 = alfa1BktX;
b222 = alfa1Bku;

alfa3u = alfa2u*AktX;
alfa2BktX = alfa2u*BktX;
alfa2Bku = alfa2u*Bku;
a23 = @(Zku)(alfa3u*Zku);
b213 = alfa2BktX;
b223 = alfa2Bku;

alfa4u = alfa3u*AktX;
alfa3BktX = alfa3u*BktX;
alfa3Bku = alfa3u*Bku;
a24 = @(Zku)(alfa4u*Zku);
b214 = alfa3BktX;
b224 = alfa3Bku;

alfa5u = alfa4u*AktX;
alfa4BktX = alfa4u*BktX;
alfa4Bku = alfa4u*Bku;
a25 = @(Zku)(alfa5u*Zku);
b215 = alfa4BktX;
b225 = alfa4Bku;

alfa6u = alfa5u*AktX;
alfa5BktX = alfa5u*BktX;
alfa5Bku = alfa5u*Bku;
a26 = @(Zku)(alfa6u*Zku);
b216 = alfa5BktX;
b226 = alfa5Bku;

alfa7u = alfa6u*AktX;
alfa6BktX = alfa6u*BktX;
alfa6Bku = alfa6u*Bku;
a27 = @(Zku)(alfa7u*Zku);
b217 = alfa6BktX;
b227 = alfa6Bku;

adU = [b211, b221, b212, b222, b213, b223, b214, b224, b215, b225,...
    b216, b226, b217, b227];
%}

%%
%Newmark:
gamaNmk = argumentos.gama_newmark;
betaNmk = argumentos.beta_newmark;
f_up_nM1 = @(u_nM1, u_n, up_n, upp_n, dt)((gamaNmk/(betaNmk*dt))*(u_nM1 - u_n) +...
    (1 - gamaNmk/betaNmk)*up_n + (1 - gamaNmk/(2*betaNmk))*dt*upp_n);
repositorio.f_up_nM1 = f_up_nM1;
f_upp_nM1 = @(u_nM1, u_n, up_n, upp_n, dt)((1/(betaNmk*dt^2))*(u_nM1 - u_n) -...
    (1/(betaNmk*dt))*up_n + (1 - 1/(2*betaNmk))*upp_n);
repositorio.f_upp_nM1 = f_upp_nM1;

%%
%Pré-integração:
alfa1Nmk = @(dtnM1)(gamaNmk*(betaNmk*dtnM1)^-1);
alfa2Nmk = @(dtnM1)((betaNmk*dtnM1^2)^-1);
fi1Nmk = @(dtnM1,un,upn,uppn)(-alfa1Nmk(dtnM1)*un +...
    (1 - gamaNmk*betaNmk^-1)*upn +...
    (1 - gamaNmk*(2*betaNmk)^-1)*dtnM1*uppn);
repositorio.fi1Nmk = fi1Nmk;
fi2Nmk = @(dtnM1,un,upn,uppn)(-alfa2Nmk(dtnM1)*un -...
    (betaNmk*dtnM1)^-1*upn +...
    (1 - (2*betaNmk)^-1)*uppn);
repositorio.fi2Nmk = fi2Nmk;
fupnM1 = @(dtnM1,unM1,un,upn,uppn)(alfa1Nmk(dtnM1)*unM1 +...
    fi1Nmk(dtnM1,un,upn,uppn));
repositorio.fupnM1 = fupnM1;
fuppnM1 = @(dtnM1,unM1,un,upn,uppn)(alfa2Nmk(dtnM1)*unM1 +...
    fi2Nmk(dtnM1,un,upn,uppn));
repositorio.fuppnM1 = fuppnM1;
fdupduTnmk = @(dtnM1)(alfa1Nmk(dtnM1));
repositorio.fdupduTnmk = fdupduTnmk;
fduppduTnmk = @(dtnM1)(alfa2Nmk(dtnM1));
repositorio.fduppduTnmk = fduppduTnmk;

Id4 = eye(4);
switch simulacao
    case 1
        %Entrada:
        Tref = 2; %Nm
        fTx = @(t)(Tref);
        %Newmark:
        KtNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mt + alfa1Nmk(dtnM1)*Ct + Kt);
    case 2
        %Entrada:
        Kf = 10^5; %N/m
        ub1Ref = 1*10^-4; %m
        fFbu = @(t,ub1)(-Kf*(ub1 - ub1Ref));
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Kbu);
    case 3
        %Obs.: Passa a passo natural:
        %(1) Passar para a variável x
        %(2) Aplicar fórmulas de Newmark
        
        Amix = Id2(:,1);
        Afix = Id2(:,2);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
    case 4
        dfipBddfiBdnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        GamaFiBdUnM1 = @(dtnM1,fiBdnM1,fiBdn,fipBdn,fippBdn)(...
            GamaFiBdU(fupnM1(dtnM1,fiBdnM1,fiBdn,fipBdn,fippBdn)));
        dGamaFiBdUdfiBdUnM1 = @(dtnM1,fiBdnM1,fiBdn,fipBdn,fippBdn)(...
            dGamaFiUdfipBdU(fupnM1(dtnM1,fiBdnM1,fiBdn,fipBdn,fippBdn))*...
            dfipBddfiBdnM1(dtnM1));
    case 5
        Id5 = eye(5);
        Atx = Id5(:,1:3);
        Amix = Id5(:,4);
        Afix = Id5(:,5);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        
        %Newmark:
        KtNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mt + alfa1Nmk(dtnM1)*Ct + Kt);
        KtNmkx = @(dtnM1)(KtNmk(dtnM1)*Atx.');
        Mtx = Mt*Atx.';
        Ctx = Ct*Atx.';
        
        %Camada limite e dinâmica interna:
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil)(fXt(Atx.'*xTil,Atx.'*xpTil));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil)(s1(t,fXt(uTil,upTil)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd)(...
            b1fn\(UcThat(t,Z,X) - KbLcT(t,Z,X,mid,fiBd)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd)(...
            UcT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        UcTx = @(t,xTil,xpTil)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = -da1fndZT;
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z)(dda1fndZT*sign(da1fn(Z)));
        dFlc1dZT = @(Z)(dabsda1fndZT(Z));
        dabsUcThatdZT = @(t,Z,X)(dUcThatdZT*sign(UcThat(t,Z,X)));
        dKlcTdZT = @(t,Z,X)((1 - D11fn)\(dFlc1dZT(Z) + D11fn*dabsUcThatdZT(t,Z,X)));
        dabsUcThatdXT = @(t,Z,X)(dUcThatdXT*sign(UcThat(t,Z,X)));
        dKlcTdXT = @(t,Z,X)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X)));
        dKbLcTdZT = @(t,Z,X)(dKlcTdZT(t,Z,X));
        dKbLcTdXT = @(t,Z,X)(dKlcTdXT(t,Z,X));
        dUcTuTildZT = @(t,Z,X,fiBd)(...
            b1fn\(dUcThatdZT - dKbLcTdZT(t,Z,X)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X)*sat(s1(t,X)/fiBd) +...
             KbLcT(t,Z,X,mid,fiBd)*dsatS1FibddXT(t,X,fiBd))));
        dXtdZtT = sigDift;
        DUcTuTilDZT = @(t,Z,X,mid,fiBd)(dUcTuTildZT(t,Z,X,fiBd) +...
            dUcTuTildXT(t,Z,X,mid,fiBd)*dXtdZtT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZdupTilT);
        dZddmid = sigDift\dXtddmid;
        dUcTdHatdmid = @(t,mid)(-da1fndZT*dZddmid);
        dabsUcTdHatdmid = @(t,mid)(dUcTdHatdmid(t,mid)*sign(UcTdHat(t,mid)));
        dKlcTddmid = @(t,mid)((1 - D11fn)\(dFlc1dZT(fZtd(t,mid))*dZddmid + D11fn*dabsUcTdHatdmid(t,mid)));
        dKbLcTdmid = @(t,mid)(-dKlcTddmid(t,mid));
        dUcTuTildmid = @(t,X,mid,fiBd)(b1fn\(-dKbLcTdmid(t,mid)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd)(b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcT(t,Z,X,mid,fiBd)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd)(dUcTuTilduTilT(t,Z,X,mid,fiBd)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxTx = @(t,xTil,xpTil)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Atx.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd)(dUcTuTildupTilT(t,Z,X,mid,fiBd)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxpTx = @(t,xTil,xpTil)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
    case 6
        %Obs.: xTil = Aux*uTil + Afix*fiBdU
        Aux = Id4(:,1:3);
        AfiUx = Id4(:,4);
        fux = @(xTil)(Aux.'*xTil);
        ffiBdUx = @(xTil)(AfiUx.'*xTil);
        
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Kbu);
        KuNmkx = @(dtnM1)(KuNmk(dtnM1)*Aux.');
        Mux = Mu*Aux.';
        Cux = Cu*Aux.';
        
        %Camada limite:
        ffipBdUx = @(xTilp)(AfiUx.'*xTilp);
        GamaFiBdUx = @(xTilp)(GamaFiBdU(ffipBdUx(xTilp)));
        GamaFiBdUxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdUx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiUx = KfiU*AfiUx.';
        
        KlcUdx = @(t)(KlcUd(t));
        KlcUdxnM1 = @(tnM1)(KlcUdx(tnM1));
        
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        
        dfipBddxTilpT = AfiUx.';
        dGamaFiUdfipBdUx = @(xTilp)(dGamaFiUdfipBdU(ffipBdUx(xTilp)));
        dGamaFiUxdxTilpT = @(xTilp)(dGamaFiUdfipBdUx(xTilp)*dfipBddxTilpT);
        dGamaFiUxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiUxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZudx = @(t)(fZud(t));
        
        %Vetores de estado:
        fZux = @(xTil,xpTil)(fZ(Aux.'*xTil,Aux.'*xpTil));
        fZuxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZux(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXux = @(xTil,xpTil)(fXu(Aux.'*xTil,Aux.'*xpTil));
        fXuxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXux(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s2uTil = @(t,uTil,upTil)(s2(t,fXu(uTil,upTil)));
        s2xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s2(tnM1,fXuxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcU = @(t,Z,X,fiBd)(...
            b2fn\(UcUhat(t,Z,X) - KbLcU(t,Z,X,fiBd)*sat(s2(t,X)/fiBd)));
        UcUuTil = @(t,uTil,upTil,fiBd)(...
            UcU(t,fZ(uTil,upTil),fXu(uTil,upTil),fiBd));
        UcUx = @(t,xTil,xpTil)(...
            UcUuTil(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        UcUxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcUx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTildxT = Aux.';
        dfiBdUdxT = AfiUx.';
        dUcUhatdZT = -da2fndZT;
        dy2dXT = Id6(1,:);
        dy2pdXT = Id6(2,:);
        dy2ppdXT = Id6(3,:);
        dy2pppdXT = Id6(4,:);
        dy2ppppdXT = Id6(5,:);
        dy2pppppdXT = Id6(6,:);
        dy6r2dXT = -5*lambda2*dy2pppppdXT -...
            10*lambda2^2*dy2ppppdXT -...
            10*lambda2^3*dy2pppdXT -...
            5*lambda2^4*dy2ppdXT -...
            lambda2^5*dy2pdXT;
        dUcUhatdXT = dy6r2dXT;
        dabsda2fndZT = @(Z)(dda2fndZT*sign(da2fn(Z)));
        dFlc2dZT = @(Z)(dabsda2fndZT(Z));
        dabsUcUhatdZT = @(t,Z,X)(dUcUhatdZT*sign(UcUhat(t,Z,X)));
        dKlcUdZT = @(t,Z,X)((1 - D22fn)\(dFlc2dZT(Z) + D22fn*dabsUcUhatdZT(t,Z,X)));
        dabsUcUhatdXT = @(t,Z,X)(dUcUhatdXT*sign(UcUhat(t,Z,X)));
        dKlcUdXT = @(t,Z,X)((1 - D22fn)\(D22fn*dabsUcUhatdXT(t,Z,X)));
        dKbLcUdZT = @(t,Z,X)(dKlcUdZT(t,Z,X));
        dKbLcUdXT = @(t,Z,X)(dKlcUdXT(t,Z,X));
        dUcUuTildZT = @(t,Z,X,fiBd)(...
            b2fn\(dUcUhatdZT - dKbLcUdZT(t,Z,X)*sat(s2(t,X)/fiBd)));
        dy5r2dXT = -5*lambda2*dy2ppppdXT -...
            10*lambda2^2*dy2pppdXT -...
            10*lambda2^3*dy2ppdXT -...
            5*lambda2^4*dy2pdXT -...
            lambda2^5*dy2dXT;
        ds2dXT = dy2pppppdXT - dy5r2dXT;
        dsatS2FibddXT = @(t,X,fiBd)(dsatdnu(s2(t,X)/fiBd)*(ds2dXT/fiBd));
        dUcUuTildXT = @(t,Z,X,fiBd)(...
            b2fn\(dUcUhatdXT -...
            (dKbLcUdXT(t,Z,X)*sat(s2(t,X)/fiBd) +...
             KbLcU(t,Z,X,fiBd)*dsatS2FibddXT(t,X,fiBd))));
        dXudZuT = sigDifu;
        DUcUuTilDZT = @(t,Z,X,fiBd)(dUcUuTildZT(t,Z,X,fiBd) +...
            dUcUuTildXT(t,Z,X,fiBd)*dXudZuT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcUuTilduTilT = @(t,Z,X,fiBd)(DUcUuTilDZT(t,Z,X,fiBd)*dZduTilT);
        dUcUuTildupTilT = @(t,Z,X,fiBd)(DUcUuTilDZT(t,Z,X,fiBd)*dZdupTilT);
        dKbLcUdfiBd = KfiU;
        dsatSFibdUdfiBd = @(t,X,fiBd)(dsatdnu(s2(t,X)/fiBd)*(-s2(t,X)/fiBd^2));
        dUcUuTildfiBd = @(t,Z,X,fiBd)(b2fn\(-(dKbLcUdfiBd*sat(s2(t,X)/fiBd) +...
            KbLcU(t,Z,X,fiBd)*dsatSFibdUdfiBd(t,X,fiBd))));
        dUcUxdxT = @(t,Z,X,fiBd)(dUcUuTilduTilT(t,Z,X,fiBd)*duTildxT +...
            dUcUuTildfiBd(t,Z,X,fiBd)*dfiBdUdxT);
        dUcUxdxTuTil = @(t,uTil,upTil,fiBd)(...
            dUcUxdxT(t,fZ(uTil,upTil),fXu(uTil,upTil),fiBd));
        dUcUxdxTx = @(t,xTil,xpTil)(...
            dUcUxdxTuTil(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        dUcUxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Aux.';
        dUcUxdxpT = @(t,Z,X,fiBd)(dUcUuTildupTilT(t,Z,X,fiBd)*dupTildxpT);
        dUcUxdxpTuTil = @(t,uTil,upTil,fiBd)(...
            dUcUxdxpT(t,fZ(uTil,upTil),fXu(uTil,upTil),fiBd));
        dUcUxdxpTx = @(t,xTil,xpTil)(...
            dUcUxdxpTuTil(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        dUcUxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcUxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcUxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
    case 7
        %Obs.: xTil = Aux*uTil + Afix*fiBdU
        Id4 = eye(4);
        Aux = Id4(:,1:3);
        AfiUx = Id4(:,4);
        fux = @(xTil)(Aux.'*xTil);
        ffiBdUx = @(xTil)(AfiUx.'*xTil);
        
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Kbu);
        KuNmkx = @(dtnM1)(KuNmk(dtnM1)*Aux.');
        Mux = Mu*Aux.';
        Cux = Cu*Aux.';
        
        %Camada limite e dinâmica interna:
        ffipBdUx = @(xTilp)(AfiUx.'*xTilp);
        GamaFiBdUxPS = @(xTilp)(GamaFiBdUps(ffipBdUx(xTilp)));
        GamaFiBdUxnM1PS = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdUxPS(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiUxPS = KfiUps*AfiUx.';
        
        KlcUdxPS = @(tal)(KlcUdPS(tal));
        KlcUdxnM1PS = @(talnM1)(KlcUdxPS(talnM1));
        
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        
        dfipBddxTilpT = AfiUx.';
        dGamaFiUdfipBdUxPS = @(xTilp)(dGamaFiUdfipBdUps(ffipBdUx(xTilp)));
        dGamaFiUxdxTilpTps = @(xTilp)(dGamaFiUdfipBdUxPS(xTilp)*dfipBddxTilpT);
        dGamaFiUxdxTilnM1Tps = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiUxdxTilpTps(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZudxPS = @(t)(fZudPS(t));
        
        %Vetores de estado:
        fZux = @(xTil,xpTil)(fZ(Aux.'*xTil,Aux.'*xpTil));
        fZuxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZux(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXuxPS = @(xTil,xpTil)(fXuPS(Aux.'*xTil,Aux.'*xpTil));
        fXuxnM1PS = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXuxPS(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s2uTilPS = @(t,uTil,upTil)(s2PS(t,fXuPS(uTil,upTil)));
        s2xnM1PS = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s2PS(tnM1,fXuxnM1PS(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcUps = @(t,Z,X,fiBd)(...
            b2fn\(UcUhatPS(t,Z,X) - KbLcUps(t,Z,X,fiBd)*sat(s2PS(t,X)/fiBd)));
        UcUuTilPS = @(t,uTil,upTil,fiBd)(...
            UcUps(t,fZ(uTil,upTil),fXuPS(uTil,upTil),fiBd));
        UcUxPS = @(t,xTil,xpTil)(...
            UcUuTilPS(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        UcUxnM1PS = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcUxPS(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTildxT = Aux.';
        dfiBdUdxT = AfiUx.';
        dUcUhatdZTps = -da2fndZTps;
        dy2dXT = Id6(1,:);
        dy2pdXT = Id6(2,:);
        dy2ppdXT = Id6(3,:);
        dy2pppdXT = Id6(4,:);
        dy2ppppdXT = Id6(5,:);
        dy2pppppdXT = Id6(6,:);
        dy6r2dXT = -5*lambda2*dy2pppppdXT -...
            10*lambda2^2*dy2ppppdXT -...
            10*lambda2^3*dy2pppdXT -...
            5*lambda2^4*dy2ppdXT -...
            lambda2^5*dy2pdXT;
        dUcUhatdXTps = dy6r2dXT;
        dabsda2fndZTps = @(Z)(dda2fndZTps*sign(da2fnPS(Z)));
        dFlc2dZTps = @(Z)(dabsda2fndZTps(Z));
        dabsUcUhatdZTps = @(t,Z,X)(dUcUhatdZTps*sign(UcUhatPS(t,Z,X)));
        dKlcUdZTps = @(t,Z,X)((1 - D22fnPS)\(dFlc2dZTps(Z) + D22fnPS*dabsUcUhatdZTps(t,Z,X)));
        dabsUcUhatdXTps = @(t,Z,X)(dUcUhatdXTps*sign(UcUhatPS(t,Z,X)));
        dKlcUdXTps = @(t,Z,X)((1 - D22fnPS)\(D22fnPS*dabsUcUhatdXTps(t,Z,X)));
        dKbLcUdZTps = @(t,Z,X)(dKlcUdZTps(t,Z,X));
        dKbLcUdXTps = @(t,Z,X)(dKlcUdXTps(t,Z,X));
        dUcUuTildZTps = @(t,Z,X,fiBd)(...
            b2fnPS\(dUcUhatdZTps - dKbLcUdZTps(t,Z,X)*sat(s2PS(t,X)/fiBd)));
        dy5r2dXTps = -5*lambda2*dy2ppppdXT -...
            10*lambda2^2*dy2pppdXT -...
            10*lambda2^3*dy2ppdXT -...
            5*lambda2^4*dy2pdXT -...
            lambda2^5*dy2dXT;
        ds2dXTps = dy2pppppdXT - dy5r2dXTps;
        dsatS2FibddXTps = @(t,X,fiBd)(dsatdnu(s2PS(t,X)/fiBd)*(ds2dXTps/fiBd));
        dUcUuTildXTps = @(t,Z,X,fiBd)(...
            b2fnPS\(dUcUhatdXTps -...
            (dKbLcUdXTps(t,Z,X)*sat(s2PS(t,X)/fiBd) +...
             KbLcUps(t,Z,X,fiBd)*dsatS2FibddXTps(t,X,fiBd))));
        dXudZuTps = sigDifuPS;
        DUcUuTilDZTps = @(t,Z,X,fiBd)(dUcUuTildZTps(t,Z,X,fiBd) +...
            dUcUuTildXTps(t,Z,X,fiBd)*dXudZuTps);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcUuTilduTilTps = @(t,Z,X,fiBd)(DUcUuTilDZTps(t,Z,X,fiBd)*dZduTilT);
        dUcUuTildupTilTps = @(t,Z,X,fiBd)(DUcUuTilDZTps(t,Z,X,fiBd)*dZdupTilT);
        dKbLcUdfiBdPS = KfiUps;
        dsatSFibdUdfiBdPS = @(t,X,fiBd)(dsatdnu(s2PS(t,X)/fiBd)*(-s2PS(t,X)/fiBd^2));
        dUcUuTildfiBdPS = @(t,Z,X,fiBd)(b2fnPS\(-(dKbLcUdfiBdPS*sat(s2PS(t,X)/fiBd) +...
            KbLcUps(t,Z,X,fiBd)*dsatSFibdUdfiBdPS(t,X,fiBd))));
        dUcUxdxTps = @(t,Z,X,fiBd)(dUcUuTilduTilTps(t,Z,X,fiBd)*duTildxT +...
            dUcUuTildfiBdPS(t,Z,X,fiBd)*dfiBdUdxT);
        dUcUxdxTuTilPS = @(t,uTil,upTil,fiBd)(...
            dUcUxdxTps(t,fZ(uTil,upTil),fXuPS(uTil,upTil),fiBd));
        dUcUxdxTxPS = @(t,xTil,xpTil)(...
            dUcUxdxTuTilPS(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        dUcUxdxTnM1PS = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUxdxTxPS(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Aux.';
        dUcUxdxpTps = @(t,Z,X,fiBd)(dUcUuTildupTilTps(t,Z,X,fiBd)*dupTildxpT);
        dUcUxdxpTuTilPS = @(t,uTil,upTil,fiBd)(...
            dUcUxdxpTps(t,fZ(uTil,upTil),fXuPS(uTil,upTil),fiBd));
        dUcUxdxpTxPS = @(t,xTil,xpTil)(...
            dUcUxdxpTuTilPS(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        dUcUxdxpTnM1PS = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUxdxpTxPS(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcUxDxnM1Tps = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUxdxTnM1PS(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcUxdxpTnM1PS(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
    case 8
        %%
        %Obs.: xTil = Aux*uTil + Afix*fiBdU
        Id4 = eye(4);
        Aux = Id4(:,1:3);
        AfiUx = Id4(:,4);
        fux = @(xTil)(Aux.'*xTil);
        ffiBdUx = @(xTil)(AfiUx.'*xTil);
        
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Kbu);
        KuNmkx = @(dtnM1)(KuNmk(dtnM1)*Aux.');
        Mux = Mu*Aux.';
        Cux = Cu*Aux.';
        
        %Camada limite e dinâmica interna:
        ffipBdUx = @(xTilp)(AfiUx.'*xTilp);
        GamaFiBdUx = @(xTilp)(GamaFiBdU(ffipBdUx(xTilp)));
        GamaFiBdUxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdUx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiUx = KfiU*AfiUx.';
        
        KlcUdx = @(t)(KlcUd(t));
        KlcUdxnM1 = @(tnM1)(KlcUdx(tnM1));
        
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        
        dfipBddxTilpT = AfiUx.';
        dGamaFiUdfipBdUx = @(xTilp)(dGamaFiUdfipBdU(ffipBdUx(xTilp)));
        dGamaFiUxdxTilpT = @(xTilp)(dGamaFiUdfipBdUx(xTilp)*dfipBddxTilpT);
        dGamaFiUxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiUxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZudx = @(t)(fZud(t));
        
        %Vetores de estado:
        fZux = @(xTil,xpTil)(fZ(Aux.'*xTil,Aux.'*xpTil));
        fZuxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZux(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXuxPS = @(xTil,xpTil)(fXuPS(Aux.'*xTil,Aux.'*xpTil));
        fXuxnM1PS = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXuxPS(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s2uTilPS = @(tal,uTil,upTil)(s2PS(tal,fXuPS(uTil,upTil)));
        s2uTil = @(t,uTil,upTil)(s2uTilPS(t/epsPS,uTil,upTil));
        s2xnM1PS = @(talnM1,dtalnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s2PS(talnM1,fXuxnM1PS(dtalnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        s2xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s2xnM1PS(tnM1/epsPS,dtnM1/epsPS,xTilnM1,xTiln,xTilpn,xTilppn));
        
        %Lei de controle:
        FufHat = @(tal,Zuf,Xuf)(UcUhatPS(tal,Zuf,Xuf));
        Fuf = @(tal,Zuf,Xuf,fiBd)(b2fnPS\(FufHat(tal,Zuf,Xuf) -...
            KbLcUps(tal,Zuf,Xuf,fiBd)*sat(s2PS(tal,Xuf)/fiBd)));
        fZuf = @(t,Zu)(Zu + Au\(BbU*y2d(t)));
        fXuf = @(t,Zu)(sigDifuPS*fZuf(t,Zu));
        UcU = @(t,Zu,fiBd)(y2d(t) + Fuf(t/epsPS,fZuf(t,Zu),fXuf(t,Zu),fiBd));
        UcUuTil = @(t,uTil,upTil,fiBd)(UcU(t,fZ(uTil,upTil),fiBd));
        UcUx = @(t,xTil,xpTil)(UcUuTil(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        UcUxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcUx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Derivadas:
        duTildxT = Aux.';
        dfiBdUdxT = AfiUx.';
        dFufHatdZfT = -da2fndZTps;
        dy2dXT = Id6(1,:);
        dy2pdXT = Id6(2,:);
        dy2ppdXT = Id6(3,:);
        dy2pppdXT = Id6(4,:);
        dy2ppppdXT = Id6(5,:);
        dy2pppppdXT = Id6(6,:);
        dy6r2dXT = -5*lambda2PS*dy2pppppdXT -...
            10*lambda2PS^2*dy2ppppdXT -...
            10*lambda2PS^3*dy2pppdXT -...
            5*lambda2PS^4*dy2ppdXT -...
            lambda2PS^5*dy2pdXT;
        dFufHatdXfT = dy6r2dXT;
        dabsda2fndZfT = @(Zf)(dda2fndZTps*sign(da2fnPS(Zf)));
        dFlc2dZfT = @(Zf)(dabsda2fndZfT(Zf));
        dabsFufHatdZfT = @(tal,Zf,Xf)(dFufHatdZfT*sign(FufHat(tal,Zf,Xf)));
        dKlcUdZfT = @(tal,Zf,Xf)((1 - D22fnPS)\(dFlc2dZfT(Zf) + D22fnPS*dabsFufHatdZfT(tal,Zf,Xf)));
        dabsFufHatdXfT = @(tal,Zf,Xf)(dFufHatdXfT*sign(FufHat(tal,Zf,Xf)));
        dKlcUdXfT = @(tal,Zf,Xf)((1 - D22fnPS)\(D22fnPS*dabsFufHatdXfT(tal,Zf,Xf)));
        dKbLcUdZfT = @(tal,Zf,Xf)(dKlcUdZfT(tal,Zf,Xf));
        dKbLcUdXfT = @(tal,Zf,Xf)(dKlcUdXfT(tal,Zf,Xf));
        dFufuTildZfT = @(tal,Zf,Xf,fiBd)(...
            b2fnPS\(dFufHatdZfT - dKbLcUdZfT(tal,Zf,Xf)*sat(s2PS(tal,Xf)/fiBd)));
        dy5r2dXfT = -5*lambda2PS*dy2ppppdXT -...
            10*lambda2PS^2*dy2pppdXT -...
            10*lambda2PS^3*dy2ppdXT -...
            5*lambda2PS^4*dy2pdXT -...
            lambda2PS^5*dy2dXT;
        ds2dXfT = dy2pppppdXT - dy5r2dXfT;
        dsatS2FibddXfT = @(tal,Xf,fiBd)(dsatdnu(s2PS(tal,Xf)/fiBd)*(ds2dXfT/fiBd));
        dFufuTildXfT = @(tal,Zf,Xf,fiBd)(...
            b2fnPS\(dFufHatdXfT -...
            (dKbLcUdXfT(tal,Zf,Xf)*sat(s2PS(tal,Xf)/fiBd) +...
             KbLcUps(tal,Zf,Xf,fiBd)*dsatS2FibddXfT(tal,Xf,fiBd))));
        dXudZuTps = sigDifuPS;
        DFufuTilDZfT = @(tal,Zf,Xf,fiBd)(dFufuTildZfT(tal,Zf,Xf,fiBd) +...
            dFufuTildXfT(tal,Zf,Xf,fiBd)*dXudZuTps);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dKbLcUdfiBdPS = KfiUps;
        dsatSFibdUdfiBdPS = @(tal,Xf,fiBd)(dsatdnu(s2PS(tal,Xf)/fiBd)*(-s2PS(tal,Xf)/fiBd^2));
        dFufuTildfiBdPS = @(tal,Zf,Xf,fiBd)(b2fnPS\(-(dKbLcUdfiBdPS*sat(s2PS(tal,Xf)/fiBd) +...
            KbLcUps(tal,Zf,Xf,fiBd)*dsatSFibdUdfiBdPS(tal,Xf,fiBd))));
        
        dZfdZuT = Id6;
        dUcUduTilT = @(t,Zu,fiBd)(DFufuTilDZfT(t/epsPS,fZuf(t,Zu),fXuf(t,Zu),fiBd)*...
            dZfdZuT*dZduTilT);
        dUcUdupTilT = @(t,Zu,fiBd)(DFufuTilDZfT(t/epsPS,fZuf(t,Zu),fXuf(t,Zu),fiBd)*...
            dZfdZuT*dZdupTilT);
        dUcUdfiBd = @(t,Zu,fiBd)(dFufuTildfiBdPS(t/epsPS,fZuf(t,Zu),fXuf(t,Zu),fiBd));
        dUcUduTuTil = @(t,uTil,upTil,fiBd)(dUcUduTilT(t,fZ(uTil,upTil),fiBd));
        dUcUdupTuTil = @(t,uTil,upTil,fiBd)(dUcUdupTilT(t,fZ(uTil,upTil),fiBd));
        dUcUdfiBduTil = @(t,uTil,upTil,fiBd)(dUcUdfiBd(t,fZ(uTil,upTil),fiBd));
        dUcUduTx = @(t,xTil,xpTil)(dUcUduTuTil(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        dUcUdupTx = @(t,xTil,xpTil)(dUcUdupTuTil(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        dUcUdfiBdx = @(t,xTil,xpTil)(dUcUdfiBduTil(t,Aux.'*xTil,Aux.'*xpTil,AfiUx.'*xTil));
        dUcUdxT = @(t,xTil,xpTil)(dUcUduTx(t,xTil,xpTil)*duTildxT +...
            dUcUdfiBdx(t,xTil,xpTil)*dfiBdUdxT);
        dupTildxpT = Aux.';
        dUcUdxpT = @(t,xTil,xpTil)(dUcUdupTx(t,xTil,xpTil)*dupTildxpT);
        dUcUdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUdxT(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dUcUdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUdxpT(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcUxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcUdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcUdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
    case 9
        %Entrada:
        fFbu = @(t)(y2d(t));
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
    case 10
        %%
        %Dinâmica longitudinal
        %Entrada:
        fFbu = @(t)(y2d(t));
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        %%
        %Dinâmica de torção
        Id5 = eye(5);
        Atx = Id5(:,1:3);
        Amix = Id5(:,4);
        Afix = Id5(:,5);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        
        %Newmark:
        KtNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mt + alfa1Nmk(dtnM1)*Ct + Kt);
        KtNmkx = @(dtnM1)(KtNmk(dtnM1)*Atx.');
        Mtx = Mt*Atx.';
        Ctx = Ct*Atx.';
        
        %Camada limite e dinâmica interna:
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil)(fXt(Atx.'*xTil,Atx.'*xpTil));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil)(s1(t,fXt(uTil,upTil)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd)(...
            b1fn\(UcThat(t,Z,X) - KbLcT(t,Z,X,mid,fiBd)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd)(...
            UcT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        UcTx = @(t,xTil,xpTil)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = -da1fndZT;
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z)(dda1fndZT*sign(da1fn(Z)));
        dFlc1dZT = @(Z)(dabsda1fndZT(Z));
        dabsUcThatdZT = @(t,Z,X)(dUcThatdZT*sign(UcThat(t,Z,X)));
        dKlcTdZT = @(t,Z,X)((1 - D11fn)\(dFlc1dZT(Z) + D11fn*dabsUcThatdZT(t,Z,X)));
        dabsUcThatdXT = @(t,Z,X)(dUcThatdXT*sign(UcThat(t,Z,X)));
        dKlcTdXT = @(t,Z,X)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X)));
        dKbLcTdZT = @(t,Z,X)(dKlcTdZT(t,Z,X));
        dKbLcTdXT = @(t,Z,X)(dKlcTdXT(t,Z,X));
        dUcTuTildZT = @(t,Z,X,fiBd)(...
            b1fn\(dUcThatdZT - dKbLcTdZT(t,Z,X)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X)*sat(s1(t,X)/fiBd) +...
             KbLcT(t,Z,X,mid,fiBd)*dsatS1FibddXT(t,X,fiBd))));
        dXtdZtT = sigDift;
        DUcTuTilDZT = @(t,Z,X,mid,fiBd)(dUcTuTildZT(t,Z,X,fiBd) +...
            dUcTuTildXT(t,Z,X,mid,fiBd)*dXtdZtT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = sigDift\dXddmid;
        dUcTdHatdmid = @(t,mid)(-da1fndZT*dZddmid);
        dabsUcTdHatdmid = @(t,mid)(dUcTdHatdmid(t,mid)*sign(UcTdHat(t,mid)));
        dKlcTddmid = @(t,mid)((1 - D11fn)\(dFlc1dZT(fZtd(t,mid))*dZddmid + D11fn*dabsUcTdHatdmid(t,mid)));
        dKbLcTdmid = @(t,mid)(-dKlcTddmid(t,mid));
        dUcTuTildmid = @(t,X,mid,fiBd)(b1fn\(-dKbLcTdmid(t,mid)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd)(b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcT(t,Z,X,mid,fiBd)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd)(dUcTuTilduTilT(t,Z,X,mid,fiBd)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxTx = @(t,xTil,xpTil)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Atx.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd)(dUcTuTildupTilT(t,Z,X,mid,fiBd)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxpTx = @(t,xTil,xpTil)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
    case 11
        %%
        %Dinâmica longitudinal
        %Entrada:
        fFbu = @(t)(y2d(t));
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        %%
        %Dinâmica de torção
        Id5 = eye(5);
        Atx = Id5(:,1:3);
        Amix = Id5(:,4);
        Afix = Id5(:,5);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        
        %Newmark:
        KtNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mt + alfa1Nmk(dtnM1)*Ct + Kt);
        KtNmkx = @(dtnM1)(KtNmk(dtnM1)*Atx.');
        Mtx = Mt*Atx.';
        Ctx = Ct*Atx.';
        FatTcx = @(xTilp,y2)(FatTc(ftx(xTilp),y2));
        FatTcxnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            FatTcx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1));
        FatTdcx = @(xTilp,y2)(FatTdc(ftx(xTilp),y2));
        FatTdcxnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            FatTdcx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1));
        
        %Camada limite e dinâmica interna:
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxTilpT = Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxTilpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT = Atx.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
    case 12
        %Dinâmica longitudinal
        %Entrada:
        fFbu = @(t)(y2d(t));
        
        %Obs.: xTil = AuAcx*uTilAc + Amix*mid + Afix*fiBd
        Id10 = eye(10);
        AuAcx = Id10(:,1:8);
        Amix = Id10(:,9);
        Afix = Id10(:,10);
        %Obs.: uTilAc = Atu*uTeta + Auu*ubU
        Atu = Id8(:,[2,6,8]);
        Auu = Id8(:,[1,3,7]);
        %Obs.: Uc = AtUc*Tx + AuUc*Fub
        Id2 = eye(2);
        AtUc = Id2(:,1);
        AuUc = Id2(:,2);
        
        fubTilAc = @(xTil)(AuAcx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTuAc = @(ubTilAc)(Atu.'*ubTilAc);
        fuUuAc = @(ubTilAc)(Auu.'*ubTilAc);
        fuTxTil = @(xTil)(fuTuAc(fubTilAc(xTil)));
        fuUxTil = @(xTil)(fuUuAc(fubTilAc(xTil)));
        
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1ac + alfa1Nmk(dtnM1)*Cac + Kac);
        Knmkx = @(dtnM1)(Knmk(dtnM1)*AuAcx.');
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fTilc(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfTilcdubT(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfTilcdubpT(fubTilAc(xTil),fubTilAc(xTilp)));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fTildc(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfTildcdubT(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfTildcdubpT(fubTilAc(xTil),fubTilAc(xTilp)));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuAcx.';
        dfTildubppTx = @(xTil)(dfTildubppT(fubTilAc(xTil)));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        M1acx = M1ac*AuAcx.';
        Cacx = Cac*AuAcx.';
        
        %Camada limite e dinâmica interna:
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(fuTxTil(xTil),fuTxTil(xpTil)));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil)(fXt(fuTxTil(xTil),fuTxTil(xpTil)));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTeta,upTeta)(s1(t,fXt(uTeta,upTeta)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd)(...
            b1fn\(UcThat(t,Z,X) - KbLcT(t,Z,X,mid,fiBd)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd)(...
            UcT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        UcTx = @(t,xTil,xpTil)(...
            UcTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTetadxT = Atu.'*AuAcx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = -da1fndZT;
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z)(dda1fndZT*sign(da1fn(Z)));
        dFlc1dZT = @(Z)(dabsda1fndZT(Z));
        dabsUcThatdZT = @(t,Z,X)(dUcThatdZT*sign(UcThat(t,Z,X)));
        dKlcTdZT = @(t,Z,X)((1 - D11fn)\(dFlc1dZT(Z) + D11fn*dabsUcThatdZT(t,Z,X)));
        dabsUcThatdXT = @(t,Z,X)(dUcThatdXT*sign(UcThat(t,Z,X)));
        dKlcTdXT = @(t,Z,X)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X)));
        dKbLcTdZT = @(t,Z,X)(dKlcTdZT(t,Z,X));
        dKbLcTdXT = @(t,Z,X)(dKlcTdXT(t,Z,X));
        dUcTuTildZT = @(t,Z,X,fiBd)(...
            b1fn\(dUcThatdZT - dKbLcTdZT(t,Z,X)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X)*sat(s1(t,X)/fiBd) +...
             KbLcT(t,Z,X,mid,fiBd)*dsatS1FibddXT(t,X,fiBd))));
        dXtdZtT = sigDift;
        DUcTuTilDZT = @(t,Z,X,mid,fiBd)(dUcTuTildZT(t,Z,X,fiBd) +...
            dUcTuTildXT(t,Z,X,mid,fiBd)*dXtdZtT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZdupTilT);
        dZddmid = sigDift\dXtddmid;
        dUcTdHatdmid = @(t,mid)(-da1fndZT*dZddmid);
        dabsUcTdHatdmid = @(t,mid)(dUcTdHatdmid(t,mid)*sign(UcTdHat(t,mid)));
        dKlcTddmid = @(t,mid)((1 - D11fn)\(dFlc1dZT(fZtd(t,mid))*dZddmid + D11fn*dabsUcTdHatdmid(t,mid)));
        dKbLcTdmid = @(t,mid)(-dKlcTddmid(t,mid));
        dUcTuTildmid = @(t,X,mid,fiBd)(b1fn\(-dKbLcTdmid(t,mid)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd)(b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcT(t,Z,X,mid,fiBd)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd)(dUcTuTilduTilT(t,Z,X,mid,fiBd)*duTetadxT +...
            dUcTuTildmid(t,X,mid,fiBd)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxTx = @(t,xTil,xpTil)(...
            dUcTxdxTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Atu.'*AuAcx.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd)(dUcTuTildupTilT(t,Z,X,mid,fiBd)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxpTx = @(t,xTil,xpTil)(...
            dUcTxdxpTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
    case 13
        %Dinâmica longitudinal
        %Entrada:
        Kf = 10^5; %N/m
        fFbuUc = @(t,ub1)(y2d(t) - Kf*(ub1 - y2d(t)/kNM1));
        fub1 = @(ubTil)(Id8(1,:)*ubTil);
        fFbuUcuTil = @(t,ubTil)(fFbuUc(t,fub1(ubTil)));
        fub1U = @(ubTil)(Id3(1,:)*ubTil);
        fFbuUcuTilU = @(t,ubTil)(fFbuUc(t,fub1U(ubTil)));
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        
        %Obs.: xTil = AuAcx*uTilAc + Amix*mid + Afix*fiBd
        Id10 = eye(10);
        AuAcx = Id10(:,1:8);
        Amix = Id10(:,9);
        Afix = Id10(:,10);
        %Obs.: uTilAc = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL
        Atu = Id8(:,[2,6,8]);
        Auu = Id8(:,[1,3,7]);
        Atx = AuAcx*Atu;
        Avu = Id8(:,4);
        Awu = Id8(:,5);
        %Obs.: Uc = AtUc*Tx + AuUc*Fub
        Id2 = eye(2);
        AtUc = Id2(:,1);
        AuUc = Id2(:,2);
        
        fubTilAc = @(xTil)(AuAcx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTuAc = @(ubTilAc)(Atu.'*ubTilAc);
        fuUuAc = @(ubTilAc)(Auu.'*ubTilAc);
        fvuAc = @(ubTilAc)(Avu.'*ubTilAc);
        fwuAc = @(ubTilAc)(Awu.'*ubTilAc);
        fuTxTil = @(xTil)(fuTuAc(fubTilAc(xTil)));
        fuUxTil = @(xTil)(fuUuAc(fubTilAc(xTil)));
        fFbuUcx = @(t,xTil)(fFbuUcuTil(t,fubTilAc(xTil)));
        
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1ac + alfa1Nmk(dtnM1)*Cac + Kac);
        Knmkx = @(dtnM1)(Knmk(dtnM1)*AuAcx.');
        fFbuUcxnM1 = @(tnM1,xTilnM1)(fFbuUcx(tnM1,xTilnM1));
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fTilc(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfTilcdubT(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfTilcdubpT(fubTilAc(xTil),fubTilAc(xTilp)));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fTildc(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfTildcdubT(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfTildcdubpT(fubTilAc(xTil),fubTilAc(xTilp)));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dFbuUcdub1 = -Kf;
        dub1dubT = Id8(1,:);
        dFbuUcdubT = dFbuUcdub1*dub1dubT;
        dFbuUcdxT = dFbuUcdubT*dubTildxT;
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuAcx.';
        dfTildubppTx = @(xTil)(dfTildubppT(fubTilAc(xTil)));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        M1acx = M1ac*AuAcx.';
        Cacx = Cac*AuAcx.';
        
        %Camada limite e dinâmica interna:
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT =  Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
    case 14
        %Dinâmica longitudinal
        %Entrada:
        Kf = 10^5; %N/m
        fFbuUc = @(t,ub1)(y2d(t) - Kf*(ub1 - y2d(t)/kNM1));
        fub1 = @(ubTil)(Id8(1,:)*ubTil);
        fFbuUcuTil = @(t,ubTil)(fFbuUc(t,fub1(ubTil)));
        %Obs.: xTil = AuAcx*uTilAc + Amix*mid + Afix*fiBd
        Id10 = eye(10);
        AuAcx = Id10(:,1:8);
        Amix = Id10(:,9);
        Afix = Id10(:,10);
        %Obs.: uTilAc = Atu*uTeta + Auu*ubU
        Atu = Id8(:,[2,6,8]);
        Auu = Id8(:,[1,3,7]);
        Atx = AuAcx*Atu;
        Avu = Id8(:,4);
        Awu = Id8(:,5);
        %Obs.: Uc = AtUc*Tx + AuUc*Fub
        Id2 = eye(2);
        AtUc = Id2(:,1);
        AuUc = Id2(:,2);
        
        fubTilAc = @(xTil)(AuAcx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTuAc = @(ubTilAc)(Atu.'*ubTilAc);
        fuUuAc = @(ubTilAc)(Auu.'*ubTilAc);
        fuTxTil = @(xTil)(fuTuAc(fubTilAc(xTil)));
        fuUxTil = @(xTil)(fuUuAc(fubTilAc(xTil)));
        fFbuUcx = @(t,xTil)(fFbuUcuTil(t,fubTilAc(xTil)));
        
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1ac + alfa1Nmk(dtnM1)*Cac + Kac);
        Knmkx = @(dtnM1)(Knmk(dtnM1)*AuAcx.');
        fFbuUcxnM1 = @(tnM1,xTilnM1)(fFbuUcx(tnM1,xTilnM1));
        %Termos não lineares:
        %Sem atrito no contato:
        fTilcx = @(xTil,xTilp,xTilpp)(fTilSmAt(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfTilSmAtdubT(fubTilAc(xTil),fubTilAc(xTilp),fubTilAc(xTilpp)));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfTilSmAtdubpT(fubTilAc(xTil),fubTilAc(xTilp)));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dFbuUcdub1 = -Kf;
        dub1dubT = Id8(1,:);
        dFbuUcdubT = dFbuUcdub1*dub1dubT;
        dFbuUcdxT = dFbuUcdubT*dubTildxT;
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuAcx.';
        dfTildubppTx = @(xTil)(dfTildubppT(fubTilAc(xTil)));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        M1acx = M1ac*AuAcx.';
        Cacx = Cac*AuAcx.';
        
        %Camada limite e dinâmica interna:
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(fuTxTil(xTil),fuTxTil(xpTil)));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil)(fXt(fuTxTil(xTil),fuTxTil(xpTil)));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTeta,upTeta)(s1(t,fXt(uTeta,upTeta)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd)(...
            b1fn\(UcThat(t,Z,X) - KbLcT(t,Z,X,mid,fiBd)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd)(...
            UcT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        UcTx = @(t,xTil,xpTil)(...
            UcTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTetadxT = Atu.'*AuAcx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = -da1fndZT;
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z)(dda1fndZT*sign(da1fn(Z)));
        dFlc1dZT = @(Z)(dabsda1fndZT(Z));
        dabsUcThatdZT = @(t,Z,X)(dUcThatdZT*sign(UcThat(t,Z,X)));
        dKlcTdZT = @(t,Z,X)((1 - D11fn)\(dFlc1dZT(Z) + D11fn*dabsUcThatdZT(t,Z,X)));
        dabsUcThatdXT = @(t,Z,X)(dUcThatdXT*sign(UcThat(t,Z,X)));
        dKlcTdXT = @(t,Z,X)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X)));
        dKbLcTdZT = @(t,Z,X)(dKlcTdZT(t,Z,X));
        dKbLcTdXT = @(t,Z,X)(dKlcTdXT(t,Z,X));
        dUcTuTildZT = @(t,Z,X,fiBd)(...
            b1fn\(dUcThatdZT - dKbLcTdZT(t,Z,X)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X)*sat(s1(t,X)/fiBd) +...
             KbLcT(t,Z,X,mid,fiBd)*dsatS1FibddXT(t,X,fiBd))));
        dXtdZtT = sigDift;
        DUcTuTilDZT = @(t,Z,X,mid,fiBd)(dUcTuTildZT(t,Z,X,fiBd) +...
            dUcTuTildXT(t,Z,X,mid,fiBd)*dXtdZtT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZdupTilT);
        dZddmid = sigDift\dXtddmid;
        dUcTdHatdmid = @(t,mid)(-da1fndZT*dZddmid);
        dabsUcTdHatdmid = @(t,mid)(dUcTdHatdmid(t,mid)*sign(UcTdHat(t,mid)));
        dKlcTddmid = @(t,mid)((1 - D11fn)\(dFlc1dZT(fZtd(t,mid))*dZddmid + D11fn*dabsUcTdHatdmid(t,mid)));
        dKbLcTdmid = @(t,mid)(-dKlcTddmid(t,mid));
        dUcTuTildmid = @(t,X,mid,fiBd)(b1fn\(-dKbLcTdmid(t,mid)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd)(b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcT(t,Z,X,mid,fiBd)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd)(dUcTuTilduTilT(t,Z,X,mid,fiBd)*duTetadxT +...
            dUcTuTildmid(t,X,mid,fiBd)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxTx = @(t,xTil,xpTil)(...
            dUcTxdxTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Atu.'*AuAcx.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd)(dUcTuTildupTilT(t,Z,X,mid,fiBd)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxpTx = @(t,xTil,xpTil)(...
            dUcTxdxpTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
    case 15
        %Obs.: uTilAc = Atu*uTeta + Auu*ubU
        Atu = Id8(:,[2,6,8]);
        Auu = Id8(:,[1,3,7]);
        Avu = Id8(:,4);
        Awu = Id8(:,5);
        %Entrada:
        Tref = 10; %Nm
        fTxUc = @(t)(Tref);
        fFbuUc = @(t)(y2d(t));
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1ac + alfa1Nmk(dtnM1)*Cac + Kac);
        fTilnM1 = @(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)(...
            fTilSmAt(ubTilnM1,...
                     fupnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln),...
                     fuppnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)));
        %Sem atrito no contato:
        dfTildubT = @(ubTil,ubpTil,ubppTil)(...
            dfTilSmAtdubT(ubTil,ubpTil,ubppTil));
        dfTildubTnM1 = @(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)(...
            dfTildubT(ubTilnM1,...
                      fupnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln),...
                      fuppnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)));
        dfTildubpT = @(ubTil,ubpTil)(...
            dfTilSmAtdubpT(ubTil,ubpTil));
        dfTildubpTnM1 = @(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)(...
            dfTildubpT(ubTilnM1,...
                       fupnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)));
        %Derivadas com relação ao vetor aceleração:
        dfTildubppTnM1 = @(ubTilnM1)(dfTildubppT(ubTilnM1));        
        %Derivadas totais:
        dubpnM1dubnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dubppnM1dubnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilDubnM1T = @(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)(...
            dfTildubTnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln) +...
            dfTildubpTnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)*dubpnM1dubnM1T(dtnM1) +...
            dfTildubppTnM1(ubTilnM1)*dubppnM1dubnM1T(dtnM1));
    case 16
        %Obs.: uTilAc = Atu*uTeta + Auu*ubU
        Atu = Id8(:,[2,6,8]);
        Auu = Id8(:,[1,3,7]);
        Avu = Id8(:,4);
        Awu = Id8(:,5);
        %Entrada:
        Tref = 10; %Nm
        fTxUc = @(t)(Tref);
        Kf = 10^5; %N/m
        fFbuUc = @(t,ub1)(y2d(t) - Kf*(ub1 - y2d(t)/kNM1));
        fub1 = @(ubTil)(Id8(1,:)*ubTil);
        fFbuUcuTil = @(t,ubTil)(fFbuUc(t,fub1(ubTil)));
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1ac + alfa1Nmk(dtnM1)*Cac + Kac);
        fTilnM1 = @(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)(...
            fTilSmAt(ubTilnM1,...
                     fupnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln),...
                     fuppnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)));
        fFbuUcnM1 = @(tnM1,ubTilnM1)(fFbuUcuTil(tnM1,ubTilnM1));
        %Sem atrito no contato:
        dfTildubT = @(ubTil,ubpTil,ubppTil)(...
            dfTilSmAtdubT(ubTil,ubpTil,ubppTil));
        dfTildubTnM1 = @(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)(...
            dfTildubT(ubTilnM1,...
                      fupnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln),...
                      fuppnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)));
        dfTildubpT = @(ubTil,ubpTil)(...
            dfTilSmAtdubpT(ubTil,ubpTil));
        dfTildubpTnM1 = @(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)(...
            dfTildubpT(ubTilnM1,...
                       fupnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)));
        dFbuUcdub1 = -Kf;
        dub1dubT = Id8(1,:);
        dFbuUcdubT = dFbuUcdub1*dub1dubT;
        %Derivadas com relação ao vetor aceleração:
        dfTildubppTnM1 = @(ubTilnM1)(dfTildubppT(ubTilnM1));        
        %Derivadas totais:
        dubpnM1dubnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dubppnM1dubnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilDubnM1T = @(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)(...
            dfTildubTnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln) +...
            dfTildubpTnM1(dtnM1,ubTilnM1,ubTiln,ubpTiln,ubppTiln)*dubpnM1dubnM1T(dtnM1) +...
            dfTildubppTnM1(ubTilnM1)*dubppnM1dubnM1T(dtnM1));
    case 17
        %Dinâmica longitudinal
        %Entrada:
        Kf = 10^5; %N/m
        fFbuUc = @(t,ub1)(y2d(t) - Kf*(ub1 - y2d(t)/kNM1));
        fub1 = @(ubTil)(Id8(1,:)*ubTil);
        fFbuUcuTil = @(t,ubTil)(fFbuUc(t,fub1(ubTil)));
        fub1U = @(ubTil)(Id3(1,:)*ubTil);
        fFbuUcuTilU = @(t,ubTil)(fFbuUc(t,fub1U(ubTil)));
        %%
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        %%
        %Dinâmica de torção
        Id5 = eye(5);
        Atx = Id5(:,1:3);
        Amix = Id5(:,4);
        Afix = Id5(:,5);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        
        %Newmark:
        KtNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mt + alfa1Nmk(dtnM1)*Ct + Kt);
        KtNmkx = @(dtnM1)(KtNmk(dtnM1)*Atx.');
        Mtx = Mt*Atx.';
        Ctx = Ct*Atx.';
        
        %Camada limite e dinâmica interna:
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil)(fXt(Atx.'*xTil,Atx.'*xpTil));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil)(s1(t,fXt(uTil,upTil)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd)(...
            b1fn\(UcThat(t,Z,X) - KbLcT(t,Z,X,mid,fiBd)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd)(...
            UcT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        UcTx = @(t,xTil,xpTil)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = -da1fndZT;
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z)(dda1fndZT*sign(da1fn(Z)));
        dFlc1dZT = @(Z)(dabsda1fndZT(Z));
        dabsUcThatdZT = @(t,Z,X)(dUcThatdZT*sign(UcThat(t,Z,X)));
        dKlcTdZT = @(t,Z,X)((1 - D11fn)\(dFlc1dZT(Z) + D11fn*dabsUcThatdZT(t,Z,X)));
        dabsUcThatdXT = @(t,Z,X)(dUcThatdXT*sign(UcThat(t,Z,X)));
        dKlcTdXT = @(t,Z,X)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X)));
        dKbLcTdZT = @(t,Z,X)(dKlcTdZT(t,Z,X));
        dKbLcTdXT = @(t,Z,X)(dKlcTdXT(t,Z,X));
        dUcTuTildZT = @(t,Z,X,fiBd)(...
            b1fn\(dUcThatdZT - dKbLcTdZT(t,Z,X)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X)*sat(s1(t,X)/fiBd) +...
             KbLcT(t,Z,X,mid,fiBd)*dsatS1FibddXT(t,X,fiBd))));
        dXtdZtT = sigDift;
        DUcTuTilDZT = @(t,Z,X,mid,fiBd)(dUcTuTildZT(t,Z,X,fiBd) +...
            dUcTuTildXT(t,Z,X,mid,fiBd)*dXtdZtT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = sigDift\dXddmid;
        dUcTdHatdmid = @(t,mid)(-da1fndZT*dZddmid);
        dabsUcTdHatdmid = @(t,mid)(dUcTdHatdmid(t,mid)*sign(UcTdHat(t,mid)));
        dKlcTddmid = @(t,mid)((1 - D11fn)\(dFlc1dZT(fZtd(t,mid))*dZddmid + D11fn*dabsUcTdHatdmid(t,mid)));
        dKbLcTdmid = @(t,mid)(-dKlcTddmid(t,mid));
        dUcTuTildmid = @(t,X,mid,fiBd)(b1fn\(-dKbLcTdmid(t,mid)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd)(b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcT(t,Z,X,mid,fiBd)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd)(dUcTuTilduTilT(t,Z,X,mid,fiBd)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxTx = @(t,xTil,xpTil)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Atx.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd)(dUcTuTildupTilT(t,Z,X,mid,fiBd)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxpTx = @(t,xTil,xpTil)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
        %Derivadas Controle Proporcional longitudinal:
        dFbuUcdub1 = -Kf;
        dub1dubT = Id3(1,:);
        dFbuUcdubT = dFbuUcdub1*dub1dubT;
    case 18
        %Entrada:
        Tref = 2; %Nm
        fTx = @(t)(Tref);
        %Newmark:
        KtNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mt + alfa1Nmk(dtnM1)*Ct + Kt);
        FatTcnM1 = @(dtnM1,uTilTnM1,uTilTn,upTilTn,uppTilTn,y2nM1)(...
            FatTc(fupnM1(dtnM1,uTilTnM1,uTilTn,upTilTn,uppTilTn),y2nM1));
        dFatTcduTilTnM1 = @(dtnM1,uTilTnM1,uTilTn,upTilTn,uppTilTn,y2nM1)(...
            dFatTcduTilT(fupnM1(dtnM1,uTilTnM1,uTilTn,upTilTn,uppTilTn),y2nM1));
        FatTdcnM1 = @(dtnM1,uTilTnM1,uTilTn,upTilTn,uppTilTn,y2nM1)(...
            FatTdc(fupnM1(dtnM1,uTilTnM1,uTilTn,upTilTn,uppTilTn),y2nM1));
        dFatTdcduTilTnM1 = @(dtnM1,uTilTnM1,uTilTn,upTilTn,uppTilTn,y2nM1)(...
            dFatTdcduTilT(fupnM1(dtnM1,uTilTnM1,uTilTn,upTilTn,uppTilTn),y2nM1));
        %Entrada:
        Kf = 10^5; %N/m
        ub1Ref = 1*10^-4; %m
        fFbu = @(t,ub1)(-Kf*(ub1 - ub1Ref));
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
    case 19
        IdNgr = SigmaReF.'*IdNg*SigmaReF;
        nIdNgr = length(IdNgr);
        %Entradas:
        Tref = 2; %Nm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fTx = @(t)(Tref);
        dTxduTilrT = zeros(1,nIdNgr);
        Kf = 10^5; %N/m
        ub1Ref = 1*10^-4; %m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fFbu = @(ub1)(-Kf*(ub1 - ub1Ref));
        fub1 = @(uTilr)(IdNgr(1,:)*uTilr);
        fFbuTilr = @(uTilr)(fFbu(fub1(uTilr)));
        dFbudub1 = -Kf;
        dub1duTilrT = IdNgr(1,:);
        dFburduTilrT = dFbudub1*dub1duTilrT;
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1restr + alfa1Nmk(dtnM1)*Crestr + Krestr);
        fcEFbrnM1 = @(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)(...
            fcEFbr(ubrnM1,...
                   fupnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn),...
                   fuppnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)));
        dfcEFbrdubrTnM1 = @(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)(...
            dfcEFbrdubrT(ubrnM1,...
                         fupnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn),...
                         fuppnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)));
        dfcEFbrdubprTnM1 = @(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)(...
            dfcEFbrdubprT(ubrnM1,...
                          fupnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)));
        fdcEFbrnM1 = @(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)(...
            fdcEFbr(ubrnM1,...
                    fupnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn),...
                    fuppnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)));
        dfdcEFbrdubrTnM1 = @(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)(...
            dfdcEFbrdubrT(ubrnM1,...
                          fupnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn),...
                          fuppnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)));
        dfdcEFbrdubprTnM1 = @(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)(...
            dfdcEFbrdubprT(ubrnM1,...
                           fupnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)));
        dfEFbrdubpprTnM1 = @(ubrnM1)(dfEFbrdubpprT(ubrnM1));
        %Derivadas totais:
        dubprnM1dubrnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dubpprnM1dubrnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfcEFbrDubrnM1T = @(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)(...
            dfcEFbrdubrTnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn) +...
            dfcEFbrdubprTnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)*dubprnM1dubrnM1T(dtnM1) +...
            dfEFbrdubpprTnM1(ubrnM1)*dubpprnM1dubrnM1T(dtnM1));
        DfdcEFbrDubrnM1T = @(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)(...
            dfdcEFbrdubrTnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn) +...
            dfdcEFbrdubprTnM1(dtnM1,ubrnM1,ubrn,ubprn,ubpprn)*dubprnM1dubrnM1T(dtnM1) +...
            dfEFbrdubpprTnM1(ubrnM1)*dubpprnM1dubrnM1T(dtnM1));
        %Vetor de seleção de variáveis ciclicas:
        IdLVC = zeros(Ng,1);
        for i = 1:N+1
            n = 6*i;
            IdLVC = IdLVC + IdNg(:,n);
        end
        IdLVCr = SigmaReF.'*IdLVC;
        %Matrizes de seleção de variáveis:
        IdU = zeros(N+1,Ng);
        IdV = zeros(N+1,Ng);
        IdTetaZ = zeros(N+1,Ng);
        IdW = zeros(N+1,Ng);
        IdTetaY = zeros(N+1,Ng);
        IdTetaX = zeros(N+1,Ng);
        for i = 1:N+1
            nU = 6*i-5; nV = 6*i-4; nTz = 6*i-3; nW = 6*i-2; nTy = 6*i-1; nTx = 6*i;
            IdU(i,:) = IdNg(nU,:);
            IdV(i,:) = IdNg(nV,:);
            IdTetaZ(i,:) = IdNg(nTz,:);
            IdW(i,:) = IdNg(nW,:);
            IdTetaY(i,:) = IdNg(nTy,:);
            IdTetaX(i,:) = IdNg(nTx,:);
        end
    case 20
        IdNgr = SigmaReF.'*IdNg*SigmaReF;
        nIdNgr = length(IdNgr);
        %Dinâmica longitudinal
        %Entrada:
        Kf = 10^5; %N/m
        fFbuUc = @(t,ub1)(y2d(t) - Kf*(ub1 - y2d(t)/kNM1));
        Aub1uTil = IdNg(1,:);
        Aub1uTilr = Aub1uTil*SigmaReF;
        fub1 = @(ubTilr)(Aub1uTilr*ubTilr);
        fFbuUcuTil = @(t,ubTilr)(fFbuUc(t,fub1(ubTilr)));
        %Obs.: xTil = AuAcx*uTilAc + Amix*mid + Afix*fiBd
        IdNgrM2 = eye(nIdNgr+2);
        Aurx = IdNgrM2(:,1:nIdNgr);
        Amix = IdNgrM2(:,nIdNgr+1);
        Afix = IdNgrM2(:,nIdNgr+2);
        %Obs.: uTilAc = Atu*uTeta + Auu*ubU
        Atu = SigmaReF.'*IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = SigmaReF.'*IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Avu = SigmaReF.'*IdNg(:,6*(N1+1)-4);
        Awu = SigmaReF.'*IdNg(:,6*(N1+1)-2);
        %Obs.: Uc = AtUc*Tx + AuUc*Fub
        Id2 = eye(2);
        AtUc = Id2(:,1);
        AuUc = Id2(:,2);
        
        fubTilr = @(xTil)(Aurx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTur = @(ubTilr)(Atu.'*ubTilr);
        fuUur = @(ubTilr)(Auu.'*ubTilr);
        fuTxTil = @(xTil)(fuTur(fubTilr(xTil)));
        fuUxTil = @(xTil)(fuUur(fubTilr(xTil)));
        fFbuUcx = @(t,xTil)(fFbuUcuTil(t,fubTilr(xTil)));
        
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1restr + alfa1Nmk(dtnM1)*Crestr + Krestr);
        Knmkx = @(dtnM1)(Knmk(dtnM1)*Aurx.');
        fFbuUcxnM1 = @(tnM1,xTilnM1)(fFbuUcx(tnM1,xTilnM1));
        %Termos não lineares:
        %Sem atrito no contato:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFbr(fubTilr(xTil),fubTilr(xTilp),fubTilr(xTilpp)));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = Aurx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbrdubrT(fubTilr(xTil),fubTilr(xTilp),fubTilr(xTilpp)));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = Aurx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbrdubprT(fubTilr(xTil),fubTilr(xTilp)));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dFbuUcdub1 = -Kf;
        dub1dubT = Aub1uTilr;
        dFbuUcdubT = dFbuUcdub1*dub1dubT;
        dFbuUcdxT = dFbuUcdubT*dubTildxT;
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = Aurx.';
        dfTildubppTx = @(xTil)(dfEFbrdubpprT(fubTilr(xTil)));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        M1restrx = M1restr*Aurx.';
        Crestrx = Crestr*Aurx.';
        
        %Camada limite e dinâmica interna:
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(fuTxTil(xTil),fuTxTil(xpTil)));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil)(fXt(fuTxTil(xTil),fuTxTil(xpTil)));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTeta,upTeta)(s1(t,fXt(uTeta,upTeta)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd)(...
            b1fn\(UcThat(t,Z,X) - KbLcT(t,Z,X,mid,fiBd)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd)(...
            UcT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        UcTx = @(t,xTil,xpTil)(...
            UcTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTetadxT = Atu.'*Aurx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = -da1fndZT;
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z)(dda1fndZT*sign(da1fn(Z)));
        dFlc1dZT = @(Z)(dabsda1fndZT(Z));
        dabsUcThatdZT = @(t,Z,X)(dUcThatdZT*sign(UcThat(t,Z,X)));
        dKlcTdZT = @(t,Z,X)((1 - D11fn)\(dFlc1dZT(Z) + D11fn*dabsUcThatdZT(t,Z,X)));
        dabsUcThatdXT = @(t,Z,X)(dUcThatdXT*sign(UcThat(t,Z,X)));
        dKlcTdXT = @(t,Z,X)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X)));
        dKbLcTdZT = @(t,Z,X)(dKlcTdZT(t,Z,X));
        dKbLcTdXT = @(t,Z,X)(dKlcTdXT(t,Z,X));
        dUcTuTildZT = @(t,Z,X,fiBd)(...
            b1fn\(dUcThatdZT - dKbLcTdZT(t,Z,X)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X)*sat(s1(t,X)/fiBd) +...
             KbLcT(t,Z,X,mid,fiBd)*dsatS1FibddXT(t,X,fiBd))));
        dXtdZtT = sigDift;
        DUcTuTilDZT = @(t,Z,X,mid,fiBd)(dUcTuTildZT(t,Z,X,fiBd) +...
            dUcTuTildXT(t,Z,X,mid,fiBd)*dXtdZtT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZdupTilT);
        dZddmid = sigDift\dXtddmid;
        dUcTdHatdmid = @(t,mid)(-da1fndZT*dZddmid);
        dabsUcTdHatdmid = @(t,mid)(dUcTdHatdmid(t,mid)*sign(UcTdHat(t,mid)));
        dKlcTddmid = @(t,mid)((1 - D11fn)\(dFlc1dZT(fZtd(t,mid))*dZddmid + D11fn*dabsUcTdHatdmid(t,mid)));
        dKbLcTdmid = @(t,mid)(-dKlcTddmid(t,mid));
        dUcTuTildmid = @(t,X,mid,fiBd)(b1fn\(-dKbLcTdmid(t,mid)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd)(b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcT(t,Z,X,mid,fiBd)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd)(dUcTuTilduTilT(t,Z,X,mid,fiBd)*duTetadxT +...
            dUcTuTildmid(t,X,mid,fiBd)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxTx = @(t,xTil,xpTil)(...
            dUcTxdxTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Atu.'*Aurx.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd)(dUcTuTildupTilT(t,Z,X,mid,fiBd)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxpTx = @(t,xTil,xpTil)(...
            dUcTxdxpTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
        
        %Vetor de seleção de variáveis ciclicas:
        IdLVC = zeros(Ng,1);
        for i = 1:N+1
            n = 6*i;
            IdLVC = IdLVC + IdNg(:,n);
        end
        IdLVCr = SigmaReF.'*IdLVC;
        %Matrizes de seleção de variáveis:
        IdU = zeros(N+1,Ng);
        IdV = zeros(N+1,Ng);
        IdTetaZ = zeros(N+1,Ng);
        IdW = zeros(N+1,Ng);
        IdTetaY = zeros(N+1,Ng);
        IdTetaX = zeros(N+1,Ng);
        for i = 1:N+1
            nU = 6*i-5; nV = 6*i-4; nTz = 6*i-3; nW = 6*i-2; nTy = 6*i-1; nTx = 6*i;
            IdU(i,:) = IdNg(nU,:);
            IdV(i,:) = IdNg(nV,:);
            IdTetaZ(i,:) = IdNg(nTz,:);
            IdW(i,:) = IdNg(nW,:);
            IdTetaY(i,:) = IdNg(nTy,:);
            IdTetaX(i,:) = IdNg(nTx,:);
        end
    case 21
        %Entradas:
        Tref = 2; %Nm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fTx = @(t)(Tref);
        dTxdubT = zeros(1,Ng);
        Kf = 10^5; %N/m
        ub1Ref = 1*10^-4; %m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fFbu = @(ub1)(-Kf*(ub1 - ub1Ref));
        fub1 = @(uTil)(IdNg(1,:)*uTil);
        fFbuTil = @(uTil)(fFbu(fub1(uTil)));
        dFbudub1 = -Kf;
        dub1dubT = IdNg(1,:);
        dFbudubT = dFbudub1*dub1dubT;
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEF + alfa1Nmk(dtnM1)*CgEF + KgEF);
        fcEFbnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            fcEFb(ubnM1,...
                  fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                  fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfcEFbdubTnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfcEFbdubT(ubnM1,...
                       fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                       fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfcEFbdubpTnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfcEFbdubpT(ubnM1,...
                        fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        fdcEFbnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            fdcEFb(ubnM1,...
                   fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                   fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfdcEFbdubTnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfdcEFbdubT(ubnM1,...
                        fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                        fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfdcEFbdubpTnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfdcEFbdubpT(ubnM1,...
                         fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfEFbdubppTnM1 = @(ubrnM1)(dfEFbdubppT(ubrnM1));
        %Derivadas totais:
        dubpnM1dubnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dubppnM1dubnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfcEFbDubnM1T = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfcEFbdubTnM1(dtnM1,ubnM1,ubn,ubpn,ubppn) +...
            dfcEFbdubpTnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)*dubpnM1dubnM1T(dtnM1) +...
            dfEFbdubppTnM1(ubnM1)*dubppnM1dubnM1T(dtnM1));
        DfdcEFbDubnM1T = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfdcEFbdubTnM1(dtnM1,ubnM1,ubn,ubpn,ubppn) +...
            dfdcEFbdubpTnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)*dubpnM1dubnM1T(dtnM1) +...
            dfEFbdubppTnM1(ubnM1)*dubppnM1dubnM1T(dtnM1));
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMr = eye(Ng+rEstr);
        Auz = IdNgMr(:,1:Ng);
        Almbz = IdNgMr(:,Ng+1:Ng+rEstr);
        Knmkz = @(dtnM1)(...
            [Knmk(dtnM1), kNmk*GamaR.'; 
             kNmk*GamaR, zeros(rEstr)]);
        fcEFbznM1 = @(dtnM1,znM1,ubn,ubpn,ubppn)(...
            [fcEFbnM1(dtnM1,Auz.'*znM1,ubn,ubpn,ubppn); 
             zeros(rEstr,1)]);
        fdcEFbznM1 = @(dtnM1,znM1,ubn,ubpn,ubppn)(...
            [fdcEFbnM1(dtnM1,Auz.'*znM1,ubn,ubpn,ubppn); 
             zeros(rEstr,1)]);
        fiNmkzn = @(dtnM1,ubn,ubpn,ubppn)(...
            [M1gEF*fi2Nmk(dtnM1,ubn,ubpn,ubppn) + CgEF*fi1Nmk(dtnM1,ubn,ubpn,ubppn); 
             zeros(rEstr,1)]);
        Bz = [Bef; 
              zeros(rEstr,2)];
        fUcz = @(tnM1,znM1)([fTx(tnM1); 
                             fFbuTil(Auz.'*znM1)]); 
        dubdzT = Auz.';
        dfcEFzdzT = @(dtnM1,znM1,ubn,ubpn,ubppn)(...
            [DfcEFbDubnM1T(dtnM1,Auz.'*znM1,ubn,ubpn,ubppn)*dubdzT; 
             zeros(rEstr,Ng+rEstr)]);
        dfdcEFzdzT = @(dtnM1,znM1,ubn,ubpn,ubppn)(...
            [DfdcEFbDubnM1T(dtnM1,Auz.'*znM1,ubn,ubpn,ubppn)*dubdzT; 
             zeros(rEstr,Ng+rEstr)]);
        dTxzdzT = dTxdubT*dubdzT;
        dFbuzdzT = dFbudubT*dubdzT;
        fdUczdzT = [dTxzdzT; 
                    dFbuzdzT];
        %Vetor de seleção de variáveis ciclicas:
        IdLVC = zeros(Ng,1);
        for i = 1:N+1
            n = 6*i;
            IdLVC = IdLVC + IdNg(:,n);
        end
        %Matrizes de seleção de variáveis:
        IdU = zeros(N+1,Ng);
        IdV = zeros(N+1,Ng);
        IdTetaZ = zeros(N+1,Ng);
        IdW = zeros(N+1,Ng);
        IdTetaY = zeros(N+1,Ng);
        IdTetaX = zeros(N+1,Ng);
        for i = 1:N+1
            nU = 6*i-5; nV = 6*i-4; nTz = 6*i-3; nW = 6*i-2; nTy = 6*i-1; nTx = 6*i;
            IdU(i,:) = IdNg(nU,:);
            IdV(i,:) = IdNg(nV,:);
            IdTetaZ(i,:) = IdNg(nTz,:);
            IdW(i,:) = IdNg(nW,:);
            IdTetaY(i,:) = IdNg(nTy,:);
            IdTetaX(i,:) = IdNg(nTx,:);
        end
    case 22
        %Dinâmica longitudinal
        %Entrada:
        Kf = 10^6; %N/m
        fFbuUc = @(t,ub1)(y2d(t) - Kf*(ub1 - y2d(t)/kNM1));
        Aub1uTil = IdNg(1,:);
        fub1 = @(ubTil)(Aub1uTil*ubTil);
        fFbuUcuTil = @(t,ubTil)(fFbuUc(t,fub1(ubTil)));
        %Obs.: xTil = AuAcx*uTilAc + Amix*mid + Afix*fiBd
        IdNgM2 = eye(Ng+2);
        Aux = IdNgM2(:,1:Ng);
        Amix = IdNgM2(:,Ng+1);
        Afix = IdNgM2(:,Ng+2);
        %Obs.: uTilAc = Atu*uTeta + Auu*ubU
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        %Obs.: Uc = AtUc*Tx + AuUc*Fub
        Id2 = eye(2);
        AtUc = Id2(:,1);
        AuUc = Id2(:,2);
        
        fubTil = @(xTil)(Aux.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTu = @(ubTil)(Atu.'*ubTil);
        fuUu = @(ubTil)(Auu.'*ubTil);
        fuTxTil = @(xTil)(fuTu(fubTil(xTil)));
        fuUxTil = @(xTil)(fuUu(fubTil(xTil)));
        fFbuUcx = @(t,xTil)(fFbuUcuTil(t,fubTil(xTil)));
        
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEF + alfa1Nmk(dtnM1)*CgEF + KgEF);
        Knmkx = @(dtnM1)(Knmk(dtnM1)*Aux.');
        fFbuUcxnM1 = @(tnM1,xTilnM1)(fFbuUcx(tnM1,xTilnM1));
        %Termos não lineares:
        %Sem atrito no contato:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(fubTil(xTil),fubTil(xTilp),fubTil(xTilpp)));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = Aux.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(fubTil(xTil),fubTil(xTilp),fubTil(xTilpp)));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = Aux.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(fubTil(xTil),fubTil(xTilp)));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dFbuUcdub1 = -Kf;
        dub1dubT = Aub1uTil;
        dFbuUcdubT = dFbuUcdub1*dub1dubT;
        dFbuUcdxT = dFbuUcdubT*dubTildxT;
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = Aux.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(fubTil(xTil)));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        M1gEFx = M1gEF*Aux.';
        CgEFx = CgEF*Aux.';
        
        %Camada limite e dinâmica interna:
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(fuTxTil(xTil),fuTxTil(xpTil)));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil)(fXt(fuTxTil(xTil),fuTxTil(xpTil)));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTeta,upTeta)(s1(t,fXt(uTeta,upTeta)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd)(...
            b1fn\(UcThat(t,Z,X) - KbLcT(t,Z,X,mid,fiBd)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd)(...
            UcT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        UcTx = @(t,xTil,xpTil)(...
            UcTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTetadxT = Atu.'*Aux.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = -da1fndZT;
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z)(dda1fndZT*sign(da1fn(Z)));
        dFlc1dZT = @(Z)(dabsda1fndZT(Z));
        dabsUcThatdZT = @(t,Z,X)(dUcThatdZT*sign(UcThat(t,Z,X)));
        dKlcTdZT = @(t,Z,X)((1 - D11fn)\(dFlc1dZT(Z) + D11fn*dabsUcThatdZT(t,Z,X)));
        dabsUcThatdXT = @(t,Z,X)(dUcThatdXT*sign(UcThat(t,Z,X)));
        dKlcTdXT = @(t,Z,X)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X)));
        dKbLcTdZT = @(t,Z,X)(dKlcTdZT(t,Z,X));
        dKbLcTdXT = @(t,Z,X)(dKlcTdXT(t,Z,X));
        dUcTuTildZT = @(t,Z,X,fiBd)(...
            b1fn\(dUcThatdZT - dKbLcTdZT(t,Z,X)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X)*sat(s1(t,X)/fiBd) +...
             KbLcT(t,Z,X,mid,fiBd)*dsatS1FibddXT(t,X,fiBd))));
        dXtdZtT = sigDift;
        DUcTuTilDZT = @(t,Z,X,mid,fiBd)(dUcTuTildZT(t,Z,X,fiBd) +...
            dUcTuTildXT(t,Z,X,mid,fiBd)*dXtdZtT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZdupTilT);
        dZddmid = sigDift\dXtddmid;
        dUcTdHatdmid = @(t,mid)(-da1fndZT*dZddmid);
        dabsUcTdHatdmid = @(t,mid)(dUcTdHatdmid(t,mid)*sign(UcTdHat(t,mid)));
        dKlcTddmid = @(t,mid)((1 - D11fn)\(dFlc1dZT(fZtd(t,mid))*dZddmid + D11fn*dabsUcTdHatdmid(t,mid)));
        dKbLcTdmid = @(t,mid)(-dKlcTddmid(t,mid));
        dUcTuTildmid = @(t,X,mid,fiBd)(b1fn\(-dKbLcTdmid(t,mid)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd)(b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcT(t,Z,X,mid,fiBd)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd)(dUcTuTilduTilT(t,Z,X,mid,fiBd)*duTetadxT +...
            dUcTuTildmid(t,X,mid,fiBd)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxTx = @(t,xTil,xpTil)(...
            dUcTxdxTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Atu.'*Aux.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd)(dUcTuTildupTilT(t,Z,X,mid,fiBd)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxpTx = @(t,xTil,xpTil)(...
            dUcTxdxpTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgM2Mr = eye(Ng+2+rEstr);
        Axz = IdNgM2Mr(:,1:Ng+2);
        Almbz = IdNgM2Mr(:,Ng+2+1:Ng+2+rEstr);
        Knmkxz = @(dtnM1)(Knmkx(dtnM1)*Axz.');
        fTilcxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        fTildcxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        GamaRz = Almbz*GamaR;
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1)(KlcTdxnM1(tnM1,Axz.'*znM1));
        fFbuUcxznM1 = @(tnM1,znM1)(fFbuUcxnM1(tnM1,Axz.'*znM1));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)(...
            M1gEFx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + CgEFx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn));
        GamaRxz = GamaR*(Aux.'*Axz.');
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1)*dxTildzT);
        dFbuUczdzT = dFbuUcdxT*dxTildzT;
        
        %Vetor de seleção de variáveis ciclicas:
        IdLVC = zeros(Ng,1);
        for i = 1:N+1
            n = 6*i;
            IdLVC = IdLVC + IdNg(:,n);
        end
        IdLVCr = SigmaReF.'*IdLVC;
        %Matrizes de seleção de variáveis:
        IdU = zeros(N+1,Ng);
        IdV = zeros(N+1,Ng);
        IdTetaZ = zeros(N+1,Ng);
        IdW = zeros(N+1,Ng);
        IdTetaY = zeros(N+1,Ng);
        IdTetaX = zeros(N+1,Ng);
        for i = 1:N+1
            nU = 6*i-5; nV = 6*i-4; nTz = 6*i-3; nW = 6*i-2; nTy = 6*i-1; nTx = 6*i;
            IdU(i,:) = IdNg(nU,:);
            IdV(i,:) = IdNg(nV,:);
            IdTetaZ(i,:) = IdNg(nTz,:);
            IdW(i,:) = IdNg(nW,:);
            IdTetaY(i,:) = IdNg(nTy,:);
            IdTetaX(i,:) = IdNg(nTx,:);
        end
    case 23
        %Obs.: xTil = AuAcx*uTilAc + Amix*mid + Afix*fiBd
        IdNgM2 = eye(Ng+2);
        Aux = IdNgM2(:,1:Ng);
        Amix = IdNgM2(:,Ng+1);
        Afix = IdNgM2(:,Ng+2);
        %Obs.: uTilAc = Atu*uTeta + Auu*ubU
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        %Obs.: Uc = AtUc*Tx + AuUc*Fub
        Id2 = eye(2);
        AtUc = Id2(:,1);
        AuUc = Id2(:,2);
        
        fubTil = @(xTil)(Aux.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTu = @(ubTil)(Atu.'*ubTil);
        fuUu = @(ubTil)(Auu.'*ubTil);
        fuTxTil = @(xTil)(fuTu(fubTil(xTil)));
        fuUxTil = @(xTil)(fuUu(fubTil(xTil)));
        
        %Matriz de ganho, na equação da dinâmica:
        BefTx = Bef(:,1);
        
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEF + alfa1Nmk(dtnM1)*CgEF + KgEF);
        Knmkx = @(dtnM1)(Knmk(dtnM1)*Aux.');
        %Termos não lineares:
        %Sem atrito no contato:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(fubTil(xTil),fubTil(xTilp),fubTil(xTilpp)));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = Aux.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(fubTil(xTil),fubTil(xTilp),fubTil(xTilpp)));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = Aux.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(fubTil(xTil),fubTil(xTilp)));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = Aux.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(fubTil(xTil)));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        M1gEFx = M1gEF*Aux.';
        CgEFx = CgEF*Aux.';
        
        %Camada limite e dinâmica interna:
        GamaMix = GamaMi*Amix.';
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil)(KlcTd(t,fmidx(xTil)));
        KlcTdxnM1 = @(tnM1,xTilnM1)(KlcTdx(tnM1,xTilnM1));
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil)(dKlcTddmidT(t,fmidx(xTil)));
        dmiddxTilT = Amix.';
        dKlcTdxdxTilT = @(t,xTil)(dKlcTddmidTx(t,xTil)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1)(dKlcTdxdxTilT(tnM1,xTilnM1));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil)(fZtd(t,fmidx(xTil)));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(fuTxTil(xTil),fuTxTil(xpTil)));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil)(fXt(fuTxTil(xTil),fuTxTil(xpTil)));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTeta,upTeta)(s1(t,fXt(uTeta,upTeta)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd)(...
            b1fn\(UcThat(t,Z,X) - KbLcT(t,Z,X,mid,fiBd)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd)(...
            UcT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        UcTx = @(t,xTil,xpTil)(...
            UcTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        duTetadxT = Atu.'*Aux.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = -da1fndZT;
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z)(dda1fndZT*sign(da1fn(Z)));
        dFlc1dZT = @(Z)(dabsda1fndZT(Z));
        dabsUcThatdZT = @(t,Z,X)(dUcThatdZT*sign(UcThat(t,Z,X)));
        dKlcTdZT = @(t,Z,X)((1 - D11fn)\(dFlc1dZT(Z) + D11fn*dabsUcThatdZT(t,Z,X)));
        dabsUcThatdXT = @(t,Z,X)(dUcThatdXT*sign(UcThat(t,Z,X)));
        dKlcTdXT = @(t,Z,X)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X)));
        dKbLcTdZT = @(t,Z,X)(dKlcTdZT(t,Z,X));
        dKbLcTdXT = @(t,Z,X)(dKlcTdXT(t,Z,X));
        dUcTuTildZT = @(t,Z,X,fiBd)(...
            b1fn\(dUcThatdZT - dKbLcTdZT(t,Z,X)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X)*sat(s1(t,X)/fiBd) +...
             KbLcT(t,Z,X,mid,fiBd)*dsatS1FibddXT(t,X,fiBd))));
        dXtdZtT = sigDift;
        DUcTuTilDZT = @(t,Z,X,mid,fiBd)(dUcTuTildZT(t,Z,X,fiBd) +...
            dUcTuTildXT(t,Z,X,mid,fiBd)*dXtdZtT);
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd)(DUcTuTilDZT(t,Z,X,mid,fiBd)*dZdupTilT);
        dZddmid = sigDift\dXtddmid;
        dUcTdHatdmid = @(t,mid)(-da1fndZT*dZddmid);
        dabsUcTdHatdmid = @(t,mid)(dUcTdHatdmid(t,mid)*sign(UcTdHat(t,mid)));
        dKlcTddmid = @(t,mid)((1 - D11fn)\(dFlc1dZT(fZtd(t,mid))*dZddmid + D11fn*dabsUcTdHatdmid(t,mid)));
        dKbLcTdmid = @(t,mid)(-dKlcTddmid(t,mid));
        dUcTuTildmid = @(t,X,mid,fiBd)(b1fn\(-dKbLcTdmid(t,mid)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd)(b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcT(t,Z,X,mid,fiBd)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd)(dUcTuTilduTilT(t,Z,X,mid,fiBd)*duTetadxT +...
            dUcTuTildmid(t,X,mid,fiBd)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxTx = @(t,xTil,xpTil)(...
            dUcTxdxTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dupTildxpT = Atu.'*Aux.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd)(dUcTuTildupTilT(t,Z,X,mid,fiBd)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXt(uTil,upTil),mid,fiBd));
        dUcTxdxpTx = @(t,xTil,xpTil)(...
            dUcTxdxpTuTil(t,fuTxTil(xTil),fuTxTil(xpTil),Amix.'*xTil,Afix.'*xTil));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpdxT(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgM2Mr = eye(Ng+2+rEstrT);
        Axz = IdNgM2Mr(:,1:Ng+2);
        Almbz = IdNgM2Mr(:,Ng+2+1:Ng+2+rEstrT);
        Knmkxz = @(dtnM1)(Knmkx(dtnM1)*Axz.');
        fTilcxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        fTildcxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        GamaRz = Almbz*GamaRt;
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1)(KlcTdxnM1(tnM1,Axz.'*znM1));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)(...
            M1gEFx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + CgEFx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn));
        GamaRxz = GamaRt*(Aux.'*Axz.');
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1)*dxTildzT);
        
        %Vetor de seleção de variáveis ciclicas:
        IdLVC = zeros(Ng,1);
        for i = 1:N+1
            n = 6*i;
            IdLVC = IdLVC + IdNg(:,n);
        end
        IdLVCr = SigmaReF.'*IdLVC;
        %Matrizes de seleção de variáveis:
        IdU = zeros(N+1,Ng);
        IdV = zeros(N+1,Ng);
        IdTetaZ = zeros(N+1,Ng);
        IdW = zeros(N+1,Ng);
        IdTetaY = zeros(N+1,Ng);
        IdTetaX = zeros(N+1,Ng);
        for i = 1:N+1
            nU = 6*i-5; nV = 6*i-4; nTz = 6*i-3; nW = 6*i-2; nTy = 6*i-1; nTx = 6*i;
            IdU(i,:) = IdNg(nU,:);
            IdV(i,:) = IdNg(nV,:);
            IdTetaZ(i,:) = IdNg(nTz,:);
            IdW(i,:) = IdNg(nW,:);
            IdTetaY(i,:) = IdNg(nTy,:);
            IdTetaX(i,:) = IdNg(nTx,:);
        end
    case 24
        %%
        %Dinâmica longitudinal
        %Entrada:
        Kf = 10^5; %N/m
        ub10 = 0.0*10^-3;
        fFbuUc = @(t,ub1)(y2d(t) - Kf*(ub1 - y2d(t)/kNM1) + kNM1*ub10);
        fub1 = @(ubTil)(Id8(1,:)*ubTil);
        fFbuUcuTil = @(t,ubTil)(fFbuUc(t,fub1(ubTil)));
        fub1U = @(ubTil)(Id3(1,:)*ubTil);
        fFbuUcuTilU = @(t,ubTil)(fFbuUc(t,fub1U(ubTil)));
        %Newmark:
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        %%
        %Dinâmica de torção
        Id5 = eye(5);
        Atx = Id5(:,1:3);
        Amix = Id5(:,4);
        Afix = Id5(:,5);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        
        %Newmark:
        KtNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mt + alfa1Nmk(dtnM1)*Ct + Kt);
        KtNmkx = @(dtnM1)(KtNmk(dtnM1)*Atx.');
        Mtx = Mt*Atx.';
        Ctx = Ct*Atx.';
        FatTcx = @(xTilp,y2)(FatTc(ftx(xTilp),y2));
        FatTcxnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            FatTcx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1));
        FatTdcx = @(xTilp,y2)(FatTdc(ftx(xTilp),y2));
        FatTdcxnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            FatTdcx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1));
        
        %Camada limite e dinâmica interna:
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxTilpT = Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxTilpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT = Atx.';
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
        %Derivadas Controle Proporcional longitudinal:
        dFbuUcdub1 = -Kf;
        dub1dubT = Id3(1,:);
        dFbuUcdubT = dFbuUcdub1*dub1dubT;
    case 25
        %Entradas (restrições):
        %{
        hTeta = @(t)(y1d(t));
        hU = @(t)(y2d(t)/kNM1);
        %}
        %
        alfaTeta = 2.5*60; %RPM/s
        wMax = 240; %RPM
        alfaTeta = alfaTeta*(pi/30); %Aceleração máxima [rad/s2]
        wMax = wMax*(pi/30); %Velocidade máxima [rad/s]
        T0 = wMax/alfaTeta;
        hTeta = @(t)(wMax*(t/2)*double(t<=T0) + wMax*(t - T0/2)*double(t>T0));
        hU = @(t)(0);
        %}
        %{
        hTeta = @(t)(y1d(t));
        hU = @(t)(0);
        %}
        %Obs.:
        %(1) A expressão de hTeta foi definida na forma de deslocamento
        %(2) Mas é deduzida a partir da expressão da velocidade:
        %    w(t) = alfaTeta*t*double(t<=T0) + wMax*double(t>T0)
        %(3) O sinal recebido pelo servo-driver na bancada é de velocidade
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEF + alfa1Nmk(dtnM1)*CgEF + KgEF);
        fcEFbnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            fcEFb(ubnM1,...
                  fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                  fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfcEFbdubTnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfcEFbdubT(ubnM1,...
                       fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                       fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfcEFbdubpTnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfcEFbdubpT(ubnM1,...
                        fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        fdcEFbnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            fdcEFb(ubnM1,...
                   fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                   fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfdcEFbdubTnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfdcEFbdubT(ubnM1,...
                        fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                        fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfdcEFbdubpTnM1 = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfdcEFbdubpT(ubnM1,...
                         fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)));
        dfEFbdubppTnM1 = @(ubrnM1)(dfEFbdubppT(ubrnM1));
        %Derivadas totais:
        dubpnM1dubnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dubppnM1dubnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfcEFbDubnM1T = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfcEFbdubTnM1(dtnM1,ubnM1,ubn,ubpn,ubppn) +...
            dfcEFbdubpTnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)*dubpnM1dubnM1T(dtnM1) +...
            dfEFbdubppTnM1(ubnM1)*dubppnM1dubnM1T(dtnM1));
        DfdcEFbDubnM1T = @(dtnM1,ubnM1,ubn,ubpn,ubppn)(...
            dfdcEFbdubTnM1(dtnM1,ubnM1,ubn,ubpn,ubppn) +...
            dfdcEFbdubpTnM1(dtnM1,ubnM1,ubn,ubpn,ubppn)*dubpnM1dubnM1T(dtnM1) +...
            dfEFbdubppTnM1(ubnM1)*dubppnM1dubnM1T(dtnM1));
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+rEstr+2);
        Auz = IdNgMrM2(:,1:Ng);
        Almbz = IdNgMrM2(:,Ng+1:Ng+rEstr+2);
        GamaMA = [IdNg(6*1-5,:); 
                  IdNg(6*1,:); 
                  GamaR];
        hTil = @(t)([hU(t); 
                     hTeta(t); 
                     zeros(rEstr,1)]);
        rEstrH = length(hTil(0));
        Knmkz = @(dtnM1)(...
            [Knmk(dtnM1), kNmk*GamaMA.'; 
             kNmk*GamaMA, zeros(rEstr+2)]);
        fcEFbznM1 = @(dtnM1,tnM1,znM1,ubn,ubpn,ubppn)(...
            [fcEFbnM1(dtnM1,Auz.'*znM1,ubn,ubpn,ubppn); 
             -kNmk*hTil(tnM1)]);
        fdcEFbznM1 = @(dtnM1,tnM1,znM1,ubn,ubpn,ubppn)(...
            [fdcEFbnM1(dtnM1,Auz.'*znM1,ubn,ubpn,ubppn); 
             -kNmk*hTil(tnM1)]);
        fiNmkzn = @(dtnM1,ubn,ubpn,ubppn)(...
            [M1gEF*fi2Nmk(dtnM1,ubn,ubpn,ubppn) + CgEF*fi1Nmk(dtnM1,ubn,ubpn,ubppn); 
             zeros(rEstrH,1)]);
        dubdzT = Auz.';
        dfcEFzdzT = @(dtnM1,znM1,ubn,ubpn,ubppn)(...
            [DfcEFbDubnM1T(dtnM1,Auz.'*znM1,ubn,ubpn,ubppn)*dubdzT; 
             zeros(rEstrH,Ng+rEstrH)]);
        dfdcEFzdzT = @(dtnM1,znM1,ubn,ubpn,ubppn)(...
            [DfdcEFbDubnM1T(dtnM1,Auz.'*znM1,ubn,ubpn,ubppn)*dubdzT; 
             zeros(rEstrH,Ng+rEstrH)]);
        %Vetor de seleção de variáveis ciclicas:
        IdLVC = zeros(Ng,1);
        for i = 1:N+1
            n = 6*i;
            IdLVC = IdLVC + IdNg(:,n);
        end
        %Matrizes de seleção de variáveis:
        IdU = zeros(N+1,Ng);
        IdV = zeros(N+1,Ng);
        IdTetaZ = zeros(N+1,Ng);
        IdW = zeros(N+1,Ng);
        IdTetaY = zeros(N+1,Ng);
        IdTetaX = zeros(N+1,Ng);
        for i = 1:N+1
            nU = 6*i-5; nV = 6*i-4; nTz = 6*i-3; nW = 6*i-2; nTy = 6*i-1; nTx = 6*i;
            IdU(i,:) = IdNg(nU,:);
            IdV(i,:) = IdNg(nV,:);
            IdTetaZ(i,:) = IdNg(nTz,:);
            IdW(i,:) = IdNg(nW,:);
            IdTetaY(i,:) = IdNg(nTy,:);
            IdTetaX(i,:) = IdNg(nTx,:);
        end
    case 26
        %Entradas (restrições):
        %hU = @(t)(0);
        %hlU = @(t)(y2d(t)/(kNM1*(nUaj+1)));
        %hU = @(t)(y2d(t)/kNM1);
        
        %Dinâmica longitudinal (modelo em PC):
        %Newmark:
        %Obs.: zlnM1 = Aulz*ubl + Almblz*lmbl
        Aulz = Id4(:,1:3);
        Almblz = Id4(:,4);
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        GamaRu = Id3(1,:);
        fiuNmkn = @(dtnM1,ubn,ubpn,ubppn)(Mu*fi2Nmk(dtnM1,ubn,ubpn,ubppn) +...
            Cu*fi1Nmk(dtnM1,ubn,ubpn,ubppn));
        
        KuzNmk = @(dtnM1)([KuNmk(dtnM1), kNmk*GamaRu.'; kNmk*GamaRu, 0]);
        fiuNmkzn = @(dtnM1,ubn,ubpn,ubppn)([fiuNmkn(dtnM1,ubn,ubpn,ubppn); 0]);
        FnzNmk = @(tnM1,zlnM1)([fNuTil(Aulz.'*zlnM1); -kNmk*hlU(tnM1)]);
        dFnzdzT = @(zlnM1)([fdNudubTilT(Aulz.'*zlnM1), zeros(3,1); zeros(1,3), 0]);
        %Ayx = Id6(:,1:5);
        
        %Modelo em Elementos Finitos:
        %Obs.: xTil = AuEFx*ubTilEF + Amix*mid + Afix*fiBd
        IdNgM2 = eye(Ng+2);
        AuEFx = IdNgM2(:,1:Ng);
        Amix = IdNgM2(:,Ng+1);
        Afix = IdNgM2(:,Ng+2);
        
        M1gEFx = M1gEF*AuEFx.';
        CgEFx = CgEF*AuEFx.';
        KgEFx = KgEF*AuEFx.';
        
        %Newmark:
        %Termos lineares:
        KxNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEFx + alfa1Nmk(dtnM1)*CgEFx + KgEFx);
        fiTilxn = @(dtnM1,xn,xpn,xppn)(...
            M1gEFx*fi2Nmk(dtnM1,xn,xpn,xppn) + CgEFx*fi1Nmk(dtnM1,xn,xpn,xppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuEFx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuEFx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(AuEFx.'*xTil,AuEFx.'*xTilp));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fdcEFb(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuEFx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfdcEFbdubT(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuEFx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfdcEFbdubpT(AuEFx.'*xTil,AuEFx.'*xTilp));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuEFx.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(AuEFx.'*xTil));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Atx = AuEFx*Atu;
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        
        fubTilEF = @(xTil)(AuEFx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTuEF = @(ubTilEF)(Atu.'*ubTilEF);
        %fuUuEF = @(ubTilEF)(Auu.'*ubTilEF);
        %fvuEF = @(ubTilEF)(Avu.'*ubTilEF);
        %fwuEF = @(ubTilEF)(Awu.'*ubTilEF);
        %fuTxTil = @(xTil)(fuTuEF(fubTilEF(xTil)));
        %fuUxTil = @(xTil)(fuUuEF(fubTilEF(xTil)));
        
        %Matriz de ganho, na equação da dinâmica:
        BefTx = Bef(:,1);
        
        %Camada limite e dinâmica interna:
        %GamaMix = GamaMi*Amix.';
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT =  Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+2+rEstrT);
        Axz = IdNgMrM2(:,1:Ng+2);
        Almbz = IdNgMrM2(:,Ng+2+1:Ng+2+rEstrT);
        GamaMF = GamaRt;
        hTilMF = @(t)(fhRestr(t));
        GamaRz = Almbz*GamaMF;
        GamaRxz = GamaMF*(AuEFx.'*Axz.');
        Knmkxz = @(dtnM1)([KxNmk(dtnM1)*Axz.' + kNmk*GamaRz.'; 
            kNmk*GamaRxz]);
        fTilcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilMF(tnM1)]);
        fTildcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilMF(tnM1)]);
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        
        fmipdAtznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        %GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)([fiTilxn(dtnM1,xTiln,xTilpn,xTilppn); zeros(rEstrT,1)]);
        Befz = [BefTx; zeros(rEstrT,1)];
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        %{
        dtnM1 = 10^-4;
        znM1 = zeros(Ng+2+rEstrT,1);
        xTiln = zeros(Ng+2,1);
        xTilpn = zeros(Ng+2,1);
        xTilppn = zeros(Ng+2,1);
        ad1 = DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn);
        ad = dfTilczdzT(dtnM1,znM1,xTiln,xTilpn,xTilppn);
        %}
        dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        fdmipdAtdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdxnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
    case 27
        %Entradas (restrições):
        alfaTeta = 2.5*60; %RPM/s
        wMax = 240; %RPM
        alfaTeta = alfaTeta*(pi/30); %Aceleração máxima [rad/s2]
        wMax = wMax*(pi/30); %Velocidade máxima [rad/s]
        T0 = wMax/alfaTeta;
        fTetaX1 = @(t)(alfaTeta*(t^2/2)*double(t<=T0) + wMax*(t - T0/2)*double(t>T0));
        fTetapX1 = @(t)(alfaTeta*t*double(t<=T0) + wMax*double(t>T0));
        fTetappX1 = @(t)(alfaTeta*double(t<=T0) + 0*double(t>T0));
        %{
        hU = @(t)(0);
        %}
        %
        hU = @(t)(y2d(t)/kNM1);
        %}
        Auru = IdNg(:,[1:5,7:Ng]);
        AtXu = IdNg(:,6);
        M1uru = Auru.'*(M1gEF*Auru);
        M1tXu = Auru.'*(M1gEF*AtXu);
        Curu = Auru.'*(CgEF*Auru);
        CtXu = Auru.'*(CgEF*AtXu);
        Kuru = Auru.'*(KgEF*Auru);
        KtXu = Auru.'*(KgEF*AtXu);
        fub = @(uTx,tX1)(Auru*uTx + AtXu*tX1);
        
        fcTil = @(uTx,upTx,uppTx,tX1,tpX1,tppX1)(...
            Auru.'*...
            fcEFb(fub(uTx,tX1),fub(upTx,tpX1),fub(uppTx,tppX1)));
        dubduTxT = Auru;
        dfcTilduTxT = @(uTx,upTx,uppTx,tX1,tpX1,tppX1)(...
            Auru.'*...
            dfcEFbdubT(fub(uTx,tX1),fub(upTx,tpX1),fub(uppTx,tppX1))*...
            dubduTxT);
        dubpdupTxT = Auru;
        dfcTildupTxT = @(uTx,upTx,tX1,tpX1)(...
            Auru.'*...
            dfcEFbdubpT(fub(uTx,tX1),fub(upTx,tpX1))*...
            dubpdupTxT);
        
        fdcTil = @(uTx,upTx,uppTx,tX1,tpX1,tppX1)(...
            Auru.'*...
            fdcEFb(fub(uTx,tX1),fub(upTx,tpX1),fub(uppTx,tppX1)));
        dfdcTilduTxT = @(uTx,upTx,uppTx,tX1,tpX1,tppX1)(...
            Auru.'*...
            dfdcEFbdubT(fub(uTx,tX1),fub(upTx,tpX1),fub(uppTx,tppX1))*...
            dubduTxT);
        dfdcTildupTxT = @(uTx,upTx,tX1,tpX1)(...
            Auru.'*...
            dfdcEFbdubpT(fub(uTx,tX1),fub(upTx,tpX1))*...
            dubpdupTxT);
        
        dubppduppTxT = Auru;
        dfTilduppTxT = @(uTx,tX1)(...
            Auru.'*...
            dfEFbdubppT(fub(uTx,tX1))*...
            dubppduppTxT);
        
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1uru + alfa1Nmk(dtnM1)*Curu + Kuru);
        fcTilnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            fcTil(ubnM1,...
                  fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                  fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                  tX1nM1,tpX1nM1,tppX1nM1));
        dfcTilduTxTnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            dfcTilduTxT(ubnM1,...
                        fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                        fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                        tX1nM1,tpX1nM1,tppX1nM1));
        dfcTildupTxTnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,ubn,ubpn,ubppn)(...
            dfcTildupTxT(ubnM1,...
                         fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                         tX1nM1,tpX1nM1));
        fdcTilnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            fdcTil(ubnM1,...
                   fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                   fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                   tX1nM1,tpX1nM1,tppX1nM1));
        dfdcTilduTxTnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            dfdcTilduTxT(ubnM1,...
                         fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                         fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                         tX1nM1,tpX1nM1,tppX1nM1));
        dfdcTildupTxTnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,ubn,ubpn,ubppn)(...
            dfdcTildupTxT(ubnM1,...
                          fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                          tX1nM1,tpX1nM1));
        dfTilduppTxTnM1 = @(ubrnM1,tX1nM1)(dfTilduppTxT(ubrnM1,tX1nM1));
        %Derivadas totais:
        dubpnM1dubnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dubppnM1dubnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfcTilDuTxnM1T = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            dfcTilduTxTnM1(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn) +...
            dfcTildupTxTnM1(dtnM1,ubnM1,tX1nM1,tpX1nM1,ubn,ubpn,ubppn)*...
            dubpnM1dubnM1T(dtnM1) +...
            dfTilduppTxTnM1(ubnM1,tX1nM1)*...
            dubppnM1dubnM1T(dtnM1));
        DfdcTilDuTxnM1T = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            dfdcTilduTxTnM1(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn) +...
            dfdcTildupTxTnM1(dtnM1,ubnM1,tX1nM1,tpX1nM1,ubn,ubpn,ubppn)*...
            dubpnM1dubnM1T(dtnM1) +...
            dfTilduppTxTnM1(ubnM1,tX1nM1)*...
            dubppnM1dubnM1T(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgm1MrM1 = eye(Ng-1+rEstr+1);
        Aurz = IdNgm1MrM1(:,1:Ng-1);
        Almbz = IdNgm1MrM1(:,Ng-1+1:Ng-1+rEstr+1);
        GamaMA = [IdNg(6*1-5,:); 
                  GamaR];
        GamaMAu = GamaMA*Auru;
        GamaMAtX = GamaMA*AtXu;
        hTil = @(t)([hU(t);
                     zeros(rEstr,1)]);
        rEstrH = length(hTil(0));
        Knmkz = @(dtnM1)(...
            [Knmk(dtnM1), kNmk*GamaMAu.'; 
             kNmk*GamaMAu, zeros(rEstrH)]);
        fcTilznM1 = @(dtnM1,tnM1,znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            [fcTilnM1(dtnM1,Aurz.'*znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn); 
             -kNmk*hTil(tnM1)]);
        fdcTilznM1 = @(dtnM1,tnM1,znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            [fdcTilnM1(dtnM1,Aurz.'*znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn); 
             -kNmk*hTil(tnM1)]);
        KnmkTxz = @(tX1,tpX1,tppX1)(...
            [M1tXu*tppX1 + CtXu*tpX1 + KtXu*tX1; 
             zeros(rEstrH,1)]);
        fiNmkzn = @(dtnM1,ubn,ubpn,ubppn)(...
            [M1uru*fi2Nmk(dtnM1,ubn,ubpn,ubppn) + Curu*fi1Nmk(dtnM1,ubn,ubpn,ubppn); 
             zeros(rEstrH,1)]);
        duTxdzT = Aurz.';
        dfcTilzdzT = @(dtnM1,znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            [DfcTilDuTxnM1T(dtnM1,Aurz.'*znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)*duTxdzT; 
             zeros(rEstrH,Ng-1+rEstrH)]);
        dfdcTilzdzT = @(dtnM1,znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            [DfdcTilDuTxnM1T(dtnM1,Aurz.'*znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)*duTxdzT; 
             zeros(rEstrH,Ng-1+rEstrH)]);
        %Vetor de seleção de variáveis ciclicas:
        IdLVC = zeros(Ng,1);
        for i = 1:N+1
            n = 6*i;
            IdLVC = IdLVC + IdNg(:,n);
        end
        %Matrizes de seleção de variáveis:
        IdU = zeros(N+1,Ng);
        IdV = zeros(N+1,Ng);
        IdTetaZ = zeros(N+1,Ng);
        IdW = zeros(N+1,Ng);
        IdTetaY = zeros(N+1,Ng);
        IdTetaX = zeros(N+1,Ng);
        for i = 1:N+1
            nU = 6*i-5; nV = 6*i-4; nTz = 6*i-3; nW = 6*i-2; nTy = 6*i-1; nTx = 6*i;
            IdU(i,:) = IdNg(nU,:);
            IdV(i,:) = IdNg(nV,:);
            IdTetaZ(i,:) = IdNg(nTz,:);
            IdW(i,:) = IdNg(nW,:);
            IdTetaY(i,:) = IdNg(nTy,:);
            IdTetaX(i,:) = IdNg(nTx,:);
        end
    case 28
        %Entradas (restrições):
        tf = argumentos.tf; %s
        alfaTeta = 20; %RPM/s
        wMax = 10; %RPM
        u10 = 2.0; %mm Obs.: Deslocamento da posição original até o contato
        vU = 1.0; %mm/s Obs.: Velocidade de deslocamento
        aU = 5.0; %mm/s2 Obs.: Aceleração no deslocamento
        
        %Conversão de unidades:
        alfaTeta = alfaTeta*(pi/30); %Aceleração máxima [rad/s2]
        wMax = wMax*(pi/30); %Velocidade máxima [rad/s]
        u10 = u10*10^-3; %[m]
        vU = vU*10^-3; %[m/s]
        aU = aU*10^-3; %[m/s2]
        
        %Faixas de tempo:
        dT0 = wMax/alfaTeta; %[s] Tempo de acelaração angular
        dT1 = tf - dT0; %[s] Tempo com velocidade constante de entrada
        
        %Principais instantes de tempo:
        T0 = 0;
        T1 = T0 + dT0;
        T2 = T1 + dT1;

        fTetaX1 = @(t)(alfaTeta*t^2/2*double((T0<=t)&&(t<T1)) +...
            (alfaTeta*T1^2/2 + wMax*(t - T1))*double(T1<=t));
        
        fTetapX1 = @(t)(alfaTeta*t*double((T0<=t)&&(t<T1)) +...
            wMax*double(T1<=t));
        
        fTetappX1 = @(t)(alfaTeta*double((T0<=t)&&(t<T1)) +...
            0*double(T1<=t));
        
        hU = @(t)(0);
        
        Auru = IdNg(:,[1:5,7:Ng]);
        AtXu = IdNg(:,6);
        M1uru = Auru.'*(M1gEF*Auru);
        M1tXu = Auru.'*(M1gEF*AtXu);
        Curu = Auru.'*(CgEF*Auru);
        CtXu = Auru.'*(CgEF*AtXu);
        Kuru = Auru.'*(KgEF*Auru);
        KtXu = Auru.'*(KgEF*AtXu);
        fub = @(uTx,tX1)(Auru*uTx + AtXu*tX1);
        
        fcTil = @(uTx,upTx,uppTx,tX1,tpX1,tppX1)(...
            Auru.'*...
            fcEFb(fub(uTx,tX1),fub(upTx,tpX1),fub(uppTx,tppX1)));
        dubduTxT = Auru;
        dfcTilduTxT = @(uTx,upTx,uppTx,tX1,tpX1,tppX1)(...
            Auru.'*...
            dfcEFbdubT(fub(uTx,tX1),fub(upTx,tpX1),fub(uppTx,tppX1))*...
            dubduTxT);
        dubpdupTxT = Auru;
        dfcTildupTxT = @(uTx,upTx,tX1,tpX1)(...
            Auru.'*...
            dfcEFbdubpT(fub(uTx,tX1),fub(upTx,tpX1))*...
            dubpdupTxT);
        
        fdcTil = @(uTx,upTx,uppTx,tX1,tpX1,tppX1)(...
            Auru.'*...
            fdcEFb(fub(uTx,tX1),fub(upTx,tpX1),fub(uppTx,tppX1)));
        dfdcTilduTxT = @(uTx,upTx,uppTx,tX1,tpX1,tppX1)(...
            Auru.'*...
            dfdcEFbdubT(fub(uTx,tX1),fub(upTx,tpX1),fub(uppTx,tppX1))*...
            dubduTxT);
        dfdcTildupTxT = @(uTx,upTx,tX1,tpX1)(...
            Auru.'*...
            dfdcEFbdubpT(fub(uTx,tX1),fub(upTx,tpX1))*...
            dubpdupTxT);
        
        dubppduppTxT = Auru;
        dfTilduppTxT = @(uTx,tX1)(...
            Auru.'*...
            dfEFbdubppT(fub(uTx,tX1))*...
            dubppduppTxT);
        
        %Newmark:
        Knmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1uru + alfa1Nmk(dtnM1)*Curu + Kuru);
        fcTilnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            fcTil(ubnM1,...
                  fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                  fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                  tX1nM1,tpX1nM1,tppX1nM1));
        dfcTilduTxTnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            dfcTilduTxT(ubnM1,...
                        fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                        fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                        tX1nM1,tpX1nM1,tppX1nM1));
        dfcTildupTxTnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,ubn,ubpn,ubppn)(...
            dfcTildupTxT(ubnM1,...
                         fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                         tX1nM1,tpX1nM1));
        fdcTilnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            fdcTil(ubnM1,...
                   fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                   fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                   tX1nM1,tpX1nM1,tppX1nM1));
        dfdcTilduTxTnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            dfdcTilduTxT(ubnM1,...
                         fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                         fuppnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                         tX1nM1,tpX1nM1,tppX1nM1));
        dfdcTildupTxTnM1 = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,ubn,ubpn,ubppn)(...
            dfdcTildupTxT(ubnM1,...
                          fupnM1(dtnM1,ubnM1,ubn,ubpn,ubppn),...
                          tX1nM1,tpX1nM1));
        dfTilduppTxTnM1 = @(ubrnM1,tX1nM1)(dfTilduppTxT(ubrnM1,tX1nM1));
        %Derivadas totais:
        dubpnM1dubnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dubppnM1dubnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfcTilDuTxnM1T = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            dfcTilduTxTnM1(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn) +...
            dfcTildupTxTnM1(dtnM1,ubnM1,tX1nM1,tpX1nM1,ubn,ubpn,ubppn)*...
            dubpnM1dubnM1T(dtnM1) +...
            dfTilduppTxTnM1(ubnM1,tX1nM1)*...
            dubppnM1dubnM1T(dtnM1));
        DfdcTilDuTxnM1T = @(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            dfdcTilduTxTnM1(dtnM1,ubnM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn) +...
            dfdcTildupTxTnM1(dtnM1,ubnM1,tX1nM1,tpX1nM1,ubn,ubpn,ubppn)*...
            dubpnM1dubnM1T(dtnM1) +...
            dfTilduppTxTnM1(ubnM1,tX1nM1)*...
            dubppnM1dubnM1T(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgm1MrM1 = eye(Ng-1+rEstr+1);
        Aurz = IdNgm1MrM1(:,1:Ng-1);
        Almbz = IdNgm1MrM1(:,Ng-1+1:Ng-1+rEstr+1);
        GamaMA = [IdNg(6*1-5,:); 
                  GamaR];
        GamaMAu = GamaMA*Auru;
        GamaMAtX = GamaMA*AtXu;
        hTil = @(t)([hU(t);
                     zeros(rEstr,1)]);
        rEstrH = length(hTil(0));
        Knmkz = @(dtnM1)(...
            [Knmk(dtnM1), kNmk*GamaMAu.'; 
             kNmk*GamaMAu, zeros(rEstrH)]);
        fcTilznM1 = @(dtnM1,tnM1,znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            [fcTilnM1(dtnM1,Aurz.'*znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn); 
             -kNmk*hTil(tnM1)]);
        fdcTilznM1 = @(dtnM1,tnM1,znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            [fdcTilnM1(dtnM1,Aurz.'*znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn); 
             -kNmk*hTil(tnM1)]);
        KnmkTxz = @(tX1,tpX1,tppX1)(...
            [M1tXu*tppX1 + CtXu*tpX1 + KtXu*tX1; 
             zeros(rEstrH,1)]);
        fiNmkzn = @(dtnM1,ubn,ubpn,ubppn)(...
            [M1uru*fi2Nmk(dtnM1,ubn,ubpn,ubppn) + Curu*fi1Nmk(dtnM1,ubn,ubpn,ubppn); 
             zeros(rEstrH,1)]);
        duTxdzT = Aurz.';
        dfcTilzdzT = @(dtnM1,znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            [DfcTilDuTxnM1T(dtnM1,Aurz.'*znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)*duTxdzT; 
             zeros(rEstrH,Ng-1+rEstrH)]);
        dfdcTilzdzT = @(dtnM1,znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)(...
            [DfdcTilDuTxnM1T(dtnM1,Aurz.'*znM1,tX1nM1,tpX1nM1,tppX1nM1,ubn,ubpn,ubppn)*duTxdzT; 
             zeros(rEstrH,Ng-1+rEstrH)]);
        %Vetor de seleção de variáveis ciclicas:
        IdLVC = zeros(Ng,1);
        for i = 1:N+1
            n = 6*i;
            IdLVC = IdLVC + IdNg(:,n);
        end
        %Matrizes de seleção de variáveis:
        IdU = zeros(N+1,Ng);
        IdV = zeros(N+1,Ng);
        IdTetaZ = zeros(N+1,Ng);
        IdW = zeros(N+1,Ng);
        IdTetaY = zeros(N+1,Ng);
        IdTetaX = zeros(N+1,Ng);
        for i = 1:N+1
            nU = 6*i-5; nV = 6*i-4; nTz = 6*i-3; nW = 6*i-2; nTy = 6*i-1; nTx = 6*i;
            IdU(i,:) = IdNg(nU,:);
            IdV(i,:) = IdNg(nV,:);
            IdTetaZ(i,:) = IdNg(nTz,:);
            IdW(i,:) = IdNg(nW,:);
            IdTetaY(i,:) = IdNg(nTy,:);
            IdTetaX(i,:) = IdNg(nTx,:);
        end
        ad = 1;
    case 29
        %%
        %Dinâmica longitudinal (modelo em PC sem acoplamento):
        %Newmark:
        %Obs.: zlnM1 = Aulz*ubl + Almblz*lmbl
        Aulz = Id4(:,1:3);
        Almblz = Id4(:,4);
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        GamaRu = Id3(1,:);
        fiuNmkn = @(dtnM1,ubn,ubpn,ubppn)(Mu*fi2Nmk(dtnM1,ubn,ubpn,ubppn) +...
            Cu*fi1Nmk(dtnM1,ubn,ubpn,ubppn));
        
        KuzNmk = @(dtnM1)([KuNmk(dtnM1), kNmk*GamaRu.'; kNmk*GamaRu, 0]);
        fiuNmkzn = @(dtnM1,ubn,ubpn,ubppn)([fiuNmkn(dtnM1,ubn,ubpn,ubppn); 0]);
        FnzNmk = @(tnM1,zlnM1)([fNuTil(Aulz.'*zlnM1); -kNmk*hlU(tnM1)]);
        dFnzdzT = @(zlnM1)([fdNudubTilT(Aulz.'*zlnM1), zeros(3,1); zeros(1,3), 0]);
        %Ayx = Id6(:,1:5);
        
        %%
        %Modelo em PC com acoplamento:
        %Newmark:
        %Termos lineares:
        KacNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1ac + alfa1Nmk(dtnM1)*Cac + Kac);
        fiAcTilun = @(dtnM1,un,upn,uppn)(...
            M1ac*fi2Nmk(dtnM1,un,upn,uppn) + Cac*fi1Nmk(dtnM1,un,upn,uppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcunM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            fTilc(uTilnM1,...
                  fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                  fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTilcdubTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcdubT(uTilnM1,...
                       fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                       fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTilcdubpTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcdubpT(uTilnM1,...
                        fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        %Modelos descontínuos de atrito:
        fTildcunM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            fTildc(uTilnM1,...
                   fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                   fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTildcdubTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfdcEFbdubT(uTilnM1,...
                        fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                        fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTildcdubpTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTildcdubpT(uTilnM1,...
                         fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dfTildubppTnM1 = @(uTilnM1)(...
            dfTildubppT(uTilnM1));        
        %Derivadas totais:
        dupnM1dunM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        duppnM1dunM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDunM1T = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcdubTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn) +...
            dfTilcdubpTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)*dupnM1dunM1T(dtnM1) +...
            dfTildubppTnM1(uTilnM1)*duppnM1dunM1T(dtnM1));
        DfTildcDunM1T = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTildcdubTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn) +...
            dfTildcdubpTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)*dupnM1dunM1T(dtnM1) +...
            dfTildubppTnM1(uTilnM1)*duppnM1dunM1T(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        nAc = 8;
        IdnAcM2 = eye(nAc+2);
        AuAcz = IdnAcM2(:,1:nAc);
        AlmbAcz = IdnAcM2(:,nAc+1:nAc+2);
        GamaAc = [Id8(1,:); Id8(2,:)];
        hU = @(t)(hlU(t));
        hTilAc = @(t)([hU(t); 0]);
        KacNmkz = @(dtnM1)([KacNmk(dtnM1)*AuAcz.' + kNmk*GamaAc.'*AlmbAcz.'; 
            kNmk*GamaAc*AuAcz.']);
        fTilcuAcznM1 = @(dtnM1,tnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [fTilcunM1(dtnM1,AuAcz.'*znM1,uTiln,uTilpn,uTilppn); 
             -kNmk*hTilAc(tnM1)]);
        fTildcuAcznM1 = @(dtnM1,tnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [fTildcxnM1(dtnM1,AuAcz.'*znM1,uTiln,uTilpn,uTilppn); 
             -kNmk*hTilAc(tnM1)]);
        fiNmkuAczn = @(dtnM1,uTiln,uTilpn,uTilppn)(...
            [fiAcTilun(dtnM1,uTiln,uTilpn,uTilppn); 
             zeros(2,1)]);
        
        duAcdzT = AuAcz.';
        dfTilcAczdzT = @(dtnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [DfTilcDunM1T(dtnM1,AuAcz.'*znM1,uTiln,uTilpn,uTilppn)*duAcdzT; 
             zeros(2,nAc+2)]);
        dfTildcAczdzT = @(dtnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [DfTildcDunM1T(dtnM1,AuAcz.'*znM1,uTiln,uTilpn,uTilppn)*duAcdzT; 
             zeros(2,nAc+2)]);
        
        %%
        %Modelo em Elementos Finitos:
        %Obs.: xTil = AuEFx*ubTilEF
        AuEFx = IdNg;
        
        M1gEFx = M1gEF*AuEFx.';
        CgEFx = CgEF*AuEFx.';
        KgEFx = KgEF*AuEFx.';
        
        %Newmark:
        %Termos lineares:
        KxNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEFx + alfa1Nmk(dtnM1)*CgEFx + KgEFx);
        fiTilxn = @(dtnM1,xn,xpn,xppn)(...
            M1gEFx*fi2Nmk(dtnM1,xn,xpn,xppn) + CgEFx*fi1Nmk(dtnM1,xn,xpn,xppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuEFx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuEFx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(AuEFx.'*xTil,AuEFx.'*xTilp));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fdcEFb(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuEFx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfdcEFbdubT(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuEFx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfdcEFbdubpT(AuEFx.'*xTil,AuEFx.'*xTilp));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuEFx.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(AuEFx.'*xTil));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Atx = AuEFx*Atu;
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        
        fubTilEF = @(xTil)(AuEFx.'*xTil);
        fuTuEF = @(ubTilEF)(Atu.'*ubTilEF);
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+rEstr+2);
        Axz = IdNgMrM2(:,1:Ng);
        Almbz = IdNgMrM2(:,Ng+1:Ng+rEstr+2);
        GamaMF = GamaRtX;
        hU = @(t)(hlU(t));
        hTilMF = @(t)([hU(t); 0; hRestr]);
        GamaRz = Almbz*GamaMF;
        GamaRxz = GamaMF*(AuEFx.'*Axz.');
        Knmkxz = @(dtnM1)([KxNmk(dtnM1)*Axz.' + kNmk*GamaRz.'; 
            kNmk*GamaRxz]);
        fTilcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilMF(tnM1)]);
        fTildcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilMF(tnM1)]);
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)(...
            [fiTilxn(dtnM1,xTiln,xTilpn,xTilppn); zeros(rEstr+2,1)]);
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstr+2,Ng+rEstr+2)]);
        dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstr+2,Ng+rEstr+2)]);
    case 30
        %%
        %Entrada:
        htX = @(t)(Y1d(t) - Y1d(0));
        fTx = @(t)(htX(t));
        
        %Modelo em PC com acoplamento:
        %Newmark:
        %Termos lineares:
        KacNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1ac + alfa1Nmk(dtnM1)*Cac + Kac);
        fiAcTilun = @(dtnM1,un,upn,uppn)(...
            M1ac*fi2Nmk(dtnM1,un,upn,uppn) + Cac*fi1Nmk(dtnM1,un,upn,uppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcunM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            fTilc(uTilnM1,...
                  fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                  fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTilcdubTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcdubT(uTilnM1,...
                       fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                       fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTilcdubpTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcdubpT(uTilnM1,...
                        fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        %Modelos descontínuos de atrito:
        fTildcunM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            fTildc(uTilnM1,...
                   fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                   fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTildcdubTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfdcEFbdubT(uTilnM1,...
                        fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                        fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTildcdubpTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTildcdubpT(uTilnM1,...
                         fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dfTildubppTnM1 = @(uTilnM1)(...
            dfTildubppT(uTilnM1));        
        %Derivadas totais:
        dupnM1dunM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        duppnM1dunM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDunM1T = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcdubTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn) +...
            dfTilcdubpTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)*dupnM1dunM1T(dtnM1) +...
            dfTildubppTnM1(uTilnM1)*duppnM1dunM1T(dtnM1));
        DfTildcDunM1T = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTildcdubTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn) +...
            dfTildcdubpTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)*dupnM1dunM1T(dtnM1) +...
            dfTildubppTnM1(uTilnM1)*duppnM1dunM1T(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        nAc = 8;
        IdnAcM2 = eye(nAc+2);
        AuAcz = IdnAcM2(:,1:nAc);
        AlmbAcz = IdnAcM2(:,nAc+1:nAc+2);
        GamaAc = [Id8(1,:); Id8(2,:)];
        hU = @(t)(hlU(t));
        hTilAc = @(t)([hU(t); htX(t)]);
        KacNmkz = @(dtnM1)([KacNmk(dtnM1)*AuAcz.' + kNmk*GamaAc.'*AlmbAcz.'; 
            kNmk*GamaAc*AuAcz.']);
        fTilcuAcznM1 = @(dtnM1,tnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [fTilcunM1(dtnM1,AuAcz.'*znM1,uTiln,uTilpn,uTilppn); 
             -kNmk*hTilAc(tnM1)]);
        fTildcuAcznM1 = @(dtnM1,tnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [fTildcxnM1(dtnM1,AuAcz.'*znM1,uTiln,uTilpn,uTilppn); 
             -kNmk*hTilAc(tnM1)]);
        fiNmkuAczn = @(dtnM1,uTiln,uTilpn,uTilppn)(...
            [fiAcTilun(dtnM1,uTiln,uTilpn,uTilppn); 
             zeros(2,1)]);
        
        duAcdzT = AuAcz.';
        dfTilcAczdzT = @(dtnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [DfTilcDunM1T(dtnM1,AuAcz.'*znM1,uTiln,uTilpn,uTilppn)*duAcdzT; 
             zeros(2,nAc+2)]);
        dfTildcAczdzT = @(dtnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [DfTildcDunM1T(dtnM1,AuAcz.'*znM1,uTiln,uTilpn,uTilppn)*duAcdzT; 
             zeros(2,nAc+2)]);
        
        %%
        %Modelo em Elementos Finitos:
        %Newmark:
        %Termos lineares:
        KefNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEF + alfa1Nmk(dtnM1)*CgEF + KgEF);
        fiTiln = @(dtnM1,un,upn,uppn)(...
            M1gEF*fi2Nmk(dtnM1,un,upn,uppn) + CgEF*fi1Nmk(dtnM1,un,upn,uppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcEF = @(uTil,uTilp,uTilpp)(fcEFb(uTil,uTilp,uTilpp));
        fTilcEFnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            fTilcEF(uTilnM1,...
                   fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                   fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTilcEFdubT = @(uTil,uTilp,uTilpp)(...
            dfcEFbdubT(uTil,uTilp,uTilpp));
        dfTilcEFdubTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcEFdubT(uTilnM1,...
                         fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                         fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTilcEFdubpT = @(uTil,uTilp)(...
            dfcEFbdubpT(uTil,uTilp));
        dfTilcEFdubpTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcEFdubpT(uTilnM1,...
                          fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        %Modelos descontínuos de atrito:
        fTildcEF = @(uTil,uTilp,uTilpp)(fdcEFb(uTil,uTilp,uTilpp));
        fTildcEFnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            fTildcEF(uTilnM1,...
                     fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                   fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTildcEFdubT = @(uTil,uTilp,uTilpp)(...
            dfdcEFbdubT(uTil,uTilp,uTilpp));
        dfTildcEFdubTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTildcEFdubT(uTilnM1,...
                          fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn),...
                          fuppnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        dfTildcEFdubpT = @(uTil,uTilp)(...
            dfdcEFbdubpT(uTil,uTilp));
        dfTildcEFdubpTnM1 = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTildcEFdubpT(uTilnM1,...
                           fupnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dfTilEFdubppT = @(uTil)(dfEFbdubppT(uTil));
        dfTilEFdubppTnM1 = @(uTilnM1)(...
            dfTilEFdubppT(uTilnM1));        
        %Derivadas totais:
        dupnM1dunM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        duppnM1dunM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcEfDubnM1T = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTilcEFdubTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn) +...
            dfTilcEFdubpTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)*dupnM1dunM1T(dtnM1) +...
            dfTilEFdubppTnM1(uTilnM1)*duppnM1dunM1T(dtnM1));
        DfTildcEfDubnM1T = @(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)(...
            dfTildcEFdubTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn) +...
            dfTildcEFdubpTnM1(dtnM1,uTilnM1,uTiln,uTilpn,uTilppn)*dupnM1dunM1T(dtnM1) +...
            dfTilEFdubppTnM1(uTilnM1)*duppnM1dunM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+rEstr+2);
        AuEFz = IdNgMrM2(:,1:Ng);
        Almbz = IdNgMrM2(:,Ng+1:Ng+rEstr+2);
        GamaMF = GamaRtX;
        hTilMF = @(t)([hU(t); htX(t); hRestr]);
        GamaRlmbz = Almbz*GamaMF;
        GamaRuz = GamaMF*(AuEFz.');
        KefNmkz = @(dtnM1)([KefNmk(dtnM1)*AuEFz.' + kNmk*GamaRlmbz.'; 
            kNmk*GamaRuz]);
        fTilcEFznM1 = @(dtnM1,tnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [fTilcEFnM1(dtnM1,AuEFz.'*znM1,uTiln,uTilpn,uTilppn); 
             -kNmk*hTilMF(tnM1)]);
        fTildcEFznM1 = @(dtnM1,tnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [fTildcnM1(dtnM1,AuEFz.'*znM1,uTiln,uTilpn,uTilppn); 
             -kNmk*hTilMF(tnM1)]);
        fiNmkuzn = @(dtnM1,uTiln,uTilpn,uTilppn)(...
            [fiTiln(dtnM1,uTiln,uTilpn,uTilppn); zeros(rEstr+2,1)]);
        
        duTildzT = AuEFz.';
        dfTilcEFzdzT = @(dtnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [DfTilcEfDubnM1T(dtnM1,AuEFz.'*znM1,uTiln,uTilpn,uTilppn)*duTildzT; 
             zeros(rEstr+2,Ng+rEstr+2)]);
        dfTildcEFzdzT = @(dtnM1,znM1,uTiln,uTilpn,uTilppn)(...
            [DfTildcEfDubnM1T(dtnM1,AuEFz.'*znM1,uTiln,uTilpn,uTilppn)*duTildzT; 
             zeros(rEstr+2,Ng+rEstr+2)]);
         
        %%
        %Modelo aproximado reduzido:
        IdNgmrEstr = eye(Ng-rEstr);
        Au1ur = IdNgmrEstr(:,1);
        AuSur = IdNgmrEstr(:,2:Ng-rEstr);
        M1s = (Pk.'*AuSur)*AuSur.'*M1restr*AuSur*(AuSur.'*Pk);
        M1s = (1/2)*(M1s + M1s.');
        Cs = (Pk.'*AuSur)*AuSur.'*Crestr*AuSur*(AuSur.'*Pk);
        Cs = (1/2)*(Cs + Cs.');
        Ks = (Pk.'*AuSur)*AuSur.'*Krestr*AuSur*(AuSur.'*Pk);
        Ks = (1/2)*(Ks + Ks.');
        fubEFr = @(ub1,ubs)(Au1ur*ub1 + AuSur*ubs);
        fKls = @(ub1,ubs)((Pk.'*AuSur)*AuSur.'*fKlr(fubEFr(ub1,ubs))*AuSur*(AuSur.'*Pk));
        fKls = @(ub1,ubs)((1/2)*(fKls(ub1,ubs) + fKls(ub1,ubs).'));
        fKs = @(ub1,ubrk)(Ks + fKls(ub1,(AuSur.'*Pk)*ubrk));
        Bs = (Pk.'*AuSur)*AuSur.'*Br(:,1);
        
        M1s1 = (Pk.'*AuSur)*AuSur.'*M1restr*Au1ur;
        Cs1 = (Pk.'*AuSur)*AuSur.'*Crestr*Au1ur;
        Ks1 = (Pk.'*AuSur)*AuSur.'*Krestr*Au1ur;
        fKls1 = @(ub1,ubs)((Pk.'*AuSur)*AuSur.'*fKlr(fubEFr(ub1,ubs))*Au1ur);
        fKs1 = @(ub1,ubrk)(Ks1 + fKls1(ub1,(AuSur.'*Pk)*ubrk));
        
        yp2dAj = @(t)(dNr0bAjdw(y1d(t))*y1pd(t));
        hUp = @(t)(yp2dAj(t)/kNM1);
        ypp2dAj = @(t)(d2Nr0bAjdw2(y1d(t))*y1pd(t) + dNr0bAjdw(y1d(t))*y1ppd(t));
        hUpp = @(t)(ypp2dAj(t)/kNM1);
        
        %Entrada:
        fub1IN = @(t)(hU(t));
        fubp1IN = @(t)(hUp(t));
        fubpp1IN = @(t)(hUpp(t));
        %Formução de Newmark:
        Khats = @(dtnM1,tnM1,ubrknM1)(alfa2Nmk(dtnM1)*M1s + alfa1Nmk(dtnM1)*Cs + fKs(fub1IN(tnM1),ubrknM1));
        fiTilskn = @(dtnM1,ubrkn,ubprkn,ubpprkn)(...
            M1s*fi2Nmk(dtnM1,ubrkn,ubprkn,ubpprkn) + Cs*fi1Nmk(dtnM1,ubrkn,ubprkn,ubpprkn));
        ad = 1;
    case 31
        %Entrada:
        Tref = 10; %Nm
        fTx = @(t)(Tref);
        %Formução de Newmark:
        KhatTx = @(dtnM1,ubrknM1)(alfa2Nmk(dtnM1)*MtX + alfa1Nmk(dtnM1)*CtX + fKtX(ubrknM1));
        fiTilrkn = @(dtnM1,ubrkn,ubprkn,ubpprkn)(...
            MtX*fi2Nmk(dtnM1,ubrkn,ubprkn,ubpprkn) + CtX*fi1Nmk(dtnM1,ubrkn,ubprkn,ubpprkn));
        %Variável de integração (znM1 = Auz*ubrknM1 + Almbz*lmb):
        IdrkM1 = eye(cPk+1);
        Auz = IdrkM1(:,1:cPk);
        Almbz = IdrkM1(:,cPk+1);
        %Equação de integração: 
        %Obs.: Knmk(dtnM1,znM1)*znM1 + fiTilrkzn(dtnM1,tnM1) - Bkz*TxnM1 = 0
        Knmk = @(dtnM1,znM1)(...
            [KhatTx(dtnM1,Auz.'*znM1), kNmk*GamaRrk.'; 
             kNmk*GamaRrk, 0]);
        fiTilrkzn = @(dtnM1,tnM1,ubrkn,ubprkn,ubpprkn)(...
            [fiTilrkn(dtnM1,ubrkn,ubprkn,ubpprkn); 
             -kNmk*hRrk(tnM1)]);
        Bkz = [BtX; 0];
    case 32
        IdNgmrEstr = eye(Ng-rEstr);
        Au1ur = IdNgmrEstr(:,1);
        AuSur = IdNgmrEstr(:,2:Ng-rEstr);
        M1s = (Pk.'*AuSur)*AuSur.'*M1restr*AuSur*(AuSur.'*Pk);
        M1s = (1/2)*(M1s + M1s.');
        Cs = (Pk.'*AuSur)*AuSur.'*Crestr*AuSur*(AuSur.'*Pk);
        Cs = (1/2)*(Cs + Cs.');
        Ks = (Pk.'*AuSur)*AuSur.'*Krestr*AuSur*(AuSur.'*Pk);
        Ks = (1/2)*(Ks + Ks.');
        fubEFr = @(ub1,ubs)(Au1ur*ub1 + AuSur*ubs);
        fKls = @(ub1,ubs)((Pk.'*AuSur)*AuSur.'*fKlr(fubEFr(ub1,ubs))*AuSur*(AuSur.'*Pk));
        fKls = @(ub1,ubs)((1/2)*(fKls(ub1,ubs) + fKls(ub1,ubs).'));
        fKs = @(ub1,ubrk)(Ks + fKls(ub1,(AuSur.'*Pk)*ubrk));
        Bs = (Pk.'*AuSur)*AuSur.'*Br(:,1);
        
        M1s1 = (Pk.'*AuSur)*AuSur.'*M1restr*Au1ur;
        Cs1 = (Pk.'*AuSur)*AuSur.'*Crestr*Au1ur;
        Ks1 = (Pk.'*AuSur)*AuSur.'*Krestr*Au1ur;
        fKls1 = @(ub1,ubs)((Pk.'*AuSur)*AuSur.'*fKlr(fubEFr(ub1,ubs))*Au1ur);
        fKs1 = @(ub1,ubrk)(Ks1 + fKls1(ub1,(AuSur.'*Pk)*ubrk));
        
        yp2dAj = @(t)(dNr0bAjdw(y1d(t))*y1pd(t));
        hUp = @(t)(yp2dAj(t)/kNM1);
        ypp2dAj = @(t)(d2Nr0bAjdw2(y1d(t))*y1pd(t) + dNr0bAjdw(y1d(t))*y1ppd(t));
        hUpp = @(t)(ypp2dAj(t)/kNM1);
        
        %Entrada:
        Tref = 10; %Nm
        fTx = @(t)(Tref);
        fub1IN = @(t)(hU(t));
        fubp1IN = @(t)(hUp(t));
        fubpp1IN = @(t)(hUpp(t));
        %Formução de Newmark:
        Khats = @(dtnM1,tnM1,ubrknM1)(alfa2Nmk(dtnM1)*M1s + alfa1Nmk(dtnM1)*Cs + fKs(fub1IN(tnM1),ubrknM1));
        fiTilskn = @(dtnM1,ubrkn,ubprkn,ubpprkn)(...
            M1s*fi2Nmk(dtnM1,ubrkn,ubprkn,ubpprkn) + Cs*fi1Nmk(dtnM1,ubrkn,ubprkn,ubpprkn));
        
        %Teste (forma normal com restrições intrínsecas):
        %
        ubs0 = [zeros(Ng-rEstr-2,1); 1];
        Ks0 = Ks + fKls(0,ubs0);
        Ks10 = Ks1 + fKls1(0,ubs0);
        Ask = [zeros(cPk), eye(cPk); 
            -M1s\Ks0 -M1s\Cs];
        Bsk = [zeros(cPk,1), zeros(cPk,1), zeros(cPk,1), zeros(cPk,1); 
            M1s\Bs, -M1s\M1s1, -M1s\Cs1, -M1s\Ks10];
        
        alfab0 = CuNM1;
        alfab1 = alfab0*Ask;
        a21 = @(Zk)(alfab1*Zk);
        b211 = alfab0*Bsk(:,1);
        b221 = alfab0*Bsk(:,2);
        b231 = alfab0*Bsk(:,3);
        b241 = alfab0*Bsk(:,4);
        
        alfab2 = alfab1*Ask;
        a22 = @(Zk)(alfab2*Zk);
        b212 = alfab1*Bsk(:,1);
        b222 = alfab1*Bsk(:,2);
        b232 = alfab1*Bsk(:,3);
        b242 = alfab1*Bsk(:,4);
        
        alfab3 = alfab2*Ask;
        a23 = @(Zk)(alfab3*Zk);
        b213 = alfab2*Bsk(:,1);
        b223 = alfab2*Bsk(:,2);
        b233 = alfab2*Bsk(:,3);
        b243 = alfab2*Bsk(:,4);
        
        alfab4 = alfab3*Ask;
        a24 = @(Zk)(alfab4*Zk);
        b214 = alfab3*Bsk(:,1);
        b224 = alfab3*Bsk(:,2);
        b234 = alfab3*Bsk(:,3);
        b244 = alfab3*Bsk(:,4);
        
        %Equação: [y2, y2p, y2p2, y2p3, y2p4]^T = AlfaY2Zk*Zk + BY2*[Tx, ubpp1, ubp1, ub1]^T
        AlfaY2Zk = [alfab0; alfab1; alfab2; alfab3; alfab4];
        BY2 = [0, 0, 0, 0; 
               b211, b221, b231, b241;
               b212, b222, b232, b242; 
               b213, b223, b233, b243; 
               b214, b224, b234, b244];
        
        %{
        aaa = [b211, b221, b231, b241, b212, b222, b232, b242, b213, b223, b233, b243, b214, b224, b234, b244];
        %}
        %ad = 1;
        %}
    case 33
        %Entrada:
        Tref = 10; %Nm
        fTx = @(t)(Tref);
        repositorio.fTx = fTx;
        %Modelo alternativo em PC (Tat aprox Treacao):
        KtTr = [ K1mm,      -K1mm,     0; 
                -K1mm,  K1mm+K2mm, -K2mm; 
                 K1mm, -K2mm-K1mm,  K2mm];
        AtTr = [zeros(3), Id3; -Mt\KtTr, -Mt\Ct];
        argumentos.AtTr = AtTr;
        BbT0 = BbT;
        argumentos.BbT0 = BbT0;
        C10 = C1;
        argumentos.C10 = C10;
    case 34
        %Entrada:
        wRef = 10; %RPM
        wRef = wRef*pi/30;
        fTetap1 = @(t)(wRef*ones(length(t),1));
        fTeta1 = @(t)(wRef*t);
        fTetapp1 = @(t)(0*ones(length(t),1));
        %Quebra do vetor de posição:
        %Obs.: uTeta = Ateta1*teta1 + Ateta23*uteta23
        Ateta1 = Id3(:,1);
        Ateta23 = Id3(:,[2,3]);
        %Modelo alternativo em PC (Tat aprox Treacao):
        KtTr = [ K1mm,      -K1mm,     0; 
                -K1mm,  K1mm+K2mm, -K2mm; 
                 K1mm, -K2mm-K1mm,  K2mm];
        %Matrizes de inércia, rigidez e amortecimento:
        Mt23 = Ateta23.'*Mt*Ateta23;
        Ct23 = Ateta23.'*Ct*Ateta23;
        Kt23 = Ateta23.'*Kt*Ateta23;
        KtTr23 = Ateta23.'*KtTr*Ateta23;
        Mt231 = Ateta23.'*Mt*Ateta1;
        Ct231 = Ateta23.'*Ct*Ateta1;
        Kt231 = Ateta23.'*Kt*Ateta1;
        KtTr231 = Ateta23.'*KtTr*Ateta1;
        %Matriz de ganho e vetor de entrada:
        Bteta23 = [-Mt231, -Ct231, -Kt231];
        Teta1 = @(t)([fTetapp1(t); fTetap1(t); fTeta1(t)]);
        repositorio.Teta1 = Teta1;
        %Equações de estado:
        At23 = [zeros(2), Id2; -Mt23\Kt23, -Mt23\Ct23];
        argumentos.At23 = At23;
        Bbteta23 = [zeros(2,3); Bteta23];
        argumentos.Bbteta23 = Bbteta23;
        AtTr23 = [zeros(2), Id2; -Mt23\KtTr23, -Mt23\Ct23];
        argumentos.AtTr23 = AtTr23;
    case 35
        %Entradas (restrições):
        %hU = @(t)(0);
        %hlU = @(t)(y2d(t)/(kNM1*(nUaj+1)));
        %hU = @(t)(y2d(t)/kNM1);
        
        %Dinâmica longitudinal (modelo em PC):
        %Newmark:
        %Modelo em Elementos Finitos:
        %Obs.: xTil = AuEFx*ubTilEF + Amix*mid + Afix*fiBd
        IdNgM2 = eye(Ng+2);
        AuEFx = IdNgM2(:,1:Ng);
        Amix = IdNgM2(:,Ng+1);
        Afix = IdNgM2(:,Ng+2);
        
        M1gEFx = M1gEF*AuEFx.';
        CgEFx = CgEF*AuEFx.';
        KgEFx = KgEF*AuEFx.';
        
        %Newmark:
        %Termos lineares:
        KxNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEFx + alfa1Nmk(dtnM1)*CgEFx + KgEFx);
        fiTilxn = @(dtnM1,xn,xpn,xppn)(...
            M1gEFx*fi2Nmk(dtnM1,xn,xpn,xppn) + CgEFx*fi1Nmk(dtnM1,xn,xpn,xppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuEFx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuEFx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(AuEFx.'*xTil,AuEFx.'*xTilp));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fdcEFb(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuEFx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfdcEFbdubT(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuEFx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfdcEFbdubpT(AuEFx.'*xTil,AuEFx.'*xTilp));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuEFx.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(AuEFx.'*xTil));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Atx = AuEFx*Atu;
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        
        fubTilEF = @(xTil)(AuEFx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTuEF = @(ubTilEF)(Atu.'*ubTilEF);
        %fuUuEF = @(ubTilEF)(Auu.'*ubTilEF);
        %fvuEF = @(ubTilEF)(Avu.'*ubTilEF);
        %fwuEF = @(ubTilEF)(Awu.'*ubTilEF);
        %fuTxTil = @(xTil)(fuTuEF(fubTilEF(xTil)));
        %fuUxTil = @(xTil)(fuUuEF(fubTilEF(xTil)));
        
        %Matriz de ganho, na equação da dinâmica:
        BefTx = Bef(:,1);
        
        %Camada limite e dinâmica interna:
        %GamaMix = GamaMi*Amix.';
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT =  Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+2+rEstrT);
        Axz = IdNgMrM2(:,1:Ng+2);
        Almbz = IdNgMrM2(:,Ng+2+1:Ng+2+rEstrT);
        GamaMF = GamaRt;
        hTilMF = @(t)(fhRestr(t));
        GamaRz = Almbz*GamaMF;
        GamaRxz = GamaMF*(AuEFx.'*Axz.');
        Knmkxz = @(dtnM1)([KxNmk(dtnM1)*Axz.' + kNmk*GamaRz.'; 
            kNmk*GamaRxz]);
        fTilcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilMF(tnM1)]);
        fTildcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilMF(tnM1)]);
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        
        fmipdAtznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        %GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)([fiTilxn(dtnM1,xTiln,xTilpn,xTilppn); zeros(rEstrT,1)]);
        Befz = [BefTx; zeros(rEstrT,1)];
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
       dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        fdmipdAtdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdxnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        ad = 1;
    case 36
        %Entradas (restrições):
        %hU = @(t)(0);
        %hlU = @(t)(y2d(t)/(kNM1*(nUaj+1)));
        %hU = @(t)(y2d(t)/kNM1);
        
        %Dinâmica longitudinal (modelo em PC):
        %Newmark:
        %Modelo em Elementos Finitos:
        %Obs.: xTil = AuEFx*ubTilEF + Amix*mid + Afix*fiBd
        IdNgM2 = eye(Ng+2);
        AuEFx = IdNgM2(:,1:Ng);
        Amix = IdNgM2(:,Ng+1);
        Afix = IdNgM2(:,Ng+2);
        
        M1gEFx = M1gEF*AuEFx.';
        CgEFx = CgEF*AuEFx.';
        KgEFx = KgEF*AuEFx.';
        
        %Newmark:
        %Termos lineares:
        KxNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1gEFx + alfa1Nmk(dtnM1)*CgEFx + KgEFx);
        fiTilxn = @(dtnM1,xn,xpn,xppn)(...
            M1gEFx*fi2Nmk(dtnM1,xn,xpn,xppn) + CgEFx*fi1Nmk(dtnM1,xn,xpn,xppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuEFx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuEFx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(AuEFx.'*xTil,AuEFx.'*xTilp));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fdcEFb(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuEFx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfdcEFbdubT(AuEFx.'*xTil,AuEFx.'*xTilp,AuEFx.'*xTilpp));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuEFx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfdcEFbdubpT(AuEFx.'*xTil,AuEFx.'*xTilp));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuEFx.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(AuEFx.'*xTil));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Atx = AuEFx*Atu;
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        
        fubTilEF = @(xTil)(AuEFx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTuEF = @(ubTilEF)(Atu.'*ubTilEF);
        %fuUuEF = @(ubTilEF)(Auu.'*ubTilEF);
        %fvuEF = @(ubTilEF)(Avu.'*ubTilEF);
        %fwuEF = @(ubTilEF)(Awu.'*ubTilEF);
        %fuTxTil = @(xTil)(fuTuEF(fubTilEF(xTil)));
        %fuUxTil = @(xTil)(fuUuEF(fubTilEF(xTil)));
        
        %Matriz de ganho, na equação da dinâmica:
        BefTx = Bef(:,1);
        
        %Camada limite e dinâmica interna:
        %GamaMix = GamaMi*Amix.';
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT =  Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+2+rEstrT);
        Axz = IdNgMrM2(:,1:Ng+2);
        Almbz = IdNgMrM2(:,Ng+2+1:Ng+2+rEstrT);
        GamaMF = GamaRt;
        hTilMF = @(t)(fhRestr(t));
        GamaRz = Almbz*GamaMF;
        GamaRxz = GamaMF*(AuEFx.'*Axz.');
        Knmkxz = @(dtnM1)([KxNmk(dtnM1)*Axz.' + kNmk*GamaRz.'; 
            kNmk*GamaRxz]);
        fTilcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilMF(tnM1)]);
        fTildcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilMF(tnM1)]);
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        
        fmipdAtznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        %GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)([fiTilxn(dtnM1,xTiln,xTilpn,xTilppn); zeros(rEstrT,1)]);
        Befz = [BefTx; zeros(rEstrT,1)];
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
       dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        fdmipdAtdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdxnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        ad = 1;
    case 37
        %Entradas (restrições):
        %hU = @(t)(0);
        %hlU = @(t)(y2d(t)/(kNM1*(nUaj+1)));
        %hU = @(t)(y2d(t)/kNM1);
        
        %Dinâmica longitudinal (modelo em PC):
        %Newmark:
        %Obs.: zlnM1 = Aulz*ubl + Almblz*lmbl
        Aulz = Id4(:,1:3);
        Almblz = Id4(:,4);
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        GamaRu = Id3(1,:);
        fiuNmkn = @(dtnM1,ubn,ubpn,ubppn)(Mu*fi2Nmk(dtnM1,ubn,ubpn,ubppn) +...
            Cu*fi1Nmk(dtnM1,ubn,ubpn,ubppn));
        
        KuzNmk = @(dtnM1)([KuNmk(dtnM1), kNmk*GamaRu.'; kNmk*GamaRu, 0]);
        fiuNmkzn = @(dtnM1,ubn,ubpn,ubppn)([fiuNmkn(dtnM1,ubn,ubpn,ubppn); 0]);
        FnzNmk = @(tnM1,zlnM1)([fNuTil(Aulz.'*zlnM1); 
            -kNmk*hU(tnM1)]);
        dFnzdzT = @(zlnM1)([fdNudubTilT(Aulz.'*zlnM1), zeros(3,1); 
            zeros(1,3), 0]);
        %Ayx = Id6(:,1:5);
        
        %Modelo em Elementos Finitos:
        %Obs.: xTil = AuEFx*ubTilEF + Amix*mid + Afix*fiBd
        Id8M2 = eye(8+2);
        AuAcx = Id8M2(:,1:8);
        Amix = Id8M2(:,8+1);
        Afix = Id8M2(:,8+2);
        
        M1acx = M1ac*AuAcx.';
        Cacx = Cac*AuAcx.';
        Kacx = Kac*AuAcx.';
        
        %Newmark:
        %Termos lineares:
        KxNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1acx + alfa1Nmk(dtnM1)*Cacx + Kacx);
        fiTilxn = @(dtnM1,xn,xpn,xppn)(...
            M1acx*fi2Nmk(dtnM1,xn,xpn,xppn) + Cacx*fi1Nmk(dtnM1,xn,xpn,xppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fTilc(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfTilcdubT(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfTilcdubpT(AuAcx.'*xTil,AuAcx.'*xTilp));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fTildc(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfTildcdubT(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfTildcdubpT(AuAcx.'*xTil,AuAcx.'*xTilp));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuAcx.';
        dfTildubppTx = @(xTil)(dfTildubppT(AuAcx.'*xTil));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = Id8(:,[2,6,8]);
        Auu = Id8(:,[1,3,7]);
        Atx = AuAcx*Atu;
        Avu = Id8(:,4);
        Awu = Id8(:,5);
        
        fubTilAc = @(xTil)(AuAcx.'*xTil);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        fuTuEF = @(ubTilAc)(Atu.'*ubTilAc);
        
        %Matriz de ganho, na equação da dinâmica:
        BacTx = Bac(:,1);
        
        %Camada limite e dinâmica interna:
        %GamaMix = GamaMi*Amix.';
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT =  Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        Id8MrM2 = eye(8+2+1);
        Axz = Id8MrM2(:,1:8+2);
        Almbz = Id8MrM2(:,8+2+1);
        GamaAc = GamaRub1Ac;
        %hTilAc = @(t)(hlU(t));
        hTilAc = @(t)(hU(t));
        GamaRz = Almbz*GamaAc;
        GamaRxz = GamaAc*(AuAcx.'*Axz.');
        Knmkxz = @(dtnM1)([KxNmk(dtnM1)*Axz.' + kNmk*GamaRz.'; 
            kNmk*GamaRxz]);
        fTilcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilAc(tnM1)]);
        fTildcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilAc(tnM1)]);
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        
        fmipdAtznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        %GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)([fiTilxn(dtnM1,xTiln,xTilpn,xTilppn); zeros(1,1)]);
        Bacz = [BacTx; zeros(1,1)];
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(1,8+2+1)]);
        dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(1,8+2+1)]);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        fdmipdAtdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdxnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
    case 38
        %Entradas (restrições):
        %hU = @(t)(0);
        %hlU = @(t)(y2d(t)/(kNM1*(nUaj+1)));
        %hU = @(t)(y2d(t)/kNM1);
        
        %Matriz de ganho, na equação da dinâmica:
        BefTx = Bef(:,1);
        
        %Dinâmica longitudinal (modelo em PC):
        %Newmark:
        %Obs.: zlnM1 = Aulz*ubl + Almblz*lmbl
        Aulz = Id4(:,1:3);
        Almblz = Id4(:,4);
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        GamaRu = Id3(1,:);
        fiuNmkn = @(dtnM1,ubn,ubpn,ubppn)(Mu*fi2Nmk(dtnM1,ubn,ubpn,ubppn) +...
            Cu*fi1Nmk(dtnM1,ubn,ubpn,ubppn));
        
        KuzNmk = @(dtnM1)([KuNmk(dtnM1), kNmk*GamaRu.'; kNmk*GamaRu, 0]);
        fiuNmkzn = @(dtnM1,ubn,ubpn,ubppn)([fiuNmkn(dtnM1,ubn,ubpn,ubppn); 0]);
        FnzNmk = @(tnM1,zlnM1)([fNuTil(Aulz.'*zlnM1); 
            -kNmk*hU(tnM1)]);
        dFnzdzT = @(zlnM1)([fdNudubTilT(Aulz.'*zlnM1), zeros(3,1); 
            zeros(1,3), 0]);
        %Ayx = Id6(:,1:5);
        
        %Modelo em Elementos Finitos:
        %Obs.: xTil = AuEFx*ubTilEF + Amix*mid + Afix*fiBd
        IdNgM2 = eye(Ng+2);
        AuAcx = IdNgM2(:,1:Ng);
        Amix = IdNgM2(:,Ng+1);
        Afix = IdNgM2(:,Ng+2);
        
        M1acx = M1gEF*AuAcx.';
        Cacx = CgEF*AuAcx.';
        Kacx = KgEF*AuAcx.';
        
        %Newmark:
        %Termos lineares:
        KxNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1acx + alfa1Nmk(dtnM1)*Cacx + Kacx);
        fiTilxn = @(dtnM1,xn,xpn,xppn)(...
            M1acx*fi2Nmk(dtnM1,xn,xpn,xppn) + Cacx*fi1Nmk(dtnM1,xn,xpn,xppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(AuAcx.'*xTil,AuAcx.'*xTilp));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fdcEFb(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfdcEFbdubT(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfdcEFbdubpT(AuAcx.'*xTil,AuAcx.'*xTilp));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuAcx.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(AuAcx.'*xTil));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Atx = AuAcx*Atu;
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        
        fubTilAc = @(xTil)(AuAcx.'*xTil);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        %fuTuEF = @(ubTilAc)(Atu.'*ubTilAc);
        
        %Matriz de ganho, na equação da dinâmica:
        BacTx = Bef(:,1);
        
        %Camada limite e dinâmica interna:
        %GamaMix = GamaMi*Amix.';
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT =  Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+2+rEstrT);
        Axz = IdNgMrM2(:,1:Ng+2);
        Almbz = IdNgMrM2(:,Ng+2+1:Ng+2+rEstrT);
        GamaAc = GamaRt;
        hTilAc = @(t)(fhRestr(t));
        GamaRz = Almbz*GamaAc;
        GamaRxz = GamaAc*(AuAcx.'*Axz.');
        Knmkxz = @(dtnM1)([KxNmk(dtnM1)*Axz.' + kNmk*GamaRz.'; 
            kNmk*GamaRxz]);
        fTilcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilAc(tnM1)]);
        fTildcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilAc(tnM1)]);
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        
        fmipdAtznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        %GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)([fiTilxn(dtnM1,xTiln,xTilpn,xTilppn); zeros(rEstrT,1)]);
        Bacz = [BefTx; zeros(rEstrT,1)];
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        fdmipdAtdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdxnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
    case 39
        %Janela de tempo para média movel de y1:
        %dt = argumentos.dt_newmark;
        %tJ = t09; %[s]
        %NJ = ceil(tJ/dt);
        
        %Matriz de ganho, na equação da dinâmica:
        BefTx = Bef(:,1);
        
        %Dinâmica longitudinal (modelo em PC):
        %Newmark:
        %Obs.: zlnM1 = Aulz*ubl + Almblz*lmbl
        Aulz = Id4(:,1:3);
        Almblz = Id4(:,4);
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        GamaRu = Id3(1,:);
        fiuNmkn = @(dtnM1,ubn,ubpn,ubppn)(Mu*fi2Nmk(dtnM1,ubn,ubpn,ubppn) +...
            Cu*fi1Nmk(dtnM1,ubn,ubpn,ubppn));
        
        KuzNmk = @(dtnM1)([KuNmk(dtnM1), kNmk*GamaRu.'; kNmk*GamaRu, 0]);
        fiuNmkzn = @(dtnM1,ubn,ubpn,ubppn)([fiuNmkn(dtnM1,ubn,ubpn,ubppn); 0]);
        FnzNmk = @(tnM1,zlnM1)([fNuTil(Aulz.'*zlnM1); 
            -kNmk*hUmm(tnM1)]);
        dFnzdzT = @(zlnM1)([fdNudubTilT(Aulz.'*zlnM1), zeros(3,1); 
            zeros(1,3), 0]);
        %Ayx = Id6(:,1:5);
        
        %Modelo em Elementos Finitos:
        %Obs.: xTil = AuEFx*ubTilEF + Amix*mid + Afix*fiBd
        IdNgM2 = eye(Ng+2);
        AuAcx = IdNgM2(:,1:Ng);
        Amix = IdNgM2(:,Ng+1);
        Afix = IdNgM2(:,Ng+2);
        
        M1acx = M1gEF*AuAcx.';
        Cacx = CgEF*AuAcx.';
        Kacx = KgEF*AuAcx.';
        
        %Newmark:
        %Termos lineares:
        KxNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1acx + alfa1Nmk(dtnM1)*Cacx + Kacx);
        fiTilxn = @(dtnM1,xn,xpn,xppn)(...
            M1acx*fi2Nmk(dtnM1,xn,xpn,xppn) + Cacx*fi1Nmk(dtnM1,xn,xpn,xppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(AuAcx.'*xTil,AuAcx.'*xTilp));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fdcEFb(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfdcEFbdubT(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfdcEFbdubpT(AuAcx.'*xTil,AuAcx.'*xTilp));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuAcx.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(AuAcx.'*xTil));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Atx = AuAcx*Atu;
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        
        fubTilAc = @(xTil)(AuAcx.'*xTil);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        %fuTuEF = @(ubTilAc)(Atu.'*ubTilAc);
        AtetaYu = IdNg(:,6*(N1+1)-1);
        AtetaZu = IdNg(:,6*(N1+1)-3);
        AtetaXu = IdNg(:,6*(N1+1));
        AtetaEuuAc = [AtetaYu, AtetaZu, AtetaXu];
        AtetaEux = AuAcx*AtetaEuuAc;
        
        %Matriz de ganho, na equação da dinâmica:
        BacTx = Bef(:,1);
        
        %Camada limite e dinâmica interna:
        %GamaMix = GamaMi*Amix.';
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT =  Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+2+rEstrT);
        Axz = IdNgMrM2(:,1:Ng+2);
        Almbz = IdNgMrM2(:,Ng+2+1:Ng+2+rEstrT);
        GamaAc = GamaRt;
        hTilAc = @(t)(fhRestr(t));
        GamaRz = Almbz*GamaAc;
        GamaRxz = GamaAc*(AuAcx.'*Axz.');
        Knmkxz = @(dtnM1)([KxNmk(dtnM1)*Axz.' + kNmk*GamaRz.'; 
            kNmk*GamaRxz]);
        fTilcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilAc(tnM1)]);
        fTildcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilAc(tnM1)]);
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        
        fmipdAtznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        %GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)([fiTilxn(dtnM1,xTiln,xTilpn,xTilppn); zeros(rEstrT,1)]);
        Bacz = [BefTx; zeros(rEstrT,1)];
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        fdmipdAtdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdxnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        
        %Ângulos de Euler para o retor intermediário:
        fiEu = @(FiTil)(Id3(1,:)*FiTil);
        tetaEu = @(FiTil)(Id3(2,:)*FiTil);
        alfaEu = @(FiTil)(Id3(3,:)*FiTil);
        
        uHat1eI = [1; 0; 0];
        uHat2eI = @(fi)([0; cos(fi); sin(fi)]);
        uHat3eI = @(fi,teta)([cos(teta); sin(fi)*sin(teta); -cos(fi)*sin(teta)]);
        uHat2eI = @(FiTil)(uHat2eI(fiEu(FiTil)));
        uHat3eI = @(FiTil)(uHat3eI(fiEu(FiTil),tetaEu(FiTil)));
        ReI = @(FiTil)([uHat1eI, uHat2eI(FiTil), uHat3eI(FiTil)]);
        
        %Ângulos de Tait–Bryan para o retor intermediário:
        tetaYc = @(TetaTil)(Id3(1,:)*TetaTil);
        tetaZc = @(TetaTil)(Id3(2,:)*TetaTil);
        tetaXc = @(TetaTil)(Id3(3,:)*TetaTil);
        
        uHat1cI = [0; 1; 0];
        uHat2cIaux = @(tetaY)([sin(tetaY); 0; cos(tetaY)]);
        uHat3cIaux = @(tetaY,tetaZ)([cos(tetaY)*cos(tetaZ); sin(tetaZ); -cos(tetaZ)*sin(tetaY)]);
        uHat2cI = @(TetaTil)(uHat2cIaux(tetaYc(TetaTil)));
        uHat3cI = @(TetaTil)(uHat3cIaux(tetaYc(TetaTil),tetaZc(TetaTil)));
        RcI = @(TetaTil)([uHat1cI, uHat2cI(TetaTil), uHat3cI(TetaTil)]);
        
        %Funções para integração:
        fWce0 = @(FiTil,FipTil,TetaTil,TetapTil)(uHat1eI*fiEu(FipTil) +...
            uHat2eI(FiTil)*tetaEu(FipTil) +...
            uHat3eI(FiTil)*alfaEu(FipTil) - RcI(TetaTil)*TetapTil);
        fFipTilnM1 = @(dtnM1,FiTilnM1,un,FipTiln,FippTiln)(...
            fupnM1(dtnM1,FiTilnM1,un,FipTiln,FippTiln));
        fWce0nM1 = @(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln,TetaTilnM1,TetapTilnM1)(...
            fWce0(FiTilnM1,fFipTilnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln),TetaTilnM1,TetapTilnM1));
        
        duHat1eIdFiT = zeros(3);
        duHat2eIdFiT = @(fi)([0, 0, 0; 
            -sin(fi), 0, 0; 
            cos(fi), 0, 0]);
        duHat3eIdFiT = @(fi,teta)([0, -sin(teta), 0; 
            cos(fi)*sin(teta), cos(teta)*sin(fi), 0; 
            sin(fi)*sin(teta), -cos(fi)*cos(teta), 0]);
        duHat2eIdFiT = @(FiTil)(duHat2eIdFiT(fiEu(FiTil)));
        duHat3eIdFiT = @(FiTil)(duHat3eIdFiT(fiEu(FiTil),tetaEu(FiTil)));
        
        dfipdFipT = Id3(1,:);
        dtetapdFipT = Id3(2,:);
        dalfapdFipT = Id3(3,:);
        dFipnM1dFinM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dfWce0dFinM1T = @(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln)(...
            duHat1eIdFiT*fiEu(fFipTilnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln)) +...
            duHat2eIdFiT(FiTilnM1)*tetaEu(fFipTilnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln)) +...
            duHat3eIdFiT(FiTilnM1)*alfaEu(fFipTilnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln)) +...
            (uHat1eI*dfipdFipT +...
             uHat2eI(FiTilnM1)*dtetapdFipT +...
             uHat3eI(FiTilnM1)*dalfapdFipT)*...
            dFipnM1dFinM1T(dtnM1));
        
        %Aceleração angular:
        duHat1cIdTetaT = zeros(3);
        duHat2cIdTetaT = @(tetaY)([cos(tetaY), 0, 0; 
            0, 0, 0;
            -sin(tetaY), 0, 0]);
        duHat3cIdTetaT = @(tetaY,tetaZ)([ -cos(tetaZ)*sin(tetaY), -cos(tetaY)*sin(tetaZ), 0;
            0, cos(tetaZ), 0;
            -cos(tetaY)*cos(tetaZ),  sin(tetaY)*sin(tetaZ), 0]);
        duHat2cIdTetaT = @(TetaTil)(duHat2cIdTetaT(tetaYc(TetaTil)));
        duHat3cIdTetaT = @(TetaTil)(duHat3cIdTetaT(tetaYc(TetaTil),tetaZc(TetaTil)));
        
        dReIdFiTFip = @(FiTil,FipTil)(...
            [duHat1eIdFiT*FipTil, duHat2eIdFiT(FiTil)*FipTil, duHat3eIdFiT(FiTil)*FipTil]);
        dRcIdTetaTTetap = @(TetaTil,TetapTil)(...
            [duHat1cIdTetaT*TetapTil, duHat2cIdTetaT(TetaTil)*TetapTil, duHat3cIdTetaT(TetaTil)*TetapTil]);
        
        fFipTil = @(FiTil,TetaTil,TetapTil)(ReI(FiTil)\(RcI(TetaTil)*TetapTil));
        dRcITetapdt = @(TetaTil,TetapTil,TetappTil)(...
            RcI(TetaTil)*TetappTil + dRcIdTetaTTetap(TetaTil,TetapTil)*TetapTil);
        fFippTil = @(FiTil,FipTil,TetaTil,TetapTil,TetappTil)(...
            ReI(FiTil)\(-dReIdFiTFip(FiTil,FipTil)*FipTil + dRcITetapdt(TetaTil,TetapTil,TetappTil)));
    case 40
        %Janela de tempo para média movel de y1:
        dt = argumentos.dt_newmark;
        tJ = t09; %[s]
        NJ = ceil(tJ/dt);
        
        %Matriz de ganho, na equação da dinâmica:
        BefTx = Bef(:,1);
        
        %Dinâmica longitudinal (modelo em PC):
        %Newmark:
        %Obs.: zlnM1 = Aulz*ubl + Almblz*lmbl
        Aulz = Id4(:,1:3);
        Almblz = Id4(:,4);
        KuNmk = @(dtnM1)(alfa2Nmk(dtnM1)*Mu + alfa1Nmk(dtnM1)*Cu + Ku);
        GamaRu = Id3(1,:);
        fiuNmkn = @(dtnM1,ubn,ubpn,ubppn)(Mu*fi2Nmk(dtnM1,ubn,ubpn,ubppn) +...
            Cu*fi1Nmk(dtnM1,ubn,ubpn,ubppn));
        
        KuzNmk = @(dtnM1)([KuNmk(dtnM1), kNmk*GamaRu.'; kNmk*GamaRu, 0]);
        fiuNmkzn = @(dtnM1,ubn,ubpn,ubppn)([fiuNmkn(dtnM1,ubn,ubpn,ubppn); 0]);
        FnzNmk = @(tnM1,zlnM1,yb1nM1)([fNuTil(Aulz.'*zlnM1); 
            -kNmk*hUmmAst(tnM1,yb1nM1)]);
        dFnzdzT = @(zlnM1)([fdNudubTilT(Aulz.'*zlnM1), zeros(3,1); 
            zeros(1,3), 0]);
        %Ayx = Id6(:,1:5);
        
        %Modelo em Elementos Finitos:
        %Obs.: xTil = AuEFx*ubTilEF + Amix*mid + Afix*fiBd
        IdNgM2 = eye(Ng+2);
        AuAcx = IdNgM2(:,1:Ng);
        Amix = IdNgM2(:,Ng+1);
        Afix = IdNgM2(:,Ng+2);
        
        M1acx = M1gEF*AuAcx.';
        Cacx = CgEF*AuAcx.';
        Kacx = KgEF*AuAcx.';
        
        %Newmark:
        %Termos lineares:
        KxNmk = @(dtnM1)(alfa2Nmk(dtnM1)*M1acx + alfa1Nmk(dtnM1)*Cacx + Kacx);
        fiTilxn = @(dtnM1,xn,xpn,xppn)(...
            M1acx*fi2Nmk(dtnM1,xn,xpn,xppn) + Cacx*fi1Nmk(dtnM1,xn,xpn,xppn));
        
        %Termos não lineares:
        %Modelos contínuos de atrito:
        fTilcx = @(xTil,xTilp,xTilpp)(fcEFb(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        fTilcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTilcdubTx = @(xTil,xTilp,xTilpp)(...
            dfcEFbdubT(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        dfTilcdxTilT = @(xTil,xTilp,xTilpp)(dfTilcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTilcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilT(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTilcdubpTx = @(xTil,xTilp)(...
            dfcEFbdubpT(AuAcx.'*xTil,AuAcx.'*xTilp));
        dfTilcdxTilpT = @(xTil,xTilp)(dfTilcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTilcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTilpT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Modelos descontínuos de atrito:
        fTildcx = @(xTil,xTilp,xTilpp)(fdcEFb(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        fTildcxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fTildcx(xTilnM1,...
                    fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                    fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubTildxT = AuAcx.';
        dfTildcdubTx = @(xTil,xTilp,xTilpp)(...
            dfdcEFbdubT(AuAcx.'*xTil,AuAcx.'*xTilp,AuAcx.'*xTilpp));
        dfTildcdxTilT = @(xTil,xTilp,xTilpp)(dfTildcdubTx(xTil,xTilp,xTilpp)*dubTildxT);
        dfTildcdxTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilT(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        dubpTildxpT = AuAcx.';
        dfTildcdubpTx = @(xTil,xTilp)(...
            dfdcEFbdubpT(AuAcx.'*xTil,AuAcx.'*xTilp));
        dfTildcdxTilpT = @(xTil,xTilp)(dfTildcdubpTx(xTil,xTilp)*dubpTildxpT);
        dfTildcdxpTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTilpT(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        %Derivadas com relação ao vetor aceleração:
        dubppTildxppT = AuAcx.';
        dfTildubppTx = @(xTil)(dfEFbdubppT(AuAcx.'*xTil));
        dfTildxTilppT = @(xTil)(dfTildubppTx(xTil)*dubppTildxppT);
        dfTildxppTnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildxTilppT(xTilnM1));        
        %Derivadas totais:
        dxpnM1dxnM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(alfa2Nmk(dtnM1));
        DfTilcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTilcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTilcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        DfTildcDxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dfTildcdxTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn) +...
            dfTildcdxpTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxpnM1dxnM1T(dtnM1) +...
            dfTildxppTnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)*dxppnM1dxnM1T(dtnM1));
        
        %Obs.: ubTilEF = Atu*uTeta + Auu*ubU + Avu*vL + Awu*wL +...
        Atu = IdNg(:,[6*1,6*(N1+1),6*(N+1)]);
        Auu = IdNg(:,[6*1-5,6*(N1+1)-5,6*(N+1)-5]);
        Atx = AuAcx*Atu;
        Avu = IdNg(:,6*(N1+1)-4);
        Awu = IdNg(:,6*(N1+1)-2);
        
        fubTilAc = @(xTil)(AuAcx.'*xTil);
        ftx = @(xTil)(Atx.'*xTil);
        fmidx = @(xTil)(Amix.'*xTil);
        ffiBdTx = @(xTil)(Afix.'*xTil);
        %fuTuEF = @(ubTilAc)(Atu.'*ubTilAc);
        AtetaYu = IdNg(:,6*(N1+1)-1);
        AtetaZu = IdNg(:,6*(N1+1)-3);
        AtetaXu = IdNg(:,6*(N1+1));
        AtetaEuuAc = [AtetaYu, AtetaZu, AtetaXu];
        AtetaEux = AuAcx*AtetaEuuAc;
        
        %Matriz de ganho, na equação da dinâmica:
        BacTx = Bef(:,1);
        
        %Camada limite e dinâmica interna:
        %GamaMix = GamaMi*Amix.';
        fmipdAtx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fmipdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fmipdAtxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        
        fmidpx = @(xTilp)(Amix.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        ffipBdTx = @(xTilp)(Afix.'*xTilp);
        GamaFiBdTx = @(xTilp)(GamaFiBdT(ffipBdTx(xTilp)));
        GamaFiBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        
        KfiTx = KfiT*Afix.';
        
        KlcTdx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            KlcTdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        KlcTdxnM1 = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dupTildxpT =  Atx.';
        dxTilpdxTilT = @(dtnM1)(alfa1Nmk(dtnM1));
        dFatTcduTilTx = @(xTilp,y2)(dFatTcduTilT(ftx(xTilp),y2));
        dFatTcdxTilpT = @(xTilp,y2)(dFatTcduTilTx(xTilp,y2)*dupTildxpT);
        DFatTcDxTilTnM1 = @(dtnM1,xTilnM1,y2nM1,xTiln,xTilpn,xTilppn)(...
            dFatTcdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1)*...
            dxTilpdxTilT(dtnM1));
        
        dmiddxTilT = Amix.';
        fdmipdAtdmidx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(...
            fdmipdAtdmid(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        fdmipdAtdxnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdmidx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dmiddxTilT);
        
        dmipddxTilpT = Amix.';
        dxTilpdxTilnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dmipddxTilnM1T = @(dtnM1)(dmipddxTilpT*dxTilpdxTilnM1T(dtnM1));
        
        dfipBddxTilpT = Afix.';
        dGamaFiTdfipBdTx = @(xTilp)(dGamaFiTdfipBdT(ffipBdTx(xTilp)));
        dGamaFiTxdxTilpT = @(xTilp)(dGamaFiTdfipBdTx(xTilp)*dfipBddxTilpT);
        dGamaFiTxdxTilnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*...
            dxTilpdxTilnM1T(dtnM1));
        
        dKlcTddmidTx = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(fdKlcTddmidT(t,fmidx(xTil),y2,yp2,y2p2,y3p2,y4p2,Z0));
        dKlcTdxdxTilT = @(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(dKlcTddmidTx(t,xTil,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxTilT);
        dKlcTdxdxTilnM1T = @(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilT(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        %Vetor de estado desejado (cinemática inversa):
        fZtdx = @(t,xTil,y2,yp2,y2p2,y3p2,Z0)(fZtdAt(t,fmidx(xTil),y2,yp2,y2p2,y3p2,Z0));
        
        %Vetores de estado:
        fZx = @(xTil,xpTil)(fZ(Atx.'*xTil,Atx.'*xpTil));
        fZxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            fZx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fXx = @(xTil,xpTil,y2,yp2,y2p2,y3p2)(fXtAt(Atx.'*xTil,Atx.'*xpTil,y2,yp2,y2p2,y3p2));
        fXxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            fXx(xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1));
        
        %Variável de escorregamento:
        s1uTil = @(t,uTil,upTil,y2,yp2,y2p2,y3p2)(s1(t,fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2)));
        s1xnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)(...
            s1(tnM1,fXxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1)));
        
        %Difeomorfismo para a dinâmica longitudinal:
        QsiDifUu = @(uTil,upTil)(QsiDifU(fZ(uTil,upTil)));
        
        %Lei de controle:
        UcT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(...
            UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2) -...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        UcTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            UcTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        UcTxnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        duTildxT = Atx.';
        dmiddxT = Amix.';
        dfiBdTdxT = Afix.';
        dUcThatdZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(-da1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dy1dXT = Id6(1,:);
        dy1pdXT = Id6(2,:);
        dy1ppdXT = Id6(3,:);
        dy1pppdXT = Id6(4,:);
        dy1ppppdXT = Id6(5,:);
        dy5r1dXT = -4*lambda1*dy1ppppdXT - 6*lambda1^2*dy1pppdXT -...
            4*lambda1^3*dy1ppdXT - lambda1^4*dy1pdXT;
        dUcThatdXT = dy5r1dXT;
        dabsda1fndZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dda1fnAtdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*...
            sign(da1fnAt(Z,y2,yp2,y2p2,y3p2,y4p2)));
        dFlc1dZT = @(Z,y2,yp2,y2p2,y3p2,y4p2)(dabsda1fndZT(Z,y2,yp2,y2p2,y3p2,y4p2));
        dabsUcThatdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(...
            dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2)*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(...
            dFlc1dZT(Z,y2,yp2,y2p2,y3p2,y4p2) + D11fn*dabsUcThatdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dabsUcThatdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dUcThatdXT*sign(UcThatAt(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKlcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)((1 - D11fn)\(D11fn*dabsUcThatdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)));
        dKbLcTdZT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dKbLcTdXT = @(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)(dKlcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2));
        dUcTuTildZT = @(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2)(...
            b1fn\(dUcThatdZT(Z,y2,yp2,y2p2,y3p2,y4p2) - dKbLcTdZT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd)));
        dy4r1dXT = -4*lambda1*dy1pppdXT - 6*lambda1^2*dy1ppdXT -...
            4*lambda1^3*dy1pdXT - lambda1^4*dy1dXT;
        ds1dXT = dy1ppppdXT - dy4r1dXT;
        dsatS1FibddXT = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(ds1dXT/fiBd));
        dUcTuTildXT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(dUcThatdXT -...
            (dKbLcTdXT(t,Z,X,y2,yp2,y2p2,y3p2,y4p2)*sat(s1(t,X)/fiBd) +...
             KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatS1FibddXT(t,X,fiBd))));
        fdXtdZtT = @(Z,y2,yp2,y2p2,y3p2)(QsiDifTx(Z,y2,yp2,y2p2,y3p2));
        DUcTuTilDZT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildZT(t,Z,X,fiBd,y2,yp2,y2p2,y3p2,y4p2) +...
            dUcTuTildXT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*fdXtdZtT(Z,y2,yp2,y2p2,y3p2));
        
        dZduTilT = AuTil;
        dZdupTilT = AvTil;
        dUcTuTilduTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZduTilT);
        dUcTuTildupTilT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            DUcTuTilDZT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dZdupTilT);
        dXddmid = AdT;
        dZddmid = @(t,mid,y2,yp2,y2p2,y3p2,Z0)(diQsiDifTxdXT(fXtd(t,mid),y2,yp2,y2p2,y3p2,Z0)*dXddmid);
        dUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            -da1fnAtdZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0));
        dabsUcTdHatdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sign(UcTdHatAt(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKlcTddmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)((1 - D11fn)\(...
            dFlc1dZT(fZtdAt(t,mid,y2,yp2,y2p2,y3p2,Z0),y2,yp2,y2p2,y3p2,y4p2)*...
            dZddmid(t,mid,y2,yp2,y2p2,y3p2,Z0) +...
            D11fn*dabsUcTdHatdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)));
        dKbLcTdmid = @(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)(-dKlcTddmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTuTildmid = @(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-dKbLcTdmid(t,mid,y2,yp2,y2p2,y3p2,y4p2,Z0)*sat(s1(t,X)/fiBd)));
        dKbLcTdfiBd = KfiT;
        dsatSFibdTdfiBd = @(t,X,fiBd)(dsatdnu(s1(t,X)/fiBd)*(-s1(t,X)/fiBd^2));
        dUcTuTildfiBd = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            b1fn\(-(dKbLcTdfiBd*sat(s1(t,X)/fiBd) +...
            KbLcTat(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dsatSFibdTdfiBd(t,X,fiBd))));
        dUcTxdxT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTilduTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*duTildxT +...
            dUcTuTildmid(t,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dmiddxT +...
            dUcTuTildfiBd(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dfiBdTdxT);
        dUcTxdxTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        
        dUcTxdxpT = @(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTuTildupTilT(t,Z,X,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)*dupTildxpT);
        dUcTxdxpTuTil = @(t,uTil,upTil,mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpT(t,fZ(uTil,upTil),fXtAt(uTil,upTil,y2,yp2,y2p2,y3p2),mid,fiBd,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTx = @(t,xTil,xpTil,y2,yp2,y2p2,y3p2,y4p2,Z0)(...
            dUcTxdxpTuTil(t,Atx.'*xTil,Atx.'*xpTil,Amix.'*xTil,Afix.'*xTil,y2,yp2,y2p2,y3p2,y4p2,Z0));
        dUcTxdxpTnM1 = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxpTx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        dxpdxT = @(dtnM1)(alfa1Nmk(dtnM1));
        DUcTxDxnM1T = @(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dUcTxdxTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0) +...
            dUcTxdxpTnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxpdxT(dtnM1));
        
        %Restrições extrínsecas (multiplicadores de Lagrange):
        IdNgMrM2 = eye(Ng+2+rEstrT);
        Axz = IdNgMrM2(:,1:Ng+2);
        Almbz = IdNgMrM2(:,Ng+2+1:Ng+2+rEstrT);
        GamaAc = GamaRt;
        hTilAc = @(t)(fhRestr(t));
        GamaRz = Almbz*GamaAc;
        GamaRxz = GamaAc*(AuAcx.'*Axz.');
        Knmkxz = @(dtnM1)([KxNmk(dtnM1)*Axz.' + kNmk*GamaRz.'; 
            kNmk*GamaRxz]);
        fTilcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTilcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilAc(tnM1)]);
        fTildcxznM1 = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [fTildcxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn); 
             -kNmk*hTilAc(tnM1)]);
        UcTxznM1 = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            UcTxnM1(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fmipdznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            fmidpxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        
        fmipdAtznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fmipdAtxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0));
        %GamaMixz = GamaMix*Axz.';
        GamaFiBdTxznM1 = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdTxnM1(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn));
        KfiTxz = KfiTx*Axz.';
        KlcTdxznM1 = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            KlcTdxnM1(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0));
        fiNmkxzn = @(dtnM1,xTiln,xTilpn,xTilppn)([fiTilxn(dtnM1,xTiln,xTilpn,xTilppn); zeros(rEstrT,1)]);
        Bacz = [BefTx; zeros(rEstrT,1)];
        
        dxTildzT = Axz.';
        dfTilczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTilcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dfTildczdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            [DfTildcDxnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT; 
             zeros(rEstrT,Ng+2+rEstrT)]);
        dUcTxzdzT = @(tnM1,dtnM1,znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            DUcTxDxnM1T(tnM1,dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        dmipddzT = @(dtnM1)(dmipddxTilnM1T(dtnM1)*dxTildzT);
        fdmipdAtdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)(...
            fdmipdAtdxnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,Z0)*dxTildzT);
        dGamaFiTxzdzT = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamaFiTxdxTilnM1T(dtnM1,Axz.'*znM1,xTiln,xTilpn,xTilppn)*dxTildzT);
        dKlcTdxzdzT = @(tnM1,znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)(...
            dKlcTdxdxTilnM1T(tnM1,Axz.'*znM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,Z0)*...
            dxTildzT);
        
        %Ângulos de Euler para o retor intermediário:
        fiEu = @(FiTil)(Id3(1,:)*FiTil);
        tetaEu = @(FiTil)(Id3(2,:)*FiTil);
        alfaEu = @(FiTil)(Id3(3,:)*FiTil);
        
        uHat1eI = [1; 0; 0];
        uHat2eI = @(fi)([0; cos(fi); sin(fi)]);
        uHat3eI = @(fi,teta)([cos(teta); sin(fi)*sin(teta); -cos(fi)*sin(teta)]);
        uHat2eI = @(FiTil)(uHat2eI(fiEu(FiTil)));
        uHat3eI = @(FiTil)(uHat3eI(fiEu(FiTil),tetaEu(FiTil)));
        ReI = @(FiTil)([uHat1eI, uHat2eI(FiTil), uHat3eI(FiTil)]);
        
        %Ângulos de Tait–Bryan para o retor intermediário:
        tetaYc = @(TetaTil)(Id3(1,:)*TetaTil);
        tetaZc = @(TetaTil)(Id3(2,:)*TetaTil);
        tetaXc = @(TetaTil)(Id3(3,:)*TetaTil);
        
        uHat1cI = [0; 1; 0];
        uHat2cIaux = @(tetaY)([sin(tetaY); 0; cos(tetaY)]);
        uHat3cIaux = @(tetaY,tetaZ)([cos(tetaY)*cos(tetaZ); sin(tetaZ); -cos(tetaZ)*sin(tetaY)]);
        uHat2cI = @(TetaTil)(uHat2cIaux(tetaYc(TetaTil)));
        uHat3cI = @(TetaTil)(uHat3cIaux(tetaYc(TetaTil),tetaZc(TetaTil)));
        RcI = @(TetaTil)([uHat1cI, uHat2cI(TetaTil), uHat3cI(TetaTil)]);
        
        %Funções para integração:
        fWce0 = @(FiTil,FipTil,TetaTil,TetapTil)(uHat1eI*fiEu(FipTil) +...
            uHat2eI(FiTil)*tetaEu(FipTil) +...
            uHat3eI(FiTil)*alfaEu(FipTil) - RcI(TetaTil)*TetapTil);
        fFipTilnM1 = @(dtnM1,FiTilnM1,un,FipTiln,FippTiln)(...
            fupnM1(dtnM1,FiTilnM1,un,FipTiln,FippTiln));
        fWce0nM1 = @(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln,TetaTilnM1,TetapTilnM1)(...
            fWce0(FiTilnM1,fFipTilnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln),TetaTilnM1,TetapTilnM1));
        
        duHat1eIdFiT = zeros(3);
        duHat2eIdFiT = @(fi)([0, 0, 0; 
            -sin(fi), 0, 0; 
            cos(fi), 0, 0]);
        duHat3eIdFiT = @(fi,teta)([0, -sin(teta), 0; 
            cos(fi)*sin(teta), cos(teta)*sin(fi), 0; 
            sin(fi)*sin(teta), -cos(fi)*cos(teta), 0]);
        duHat2eIdFiT = @(FiTil)(duHat2eIdFiT(fiEu(FiTil)));
        duHat3eIdFiT = @(FiTil)(duHat3eIdFiT(fiEu(FiTil),tetaEu(FiTil)));
        
        dfipdFipT = Id3(1,:);
        dtetapdFipT = Id3(2,:);
        dalfapdFipT = Id3(3,:);
        dFipnM1dFinM1T = @(dtnM1)(alfa1Nmk(dtnM1));
        dfWce0dFinM1T = @(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln)(...
            duHat1eIdFiT*fiEu(fFipTilnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln)) +...
            duHat2eIdFiT(FiTilnM1)*tetaEu(fFipTilnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln)) +...
            duHat3eIdFiT(FiTilnM1)*alfaEu(fFipTilnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln)) +...
            (uHat1eI*dfipdFipT +...
             uHat2eI(FiTilnM1)*dtetapdFipT +...
             uHat3eI(FiTilnM1)*dalfapdFipT)*...
            dFipnM1dFinM1T(dtnM1));
        
        %Aceleração angular:
        duHat1cIdTetaT = zeros(3);
        duHat2cIdTetaT = @(tetaY)([cos(tetaY), 0, 0; 
            0, 0, 0;
            -sin(tetaY), 0, 0]);
        duHat3cIdTetaT = @(tetaY,tetaZ)([ -cos(tetaZ)*sin(tetaY), -cos(tetaY)*sin(tetaZ), 0;
            0, cos(tetaZ), 0;
            -cos(tetaY)*cos(tetaZ),  sin(tetaY)*sin(tetaZ), 0]);
        duHat2cIdTetaT = @(TetaTil)(duHat2cIdTetaT(tetaYc(TetaTil)));
        duHat3cIdTetaT = @(TetaTil)(duHat3cIdTetaT(tetaYc(TetaTil),tetaZc(TetaTil)));
        
        dReIdFiTFip = @(FiTil,FipTil)(...
            [duHat1eIdFiT*FipTil, duHat2eIdFiT(FiTil)*FipTil, duHat3eIdFiT(FiTil)*FipTil]);
        dRcIdTetaTTetap = @(TetaTil,TetapTil)(...
            [duHat1cIdTetaT*TetapTil, duHat2cIdTetaT(TetaTil)*TetapTil, duHat3cIdTetaT(TetaTil)*TetapTil]);
        
        fFipTil = @(FiTil,TetaTil,TetapTil)(ReI(FiTil)\(RcI(TetaTil)*TetapTil));
        dRcITetapdt = @(TetaTil,TetapTil,TetappTil)(...
            RcI(TetaTil)*TetappTil + dRcIdTetaTTetap(TetaTil,TetapTil)*TetapTil);
        fFippTil = @(FiTil,FipTil,TetaTil,TetapTil,TetappTil)(...
            ReI(FiTil)\(-dReIdFiTFip(FiTil,FipTil)*FipTil + dRcITetapdt(TetaTil,TetapTil,TetappTil)));
    otherwise
        disp('other value')
end

%%
%Integração:
t0 = argumentos.t0;
dt = argumentos.dt_newmark;
T = argumentos.T;
nT = argumentos.nT;

tn = t0;
switch simulacao
    case 1
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        Tx0 = fTx(tn);
        uTilppT0 = Mt\(Bt*Tx0 - Ct*uTilpT0 - Kt*uTilT0);
        ZZ = zeros(nT,3*3+1);
        ZZ(1,:) = [uTilT0.', uTilpT0.', uTilppT0.', Tx0.'];
        uTilTn = uTilT0;
        uTilpTn = uTilpT0;
        uTilppTn = uTilppT0;
        %
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            uTilTnM10 = uTilTn + dtnM1*uTilpTn + 0.5*(dtnM1*uTilppTn)*dtnM1;
            %Integração:
            Eta0 = uTilTnM10;
            Psi = @(eta)(KtNmk(dtnM1)*eta - Bt*fTx(tnM1) +...
                Mt*fi2Nmk(dtnM1,uTilTn,uTilpTn,uTilppTn) +...
                Ct*fi1Nmk(dtnM1,uTilTn,uTilpTn,uTilppTn));
            dPsidEtaT = @(eta)(KtNmk(dtnM1));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            uTilTnM1 = Eta_conv;
            uTilpTnM1 = fupnM1(dtnM1,uTilTnM1,uTilTn,uTilpTn,uTilppTn);
            uTilppTnM1 = fuppnM1(dtnM1,uTilTnM1,uTilTn,uTilpTn,uTilppTn);
            TxnM1 = fTx(tnM1);
            %Limitação de variáveis cíclicas:
            if ind_LVC
                %Limitação das variáveis cíclicas:
                teta1 = Id3(1,:)*uTilTnM1;
                uTilTnM1 = uTilTnM1 - 2*pi*ones(3,1)*double(teta1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [uTilTnM1.', uTilpTnM1.', uTilppTnM1.', TxnM1.'];
            %Atualização:
            uTilTn = uTilTnM1;
            uTilpTn = uTilpTnM1;
            uTilppTn = uTilppTnM1;
        end
        %}
        %Integração por ode45:
        if ODEativo
            odefun = @(t,Z)(At*Z + BbT*fTx(t));
            tspan = T;
            Z0 = zeros(6,1);
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [Tode,ZZode] = ode45(odefun,tspan,Z0,options);
        end
    case 2
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        fub1 = @(ubTil)(Id3(1,:)*ubTil);
        dub1bubTilT = Id3(1,:);
        Fbu0 = fFbu(tn,fub1(ubTilU0));
        dFbudub1 = -Kf;
        dFbudubTilT = dFbudub1*dub1bubTilT;
        ubTilppU0 = Mu\(Bu*Fbu0 - Cu*ubTilpU0 - Kbu*ubTilU0);
        ZZ = zeros(nT,3*3+1);
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', Fbu0.'];
        ZZpd = zeros(nT-nPD,3*3+1);
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            %Integração:
            Eta0 = ubTilUnM10;
            Psi = @(eta)(KuNmk(dtnM1)*eta - Bu*fFbu(tnM1,fub1(eta)) +...
                Mu*fi2Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn) +...
                Cu*fi1Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsidEtaT = @(eta)(KuNmk(dtnM1) - Bu*dFbudubTilT);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            ubTilUnM1 = Eta_conv;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ub1nM1 = fub1(ubTilUnM1);
            FbunM1 = fFbu(tnM1,ub1nM1);
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', FbunM1.'];
            if i>=nPD && MEANativo
                ZZpd(i-nPD+1,:) = mean(ZZ(i-nPD+1:i,:));
            end
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
        end
        %Integração por ode45:
        if ODEativo
            odefun = @(t,Z)(Au*Z + BbU*fFbu(t,fub1(fubTil(Z))));
            tspan = T;
            Z0 = zeros(6,1);
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [Tode,ZZode] = ode45(odefun,tspan,Z0,options);
        end
    case 3
        y1dTil0 = y1dTil(t0);
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        xTil0 = [mid0; fiBdT0];
        xTilp0 = [mipd0; fipBdT0];
        xTilpp0 = [mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        ZZ = zeros(nT,3*length(xTil0) + length(y1dTil0) + length(Zd0));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.'];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            Eta0 = xTilnM10;
            Psi = @(eta)([fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMix*eta - SigmaMi*y1dTil(tnM1);
                GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta)]);
            dPsidEtaT = @(eta)([dmipddxTilnM1T(dtnM1) - GamaMix;
                dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta)]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            fipBdTnM1 = Afix.'*xTilpnM1;
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1 = [mippdnM1; fippBdTnM1];
            %}
            %Limitação de variáveis cíclicas:
            if ind_LVC
                midnM1 = Amix.'*xTilnM1;
                xTilnM1 = xTilnM1 - 2*pi*[1;0]*double(midnM1>=2*pi);
            end
            %Trajetória desejada:
            y1dTilnM1 = y1dTil(tnM1);
            %Vetor de estados desejado:
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.'];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 4
        y2dTil0 = y2dTil(t0);
        fiBdU0 = KfiU\KlcUd(t0);
        fipBdU0 = 0;
        fippBdU0 = ffippBdU(t0,fipBdU0);
        Zd0 = fZud(t0);
        ZZ = zeros(nT,3*length(fiBdU0) + length(y2dTil0) + length(Zd0));
        ZZ(1,:) = [fiBdU0, fipBdU0, fippBdU0, y2dTil0.', Zd0.'];
        fiBdUn = fiBdU0;
        fipBdUn = fipBdU0;
        fippBdUn = fippBdU0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            fiBdUnM10 = fiBdUn + dtnM1*fipBdUn + 0.5*(dtnM1*fippBdUn)*dtnM1;
            Eta0 = fiBdUnM10;
            Psi = @(eta)(GamaFiBdUnM1(dtnM1,eta,fiBdUn,fipBdUn,fippBdUn) + KfiU*eta - KlcUd(tnM1));
            dPsidEtaT = @(eta)(dGamaFiBdUdfiBdUnM1(dtnM1,eta,fiBdUn,fipBdUn,fippBdUn) + KfiU);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            fiBdUnM1 = Eta_conv;
            fipBdUnM1 = fupnM1(dtnM1,fiBdUnM1,fiBdUn,fipBdUn,fippBdUn);
            %xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            fippBdUnM1 = ffippBdU(tnM1,fipBdUnM1);
            %}
            %Trajetória desejada:
            y2dTilnM1 = y2dTil(tnM1);
            %Vetor de estados desejado:
            ZdnM1 = fZud(tnM1);
            %Armazenamento:
            ZZ(i,:) = [fiBdUnM1, fipBdUnM1, fippBdUnM1, y2dTilnM1.', ZdnM1.'];
            %Atualização:
            fiBdUn = fiBdUnM1;
            fipBdUn = fipBdUnM1;
            fippBdUn = fippBdUnM1;
        end
    case 5
        y1dTil0 = y1dTil(t0);
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        xTil0 = [uTilT0; mid0; fiBdT0];
        xTilp0 = [uTilpT0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0);
        uTilppT0 = Mt\(Bt*Tx0 - Ct*uTilpT0 - Kt*uTilT0);
        xTilpp0 = [uTilppT0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        s10 = s1uTil(t0,uTilT0,uTilpT0);
        ZZ = zeros(nT,3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            %Integração:
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [KtNmkx(dtnM1)*eta - Bt*UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn) + Mtx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Ctx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMix*eta - SigmaMi*y1dTil(tnM1);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta)]);
            dPsidEtaT = @(eta)(...
                [KtNmkx(dtnM1) - Bt*DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn);
                 dmipddxTilnM1T(dtnM1) - GamaMix;
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta)]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([4,5],1) = [mippdnM1; fippBdTnM1];
            %}
            %Limitação de variáveis cíclicas:
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[ones(3,1);1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 6
        y2dTil0 = y2dTil(t0);
        uTilU0 = zeros(3,1);
        uTilpU0 = zeros(3,1);
        Z0 = [uTilU0; uTilpU0];
        fiBdU0 = KfiU\KlcUd(t0);
        fipBdU0 = 0;
        fippBdU0 = ffippBdU(t0,fipBdU0);
        xTil0 = [uTilU0; fiBdU0];
        xTilp0 = [uTilpU0; fipBdU0];
        Fu0 = UcUx(t0,xTil0,xTilp0);
        uTilppU0 = Mu\(Bu*Fu0 - Cu*uTilpU0 - Kbu*uTilU0);
        xTilpp0 = [uTilppU0; fippBdU0];
        Zd0 = fZudx(t0);
        s20 = s2uTil(t0,uTilU0,uTilpU0);
        ZZ = zeros(nT,3*length(xTil0) + length(y2dTil0) + length(Zd0) + length(Fu0) + length(s20));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y2dTil0.', Zd0.', Fu0, s20];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            %Integração:
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [KuNmkx(dtnM1)*eta - Bu*UcUxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn) + Mux*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cux*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 GamaFiBdUxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiUx*eta - KlcUdxnM1(tnM1)]);
            dPsidEtaT = @(eta)(...
                [KuNmkx(dtnM1) - Bu*DUcUxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn);
                 dGamaFiUxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiUx]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            fipBdUnM1 = AfiUx.'*xTilpnM1;
            fippBdUnM1 = ffippBdU(tnM1,fipBdUnM1);
            xTilppnM1(4,1) = fippBdUnM1;
            %}
            FunM1 = UcUxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            y2dTilnM1 = y2dTil(tnM1);
            ZdnM1 = fZudx(tnM1);
            s2nM1 = s2xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y2dTilnM1.', ZdnM1.', FunM1, s2nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 7
        y2dTil0 = y2dTilPS(t0);
        uTilU0 = zeros(3,1);
        uTilpU0 = zeros(3,1);
        Z0 = [uTilU0; uTilpU0];
        fiBdU0 = KfiUps\KlcUdPS(t0);
        fipBdU0 = 0;
        fippBdU0 = ffippBdUps(t0,fipBdU0);
        xTil0 = [uTilU0; fiBdU0];
        xTilp0 = [uTilpU0; fipBdU0];
        Fu0 = UcUxPS(t0,xTil0,xTilp0);
        uTilppU0 = Mu\(Bu*Fu0 - Cu*uTilpU0 - Kbu*uTilU0);
        xTilpp0 = [uTilppU0; fippBdU0];
        Zd0 = fZudxPS(t0);
        s20 = s2uTilPS(t0,uTilU0,uTilpU0);
        ZZ = zeros(nT,3*length(xTil0) + length(y2dTil0) + length(Zd0) + length(Fu0) + length(s20));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y2dTil0.', Zd0.', Fu0, s20];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            %Integração:
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [KuNmkx(dtnM1)*eta - Bu*UcUxnM1PS(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn) + Mux*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cux*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 GamaFiBdUxnM1PS(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiUxPS*eta - KlcUdxnM1PS(tnM1)]);
            dPsidEtaT = @(eta)(...
                [KuNmkx(dtnM1) - Bu*DUcUxDxnM1Tps(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn);
                 dGamaFiUxdxTilnM1Tps(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiUxPS]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            fipBdUnM1 = AfiUx.'*xTilpnM1;
            fippBdUnM1 = ffippBdUps(tnM1,fipBdUnM1);
            xTilppnM1(4,1) = fippBdUnM1;
            %}
            FunM1 = UcUxnM1PS(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            y2dTilnM1 = y2dTilPS(tnM1);
            ZdnM1 = fZudxPS(tnM1);
            s2nM1 = s2xnM1PS(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y2dTilnM1.', ZdnM1.', FunM1, s2nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 8
        y2dTil0 = y2dTil(t0);
        uTilU0 = zeros(3,1);
        uTilpU0 = zeros(3,1);
        Z0 = [uTilU0; uTilpU0];
        fiBdU0 = KfiU\KlcUd(t0);
        fipBdU0 = 0;
        fippBdU0 = ffippBdU(t0,fipBdU0);
        xTil0 = [uTilU0; fiBdU0];
        xTilp0 = [uTilpU0; fipBdU0];
        Fu0 = UcUx(t0,xTil0,xTilp0);
        uTilppU0 = Mu\(Bu*Fu0 - Cu*uTilpU0 - Kbu*uTilU0);
        xTilpp0 = [uTilppU0; fippBdU0];
        Zd0 = fZudx(t0);
        s20 = s2uTil(t0,uTilU0,uTilpU0);
        ZZ = zeros(nT,3*length(xTil0) + length(y2dTil0) + length(Zd0) + length(Fu0) + length(s20));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y2dTil0.', Zd0.', Fu0, s20];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            %Integração:
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [KuNmkx(dtnM1)*eta - Bu*UcUxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn) + Mux*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cux*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 GamaFiBdUxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiUx*eta - KlcUdxnM1(tnM1)]);
            dPsidEtaT = @(eta)(...
                [KuNmkx(dtnM1) - Bu*DUcUxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn);
                 dGamaFiUxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiUx]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            fipBdUnM1 = AfiUx.'*xTilpnM1;
            fippBdUnM1 = ffippBdU(tnM1,fipBdUnM1);
            xTilppnM1(4,1) = fippBdUnM1;
            %}
            FunM1 = UcUxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            y2dTilnM1 = y2dTil(tnM1);
            ZdnM1 = fZudx(tnM1);
            s2nM1 = s2xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y2dTilnM1.', ZdnM1.', FunM1, s2nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 9
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        fub1 = @(ubTil)(Id3(1,:)*ubTil);
        dub1bubTilT = Id3(1,:);
        Fbu0 = fFbu(tn);
        ubTilppU0 = Mu\(Bu*Fbu0 - Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        ZZ = zeros(nT,3*3+1);
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', Fbu0.'];
        ZZpd = zeros(nT-nPD,3*3+1);
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            %Integração:
            Eta0 = ubTilUnM10;
            Psi = @(eta)(KuNmk(dtnM1)*eta + fNuTil(eta) - Bu*fFbu(tnM1) +...
                Mu*fi2Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn) +...
                Cu*fi1Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsidEtaT = @(eta)(KuNmk(dtnM1));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            ubTilUnM1 = Eta_conv;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            FbunM1 = fFbu(tnM1);
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', FbunM1.'];
            if i>=nPD && MEANativo
                ZZpd(i-nPD+1,:) = mean(ZZ(i-nPD+1:i,:));
            end
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
        end
        %Integração por ode45:
        if ODEativo
            odefun = @(t,Z)(fAbu(Z)*Z + BbU*fFbu(t));
            tspan = T;
            Z0 = zeros(6,1);
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [Tode,ZZode] = ode45(odefun,tspan,Z0,options);
        end
    case 10
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        fub1 = @(ubTil)(Id3(1,:)*ubTil);
        dub1bubTilT = Id3(1,:);
        Fbu0 = fFbu(tn);
        ubTilppU0 = Mu\(Bu*Fbu0 - Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        
        y1dTil0 = y1dTil(t0);
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        xTil0 = [uTilT0; mid0; fiBdT0];
        xTilp0 = [uTilpT0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0);
        uTilppT0 = Mt\(Bt*Tx0 - Ct*uTilpT0 - Kt*uTilT0);
        xTilpp0 = [uTilppT0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        s10 = s1uTil(t0,uTilT0,uTilpT0);
        
        ZZ = zeros(nT,3*3+1 + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', Fbu0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            
            %Integração:
            EtaU0 = ubTilUnM10;
            PsiU = @(eta)(KuNmk(dtnM1)*eta + fNuTil(eta) - Bu*fFbu(tnM1) +...
                Mu*fi2Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn) +...
                Cu*fi1Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuNmk(dtnM1));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            ubTilUnM1 = EtaConvU;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            FbunM1 = fFbu(tnM1);
            
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [KtNmkx(dtnM1)*eta - Bt*UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn) + Mtx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Ctx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMix*eta - SigmaMi*y1dTil(tnM1);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta)]);
            dPsidEtaT = @(eta)(...
                [KtNmkx(dtnM1) - Bt*DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn);
                 dmipddxTilnM1T(dtnM1) - GamaMix;
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = EtaConv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([4,5],1) = [mippdnM1; fippBdTnM1];
            %}
            %Limitação de variáveis cíclicas:
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[ones(3,1);1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', FbunM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1];
            
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 11
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        fub1 = @(ubTil)(Id3(1,:)*ubTil);
        dub1bubTilT = Id3(1,:);
        Fbu0 = fFbu(tn);
        ubTilppU0 = Mu\(Bu*Fbu0 - Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        Xu0 = QsiDifUu(ubTilU0,ubTilpU0);
        y20 = y2(Xu0);
        yp20 = y2p(Xu0);
        y2p20 = y2pp(Xu0);
        y3p20 = y2ppp(Xu0);
        y4p20 = y2pppp(Xu0);
        
        y1dTil0 = y1dTil(t0);
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        xTil0 = [uTilT0; mid0; fiBdT0];
        xTilp0 = [uTilpT0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        uTilppT0 = Mt\(Bt*Tx0 - Ct*uTilpT0 - Kt*uTilT0);
        xTilpp0 = [uTilppT0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        
        ZZ = zeros(nT,3*3+1 + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', Fbu0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        %Zn = Z0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            %Integração:
            EtaU0 = ubTilUnM10;
            PsiU = @(eta)(KuNmk(dtnM1)*eta + fNuTil(eta) - Bu*fFbu(tnM1) +...
                Mu*fi2Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn) +...
                Cu*fi1Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuNmk(dtnM1) + fdNudubTilT(eta));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            ubTilUnM1 = EtaConvU;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            FbunM1 = fFbu(tnM1);
            
            XunM1 = QsiDifUu(ubTilUnM1,ubTilpUnM1);
            y2nM1 = y2(XunM1);
            yp2nM1 = y2p(XunM1);
            y2p2nM1 = y2pp(XunM1);
            y3p2nM1 = y2ppp(XunM1);
            y4p2nM1 = y2pppp(XunM1);
            
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [KtNmkx(dtnM1)*eta + FatTcxnM1(dtnM1,eta,y2nM1,xTiln,xTilpn,xTilppn) - Bt*UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + Mtx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Ctx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtxnM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [KtNmkx(dtnM1) + DFatTcDxTilTnM1(dtnM1,eta,y2nM1,xTiln,xTilpn,xTilppn) - Bt*DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                 dmipddxTilnM1T(dtnM1) - fdmipdAtdxnM1T(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = EtaConv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipdAt(tnM1,midnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            mippdnM1 = fmippdAt(tnM1,midnM1,mipdnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([4,5],1) = [mippdnM1; fippBdTnM1];
            %}
            %Limitação de variáveis cíclicas:
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[ones(3,1);1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', FbunM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1];
            
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 12
        y1dTil0 = y1dTil(t0);
        uTilT0 = zeros(3,1);
        upTilT0 = zeros(3,1);
        Z0 = [uTilT0; upTilT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        ubTilU0 = zeros(3,1);
        ubTilAc0 = Atu*uTilT0 + Auu*ubTilU0;
        ubpTilU0 = zeros(3,1);
        ubpTilAc0 = Atu*upTilT0 + Auu*ubpTilU0;
        xTil0 = [ubTilAc0; mid0; fiBdT0];
        xTilp0 = [ubpTilAc0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0);
        uppTilT0 = Mt\(Bt*Tx0 - Ct*upTilT0 - Kt*uTilT0);
        Fbu0 = fFbu(tn);
        ubppTilU0 = Mu\(Bu*Fbu0 - Cu*ubpTilU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        ubppTilAc0 = Atu*uppTilT0 + Auu*ubppTilU0;
        xTilpp0 = [ubppTilAc0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        s10 = s1uTil(t0,uTilT0,upTilT0);
        ZZ = zeros(nT,3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(Fbu0) + length(s10));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, Fbu0, s10];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            %Integração:
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [Knmkx(dtnM1)*eta + fTilcxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bac*[UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn); fFbu(tnM1)] + M1acx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cacx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMix*eta - SigmaMi*y1dTil(tnM1);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta)]);
            dPsidEtaT = @(eta)(...
                [Knmkx(dtnM1) + DfTilcDxnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bac*[DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn); zeros(1,10)];
                 dmipddxTilnM1T(dtnM1) - GamaMix;
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta)]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([9,10],1) = [mippdnM1; fippBdTnM1];
            %}
            %Limitação de variáveis cíclicas:
            if ind_LVC
                ubTilAcnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id6(1,:)*ubTilAcnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            FbunM1 = fFbu(tnM1);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, FbunM1, s1nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 13
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        y1dTil0 = y1dTil(t0);
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        %fub1 = @(ubTil)(Id3(1,:)*ubTil);
        %dub1bubTilT = Id3(1,:);
        ubTilAc0 = Atu*uTilT0 + Auu*ubTilU0;
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        Fbu0 = fFbuUcuTil(t0,ubTilAc0);
        
        ubTilppU0 = Mu\(Bu*Fbu0 - Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        Xu0 = QsiDifUu(ubTilU0,ubTilpU0);
        y20 = y2(Xu0);
        yp20 = y2p(Xu0);
        y2p20 = y2pp(Xu0);
        y3p20 = y2ppp(Xu0);
        y4p20 = y2pppp(Xu0);
        
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        ubTilU0 = zeros(3,1);
        ubTilAc0 = Atu*uTilT0 + Auu*ubTilU0;
        ubpTilU0 = zeros(3,1);
        ubpTilAc0 = Atu*uTilpT0 + Auu*ubpTilU0;
        
        xTil0 = [ubTilAc0; mid0; fiBdT0];
        xTilp0 = [ubpTilAc0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        uppTilT0 = Mt\(Bt*Tx0 - Ct*uTilpT0 - Kt*uTilT0);
        %uppTilT0 = zeros(3,1);
        ubppTilAc0 = Atu*uppTilT0 + Auu*ubTilppU0;
        ubppTilAc0 = M1ac\(Bac*[Tx0; Fbu0] - Cac*ubpTilAc0 - Kac*ubTilAc0 -...
            fTilc(ubTilAc0,ubpTilAc0,ubppTilAc0));
        xTilpp0 = [ubppTilAc0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        
        ZZ = zeros(nT,3*length(ubTilU0) + length(Fbu0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, Fbu0, s10];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        %Zn = Z0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            %Integração:
            epsilon = argumentos.epsNewtRaph;
            argumentos.epsNewtRaph = epsilon*10^-12;
            EtaU0 = ubTilUnM10;
            PsiU = @(eta)(KuNmk(dtnM1)*eta + fNuTil(eta) - Bu*fFbuUcuTilU(tnM1,eta) +...
                Mu*fi2Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn) +...
                Cu*fi1Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuNmk(dtnM1));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            ubTilUnM1 = EtaConvU;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            argumentos.epsNewtRaph = epsilon;
            %FbunM1 = fFbu(tnM1);
            
            XunM1 = QsiDifUu(ubTilUnM1,ubTilpUnM1);
            y2nM1 = y2(XunM1);
            yp2nM1 = y2p(XunM1);
            y2p2nM1 = y2pp(XunM1);
            y3p2nM1 = y2ppp(XunM1);
            y4p2nM1 = y2pppp(XunM1);
            
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [Knmkx(dtnM1)*eta + fTilcxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bac*[UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10); fFbuUcxnM1(tnM1,eta)] + M1acx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cacx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtxnM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [Knmkx(dtnM1) + DfTilcDxnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bac*[DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10); dFbuUcdxT];
                 dmipddxTilnM1T(dtnM1) - fdmipdAtdxnM1T(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([9,10],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            FbunM1 = fFbuUcxnM1(tnM1,xTilnM1);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, FbunM1, s1nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 14
        y1dTil0 = y1dTil(t0);
        uTilT0 = zeros(3,1);
        upTilT0 = zeros(3,1);
        Z0 = [uTilT0; upTilT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        ubTilU0 = zeros(3,1);
        ubTilAc0 = Atu*uTilT0 + Auu*ubTilU0;
        ubpTilU0 = zeros(3,1);
        ubpTilAc0 = Atu*upTilT0 + Auu*ubpTilU0;
        xTil0 = [ubTilAc0; mid0; fiBdT0];
        xTilp0 = [ubpTilAc0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0);
        uppTilT0 = Mt\(Bt*Tx0 - Ct*upTilT0 - Kt*uTilT0);
        Fbu0 = fFbuUcuTil(t0,ubTilAc0);
        ubppTilU0 = Mu\(Bu*Fbu0 - Cu*ubpTilU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        ubppTilAc0 = Atu*uppTilT0 + Auu*ubppTilU0;
        xTilpp0 = [ubppTilAc0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        s10 = s1uTil(t0,uTilT0,upTilT0);
        ZZ = zeros(nT,3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(Fbu0) + length(s10));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, Fbu0, s10];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            %Integração:
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [Knmkx(dtnM1)*eta + fTilcxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bac*[UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn); fFbuUcxnM1(tnM1,eta)] + M1acx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cacx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMix*eta - SigmaMi*y1dTil(tnM1);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta)]);
            dPsidEtaT = @(eta)(...
                [Knmkx(dtnM1) + DfTilcDxnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bac*[DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn); dFbuUcdxT];
                 dmipddxTilnM1T(dtnM1) - GamaMix;
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta)]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([9,10],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            FbunM1 = fFbuUcx(tnM1,xTilnM1);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, FbunM1, s1nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 15
        Tx0 = fTxUc(t0);
        Fbu0 = fFbuUc(t0);
        ubTilAc0 = zeros(8,1);
        ubpTilAc0 = zeros(8,1);
        ubppTilAc0 = zeros(8,1);
        ubppTilAc0 = M1ac\(Bac*[Tx0; Fbu0] - Cac*ubpTilAc0 - Kac*ubTilAc0 - fTilSmAt(ubTilAc0,ubpTilAc0,ubppTilAc0));
        ZZ = zeros(nT,3*length(ubTilAc0) + length(Tx0) + length(Fbu0));
        ZZ(1,:) = [ubTilAc0.', ubpTilAc0.', ubppTilAc0.', Tx0, Fbu0];
        ubTilAcn = ubTilAc0;
        ubpTilAcn = ubpTilAc0;
        ubppTilAcn = ubppTilAc0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilAcnM10 = ubTilAcn + dtnM1*ubpTilAcn + 0.5*(dtnM1*ubppTilAcn)*dtnM1;
            %Integração:
            Eta0 = ubTilAcnM10;
            Psi = @(eta)(...
                Knmk(dtnM1)*eta + fTilnM1(dtnM1,eta,ubTilAcn,ubpTilAcn,ubppTilAcn) - Bac*[fTxUc(tnM1); fFbuUc(tnM1)] +...
                M1ac*fi2Nmk(dtnM1,ubTilAcn,ubpTilAcn,ubppTilAcn) +...
                Cac*fi1Nmk(dtnM1,ubTilAcn,ubpTilAcn,ubppTilAcn));
            dPsidEtaT = @(eta)(...
                Knmk(dtnM1) + DfTilDubnM1T(dtnM1,eta,ubTilAcn,ubpTilAcn,ubppTilAcn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            ubTilAcnM1 = Eta_conv;
            ubpTilAcnM1 = fupnM1(dtnM1,ubTilAcnM1,ubTilAcn,ubpTilAcn,ubppTilAcn);
            ubppTilAcnM1 = fuppnM1(dtnM1,ubTilAcnM1,ubTilAcn,ubpTilAcn,ubppTilAcn);
            TxnM1 = fTxUc(tnM1);
            FbunM1 = fFbuUc(tnM1);
            %Armazenamento:
            ZZ(i,:) = [ubTilAcnM1.', ubpTilAcnM1.', ubppTilAcnM1.', TxnM1, FbunM1];
            %Atualização:
            ubTilAcn = ubTilAcnM1;
            ubpTilAcn = ubpTilAcnM1;
            ubppTilAcn = ubppTilAcnM1;
        end
    case 16
        Tx0 = fTxUc(t0);
        ubTilAc0 = zeros(8,1);
        ubpTilAc0 = zeros(8,1);
        ubppTilAc0 = zeros(8,1);
        Fbu0 = fFbuUcuTil(t0,ubTilAc0);
        ubppTilAc0 = M1ac\(Bac*[Tx0; Fbu0] - Cac*ubpTilAc0 - Kac*ubTilAc0 - fTilSmAt(ubTilAc0,ubpTilAc0,ubppTilAc0));
        ZZ = zeros(nT,3*length(ubTilAc0) + length(Tx0) + length(Fbu0));
        ZZ(1,:) = [ubTilAc0.', ubpTilAc0.', ubppTilAc0.', Tx0, Fbu0];
        ubTilAcn = ubTilAc0;
        ubpTilAcn = ubpTilAc0;
        ubppTilAcn = ubppTilAc0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilAcnM10 = ubTilAcn + dtnM1*ubpTilAcn + 0.5*(dtnM1*ubppTilAcn)*dtnM1;
            %Integração:
            Eta0 = ubTilAcnM10;
            Psi = @(eta)(...
                Knmk(dtnM1)*eta + fTilnM1(dtnM1,eta,ubTilAcn,ubpTilAcn,ubppTilAcn) - Bac*[fTxUc(tnM1); fFbuUcuTil(tnM1,eta)] +...
                M1ac*fi2Nmk(dtnM1,ubTilAcn,ubpTilAcn,ubppTilAcn) +...
                Cac*fi1Nmk(dtnM1,ubTilAcn,ubpTilAcn,ubppTilAcn));
            dPsidEtaT = @(eta)(...
                Knmk(dtnM1) + DfTilDubnM1T(dtnM1,eta,ubTilAcn,ubpTilAcn,ubppTilAcn) - Bac*[zeros(1,8); dFbuUcdubT]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            ubTilAcnM1 = Eta_conv;
            ubpTilAcnM1 = fupnM1(dtnM1,ubTilAcnM1,ubTilAcn,ubpTilAcn,ubppTilAcn);
            ubppTilAcnM1 = fuppnM1(dtnM1,ubTilAcnM1,ubTilAcn,ubpTilAcn,ubppTilAcn);
            TxnM1 = fTxUc(tnM1);
            FbunM1 = fFbuUcuTil(tnM1,ubTilAcnM1);
            if ind_LVC
                teta1nM1 = Id6(1,:)*ubTilAcnM1;
                ubTilAcnM1 = ubTilAcnM1 - 2*pi*[0;1;0;0;0;1;0;1]*double(teta1nM1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [ubTilAcnM1.', ubpTilAcnM1.', ubppTilAcnM1.', TxnM1, FbunM1];
            %Atualização:
            ubTilAcn = ubTilAcnM1;
            ubpTilAcn = ubpTilAcnM1;
            ubppTilAcn = ubppTilAcnM1;
        end
    case 17
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        fub1 = @(ubTil)(Id3(1,:)*ubTil);
        dub1bubTilT = Id3(1,:);
        Fbu0 = fFbuUcuTilU(t0,ubTilU0);
        ubTilppU0 = Mu\(Bu*Fbu0 - Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        
        y1dTil0 = y1dTil(t0);
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        xTil0 = [uTilT0; mid0; fiBdT0];
        xTilp0 = [uTilpT0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0);
        uTilppT0 = Mt\(Bt*Tx0 - Ct*uTilpT0 - Kt*uTilT0);
        xTilpp0 = [uTilppT0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        s10 = s1uTil(t0,uTilT0,uTilpT0);
        
        ZZ = zeros(nT,3*3+1 + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', Fbu0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            
            %Integração:fFbuUcuTilU(t0,ubTilU0)
            EtaU0 = ubTilUnM10;
            PsiU = @(eta)(KuNmk(dtnM1)*eta + fNuTil(eta) - Bu*fFbuUcuTilU(tnM1,eta) +...
                Mu*fi2Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn) +...
                Cu*fi1Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuNmk(dtnM1) + fdNudubTilT(eta) - Bu*dFbuUcdubT);
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            ubTilUnM1 = EtaConvU;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            FbunM1 = fFbuUcuTilU(tnM1,ubTilUnM1);
            
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [KtNmkx(dtnM1)*eta - Bt*UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn) + Mtx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Ctx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMix*eta - SigmaMi*y1dTil(tnM1);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta)]);
            dPsidEtaT = @(eta)(...
                [KtNmkx(dtnM1) - Bt*DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn);
                 dmipddxTilnM1T(dtnM1) - GamaMix;
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = EtaConv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([4,5],1) = [mippdnM1; fippBdTnM1];
            %}
            %Limitação de variáveis cíclicas:
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[ones(3,1);1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', FbunM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1];
            
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 18
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        fub1 = @(ubTil)(Id3(1,:)*ubTil);
        dub1bubTilT = Id3(1,:);
        Fbu0 = fFbu(tn,fub1(ubTilU0));
        dFbudub1 = -Kf;
        dFbudubTilT = dFbudub1*dub1bubTilT;
        ubTilppU0 = Mu\(Bu*Fbu0 - Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        ub30 = Id3(1,:)*ubTilU0;
        y20 = fNub3(ub30);
        Tx0 = fTx(tn);
        uTilppT0 = Mt\(Bt*Tx0 - Ct*uTilpT0 - Kt*uTilT0 - FatTc(uTilpT0,y20));
        ZZ = zeros(nT,3*3+1+3*3+1);
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', Fbu0.',uTilT0.', uTilpT0.', uTilppT0.', Tx0.'];
        uTilTn = uTilT0;
        uTilpTn = uTilpT0;
        uTilppTn = uTilppT0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            uTilTnM10 = uTilTn + dtnM1*uTilpTn + 0.5*(dtnM1*uTilppTn)*dtnM1;
            %Integração:
            Eta0 = ubTilUnM10;
            Psi = @(eta)(KuNmk(dtnM1)*eta + fNuTil(eta) - Bu*fFbu(tnM1,fub1(eta)) +...
                Mu*fi2Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn) +...
                Cu*fi1Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsidEtaT = @(eta)(KuNmk(dtnM1) + fdNudubTilT(eta) - Bu*dFbudubTilT);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            ubTilUnM1 = Eta_conv;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ub1nM1 = fub1(ubTilUnM1);
            FbunM1 = fFbu(tnM1,ub1nM1);
            
            ub3nM1 = fub3(ubTilUnM1);
            y2nM1 = fNub3(ub3nM1);
            
            Eta0 = uTilTnM10;
            %Integração com modelo de atrito contínuo na origem:
            %{
            Psi = @(eta)(KtNmk(dtnM1)*eta + FatTcnM1(dtnM1,eta,uTilTn,uTilpTn,uTilppTn,y2nM1) - Bt*fTx(tnM1) +...
                Mt*fi2Nmk(dtnM1,uTilTn,uTilpTn,uTilppTn) +...
                Ct*fi1Nmk(dtnM1,uTilTn,uTilpTn,uTilppTn));
            dPsidEtaT = @(eta)(KtNmk(dtnM1) + dFatTcduTilTnM1(dtnM1,eta,uTilTn,uTilpTn,uTilppTn,y2nM1));
            %}
            %Integração com modelo de atrito descontínuo na origem:
            %
            Psi = @(eta)(KtNmk(dtnM1)*eta + FatTdcnM1(dtnM1,eta,uTilTn,uTilpTn,uTilppTn,y2nM1) - Bt*fTx(tnM1) +...
                Mt*fi2Nmk(dtnM1,uTilTn,uTilpTn,uTilppTn) +...
                Ct*fi1Nmk(dtnM1,uTilTn,uTilpTn,uTilppTn));
            dPsidEtaT = @(eta)(KtNmk(dtnM1) + dFatTdcduTilTnM1(dtnM1,eta,uTilTn,uTilpTn,uTilppTn,y2nM1));
            %}
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            uTilTnM1 = Eta_conv;
            uTilpTnM1 = fupnM1(dtnM1,uTilTnM1,uTilTn,uTilpTn,uTilppTn);
            uTilppTnM1 = fuppnM1(dtnM1,uTilTnM1,uTilTn,uTilpTn,uTilppTn);
            TxnM1 = fTx(tnM1);
            %Limitação de variáveis cíclicas:
            if ind_LVC
                %Limitação das variáveis cíclicas:
                teta1 = Id3(1,:)*uTilTnM1;
                uTilTnM1 = uTilTnM1 - 2*pi*ones(3,1)*double(teta1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', FbunM1.',uTilTnM1.', uTilpTnM1.', uTilppTnM1.', TxnM1.'];
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            uTilTn = uTilTnM1;
            uTilpTn = uTilpTnM1;
            uTilppTn = uTilppTnM1;
        end
    case 19
        ubTilr0 = zeros(nIdNgr,1);
        ubTilpr0 = zeros(nIdNgr,1);
        ubTilppr0 = zeros(nIdNgr,1);
        Tx0 = fTx(tn);
        Fbu0 = fFbuTilr(ubTilr0);
        Uc0 = [Tx0; Fbu0];
        ubTilppr0 = M1restr\(Br*Uc0 - Crestr*ubTilpr0 - Krestr*ubTilr0 -...
            fcEFbr(ubTilr0,ubTilpr0,ubTilppr0));
        ZZ = zeros(nT,3*nIdNgr+1+1);
        ZZ(1,:) = [ubTilr0.', ubTilpr0.', ubTilppr0.', Tx0, Fbu0];
        ubTilrn = ubTilr0;
        ubTilprn = ubTilpr0;
        ubTilpprn = ubTilppr0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilrnM10 = ubTilrn + dtnM1*ubTilprn + 0.5*(dtnM1*ubTilpprn)*dtnM1;
            %Integração:
            Eta0 = ubTilrnM10;
            Psi = @(eta)(Knmk(dtnM1)*eta + fcEFbrnM1(dtnM1,eta,ubTilrn,ubTilprn,ubTilpprn) -...
                Br*[fTx(tnM1); fFbuTilr(eta)] +...
                M1restr*fi2Nmk(dtnM1,ubTilrn,ubTilprn,ubTilpprn) +...
                Crestr*fi1Nmk(dtnM1,ubTilrn,ubTilprn,ubTilpprn));
            dPsidEtaT = @(eta)(Knmk(dtnM1) +...
                DfcEFbrDubrnM1T(dtnM1,eta,ubTilrn,ubTilprn,ubTilpprn) -...
                Br*[dTxduTilrT; dFburduTilrT]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            ubTilrnM1 = Eta_conv;
            ubTilprnM1 = fupnM1(dtnM1,ubTilrnM1,ubTilrn,ubTilprn,ubTilpprn);
            ubTilpprnM1 = fuppnM1(dtnM1,ubTilrnM1,ubTilrn,ubTilprn,ubTilpprn);
            TxnM1 = fTx(tnM1);
            FbunM1 = fFbuTilr(ubTilrnM1);
            %Limitação de variáveis cíclicas:
            if ind_LVC
                %Limitação das variáveis cíclicas:
                tetaX1 = IdNgr(6*1-4,:)*ubTilrnM1;
                ubTilrnM1 = ubTilrnM1 - 2*pi*IdLVCr*double(tetaX1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [ubTilrnM1.', ubTilprnM1.', ubTilpprnM1.', TxnM1, FbunM1];
            %Atualização:
            ubTilrn = ubTilrnM1;
            ubTilprn = ubTilprnM1;
            ubTilpprn = ubTilpprnM1;
        end
    case 20
        ubTilr0 = zeros(nIdNgr,1);
        ubTilpr0 = zeros(nIdNgr,1);
        ubTilppr0 = zeros(nIdNgr,1);
        
        y1dTil0 = y1dTil(t0);
        uTilT0 = IdTetaX([1, N1+1, N+1],:)*(SigmaReF*ubTilr0);
        upTilT0 = IdTetaX([1, N1+1, N+1],:)*(SigmaReF*ubTilpr0);
        Z0 = [uTilT0; upTilT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTil0 = [ubTilr0; mid0; fiBdT0];
        xTilp0 = [ubTilpr0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0);
        Fbu0 = fFbuUcuTil(t0,ubTilr0);
        Uc0 = [Tx0; Fbu0];
        ubTilppr0 = M1restr\(Br*Uc0 - Crestr*ubTilpr0 - Krestr*ubTilr0 -...
            fcEFbr(ubTilr0,ubTilpr0,ubTilppr0));
        xTilpp0 = [ubTilppr0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        s10 = s1uTil(t0,uTilT0,upTilT0);
        ZZ = zeros(nT,3*length(xTil0) + length(y1dTil0) + length(Zd0) +...
            length(Tx0) + length(Fbu0) + length(s10));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, Fbu0, s10];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            %Integração:
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [Knmkx(dtnM1)*eta + fTilcxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - Br*[UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn); fFbuUcxnM1(tnM1,eta)] + M1restrx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Crestrx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMix*eta - SigmaMi*y1dTil(tnM1);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta)]);
            dPsidEtaT = @(eta)(...
                [Knmkx(dtnM1) + DfTilcDxnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) - Br*[DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn); dFbuUcdxT];
                 dmipddxTilnM1T(dtnM1) - GamaMix;
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta)]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([nIdNgr+1,nIdNgr+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                ubTilrnM1 = Aurx.'*xTilnM1;
                tetaX1 = IdTetaX(1,:)*(SigmaReF*ubTilrnM1);
                xTilnM1 = xTilnM1 - 2*pi*[IdLVCr; 1; 0]*double(tetaX1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            FbunM1 = fFbuUcx(tnM1,xTilnM1);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, FbunM1, s1nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 21
        ubTil0 = zeros(Ng,1);
        ubTilp0 = zeros(Ng,1);
        ubTilpp0 = zeros(Ng,1);
        Tx0 = fTx(t0);
        Fbu0 = fFbuTil(ubTil0);
        Uc0 = [Tx0; Fbu0];
        ubTilpp0 = M1gEF\(Bef*Uc0 - CgEF*ubTilp0 - KgEF*ubTil0 -...
            fcEFb(ubTil0,ubTilp0,ubTilpp0));
        lmb0 = zeros(rEstr,1);
        ZZ = zeros(nT,3*Ng+rEstr+1+1);
        ZZ(1,:) = [ubTil0.', ubTilp0.', ubTilpp0.', lmb0.', Tx0, Fbu0];
        ubTiln = ubTil0;
        ubTilpn = ubTilp0;
        ubTilppn = ubTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilnM10 = ubTiln + dtnM1*ubTilpn + 0.5*(dtnM1*ubTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [ubTilnM10; lmbnM10];
            %Integração:
            Eta0 = znM10;
            Psi = @(eta)(Knmkz(dtnM1)*eta + fcEFbznM1(dtnM1,eta,ubTiln,ubTilpn,ubTilppn) -...
                Bz*fUcz(tnM1,eta) + fiNmkzn(dtnM1,ubTiln,ubTilpn,ubTilppn));
            dPsidEtaT = @(eta)(Knmkz(dtnM1) + dfcEFzdzT(dtnM1,eta,ubTiln,ubTilpn,ubTilppn) -...
                Bz*fdUczdzT);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            ubTilnM1 = Auz.'*znM1;
            lmbnM1 = Almbz.'*znM1;
            ubTilpnM1 = fupnM1(dtnM1,ubTilnM1,ubTiln,ubTilpn,ubTilppn);
            ubTilppnM1 = fuppnM1(dtnM1,ubTilnM1,ubTiln,ubTilpn,ubTilppn);
            TxnM1 = fTx(tnM1);
            FbunM1 = fFbuTil(ubTilnM1);
            %Limitação de variáveis cíclicas:
            if ind_LVC
                %Limitação das variáveis cíclicas:
                tetaX1 = IdNg(6*1,:)*ubTilnM1;
                ubTilnM1 = ubTilnM1 - 2*pi*IdLVC*double(tetaX1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [ubTilnM1.', ubTilpnM1.', ubTilppnM1.', lmbnM1.', TxnM1, FbunM1];
            %Atualização:
            ubTiln = ubTilnM1;
            ubTilpn = ubTilpnM1;
            ubTilppn = ubTilppnM1;
        end
    case 22
        ubTil0 = zeros(Ng,1);
        ubTilp0 = zeros(Ng,1);
        ubTilpp0 = zeros(Ng,1);
        lmb0 = zeros(rEstr,1);
        
        y1dTil0 = y1dTil(t0);
        uTilT0 = IdTetaX([1, N1+1, N+1],:)*(ubTil0);
        upTilT0 = IdTetaX([1, N1+1, N+1],:)*(ubTilp0);
        Z0 = [uTilT0; upTilT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTil0 = [ubTil0; mid0; fiBdT0];
        xTilp0 = [ubTilp0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0);
        Fbu0 = fFbuUcuTil(t0,ubTil0);
        Uc0 = [Tx0; Fbu0];
        ubTilpp0 = M1gEF\(Bef*Uc0 - CgEF*ubTilp0 - KgEF*ubTil0 -...
            fcEFb(ubTil0,ubTilp0,ubTilpp0));
        xTilpp0 = [ubTilpp0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        s10 = s1uTil(t0,uTilT0,upTilT0);
        ZZ = zeros(nT,3*length(xTil0) + length(lmb0) + length(y1dTil0) + length(Zd0) +...
            length(Tx0) + length(Fbu0) + length(s10));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', lmb0.', y1dTil0.', Zd0.', Tx0, Fbu0, s10];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [xTilnM10; lmbnM10];
            %Integração:
            Eta0 = znM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + kNmk*GamaRz.'*eta - Bef*[UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn); fFbuUcxznM1(tnM1,eta)] + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMixz*eta - SigmaMi*y1dTil(tnM1);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta);
                 kNmk*GamaRxz*eta]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + kNmk*GamaRz.' - Bef*[dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn); dFbuUczdzT];
                 dmipddzT(dtnM1) - GamaMixz;
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta);
                 kNmk*GamaRxz]);
            %{
            ad = dPsidEtaT(znM10);
            %}
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            xTilnM1 = Axz.'*znM1;
            lmbnM1 = Almbz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                ubTilnM1 = Aux.'*xTilnM1;
                tetaX1 = IdTetaX(1,:)*(ubTilnM1);
                xTilnM1 = xTilnM1 - 2*pi*[IdLVC; 1; 0]*double(tetaX1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            FbunM1 = fFbuUcx(tnM1,xTilnM1);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', lmbnM1.', y1dTilnM1.', ZdnM1.', TxnM1, FbunM1, s1nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
        end
    case 23
        ubTil0 = zeros(Ng,1);
        ubTilp0 = zeros(Ng,1);
        ubTilpp0 = zeros(Ng,1);
        lmb0 = zeros(rEstrT,1);
        
        y1dTil0 = y1dTil(t0);
        uTilT0 = IdTetaX([1, N1+1, N+1],:)*(ubTil0);
        upTilT0 = IdTetaX([1, N1+1, N+1],:)*(ubTilp0);
        Z0 = [uTilT0; upTilT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTil0 = [ubTil0; mid0; fiBdT0];
        xTilp0 = [ubTilp0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0);
        ubTilpp0 = M1gEF\(BefTx*Tx0 - CgEF*ubTilp0 - KgEF*ubTil0 -...
            fcEFb(ubTil0,ubTilp0,ubTilpp0));
        xTilpp0 = [ubTilpp0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0);
        s10 = s1uTil(t0,uTilT0,upTilT0);
        ZZ = zeros(nT,3*length(xTil0) + length(lmb0) + length(y1dTil0) + length(Zd0) +...
            length(Tx0) + length(s10));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', lmb0.', y1dTil0.', Zd0.', Tx0, s10];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [xTilnM10; lmbnM10];
            %Integração:
            Eta0 = znM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + kNmk*GamaRz.'*eta - BefTx*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMixz*eta - SigmaMi*y1dTil(tnM1);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta);
                 kNmk*(GamaRxz*eta - fhRestr(tnM1))]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + kNmk*GamaRz.' - BefTx*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn);
                 dmipddzT(dtnM1) - GamaMixz;
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta);
                 kNmk*GamaRxz]);
            %{
            ad = dPsidEtaT(znM10);
            %}
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            xTilnM1 = Axz.'*znM1;
            lmbnM1 = Almbz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                ubTilnM1 = Aux.'*xTilnM1;
                tetaX1 = IdTetaX(1,:)*(ubTilnM1);
                xTilnM1 = xTilnM1 - 2*pi*[IdLVC; 1; 0]*double(tetaX1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', lmbnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
        end
    case 24
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        fub1 = @(ubTil)(Id3(1,:)*ubTil);
        dub1bubTilT = Id3(1,:);
        Fbu0 = fFbuUcuTilU(t0,ubTilU0);
        ubTilppU0 = Mu\(Bu*Fbu0 - Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        Xu0 = QsiDifUu(ubTilU0,ubTilpU0);
        y20 = y2(Xu0);
        yp20 = y2p(Xu0);
        y2p20 = y2pp(Xu0);
        y3p20 = y2ppp(Xu0);
        y4p20 = y2pppp(Xu0);
        
        y1dTil0 = y1dTil(t0);
        uTilT0 = zeros(3,1);
        uTilpT0 = zeros(3,1);
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        xTil0 = [uTilT0; mid0; fiBdT0];
        xTilp0 = [uTilpT0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        uTilppT0 = Mt\(Bt*Tx0 - Ct*uTilpT0 - Kt*uTilT0);
        xTilpp0 = [uTilppT0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        
        ZZ = zeros(nT,3*3+1 + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', Fbu0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        %Zn = Z0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            %Integração:
            EtaU0 = ubTilUnM10;
            PsiU = @(eta)(KuNmk(dtnM1)*eta + fNuTil(eta) - Bu*fFbuUcuTilU(tnM1,eta) +...
                Mu*fi2Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn) +...
                Cu*fi1Nmk(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuNmk(dtnM1) + fdNudubTilT(eta) - Bu*dFbuUcdubT);
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            ubTilUnM1 = EtaConvU;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            FbunM1 = fFbuUcuTilU(tnM1,ubTilUnM1);
            
            XunM1 = QsiDifUu(ubTilUnM1,ubTilpUnM1);
            y2nM1 = y2(XunM1);
            yp2nM1 = y2p(XunM1);
            y2p2nM1 = y2pp(XunM1);
            y3p2nM1 = y2ppp(XunM1);
            y4p2nM1 = y2pppp(XunM1);
            
            Eta0 = xTilnM10;
            Psi = @(eta)(...
                [KtNmkx(dtnM1)*eta + FatTcxnM1(dtnM1,eta,y2nM1,xTiln,xTilpn,xTilppn) - Bt*UcTxnM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + Mtx*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Ctx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                 fmidpxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtxnM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 GamaFiBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx*eta - KlcTdxnM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [KtNmkx(dtnM1) + DFatTcDxTilTnM1(dtnM1,eta,y2nM1,xTiln,xTilpn,xTilppn) - Bt*DUcTxDxnM1T(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                 dmipddxTilnM1T(dtnM1) - fdmipdAtdxnM1T(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 dGamaFiTxdxTilnM1T(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTx - dKlcTdxdxTilnM1T(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = EtaConv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipdAt(tnM1,midnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            mippdnM1 = fmippdAt(tnM1,midnM1,mipdnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([4,5],1) = [mippdnM1; fippBdTnM1];
            %}
            %Limitação de variáveis cíclicas:
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[ones(3,1);1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', FbunM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1];
            
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 25
        ubTil0 = zeros(Ng,1);
        ubTilp0 = zeros(Ng,1);
        ubTilpp0 = zeros(Ng,1);
        lmb0 = zeros(rEstr+2,1);
        ZZ = zeros(nT,3*Ng+rEstr+2);
        ZZ(1,:) = [ubTil0.', ubTilp0.', ubTilpp0.', lmb0.'];
        ubTiln = ubTil0;
        ubTilpn = ubTilp0;
        ubTilppn = ubTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilnM10 = ubTiln + dtnM1*ubTilpn + 0.5*(dtnM1*ubTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [ubTilnM10; lmbnM10];
            %Integração:
            Eta0 = znM10;
            Psi = @(eta)(Knmkz(dtnM1)*eta + fcEFbznM1(dtnM1,tnM1,eta,ubTiln,ubTilpn,ubTilppn) +...
                fiNmkzn(dtnM1,ubTiln,ubTilpn,ubTilppn));
            dPsidEtaT = @(eta)(Knmkz(dtnM1) + dfcEFzdzT(dtnM1,eta,ubTiln,ubTilpn,ubTilppn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            ubTilnM1 = Auz.'*znM1;
            lmbnM1 = Almbz.'*znM1;
            ubTilpnM1 = fupnM1(dtnM1,ubTilnM1,ubTiln,ubTilpn,ubTilppn);
            ubTilppnM1 = fuppnM1(dtnM1,ubTilnM1,ubTiln,ubTilpn,ubTilppn);
            %Limitação de variáveis cíclicas:
            if ind_LVC
                %Limitação das variáveis cíclicas:
                tetaX1 = IdNg(6*1,:)*ubTilnM1;
                ubTilnM1 = ubTilnM1 - 2*pi*IdLVC*double(tetaX1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [ubTilnM1.', ubTilpnM1.', ubTilppnM1.', lmbnM1.'];
            %Atualização:
            ubTiln = ubTilnM1;
            ubTilpn = ubTilpnM1;
            ubTilppn = ubTilppnM1;
        end
    case 26
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        ubTilppU0 = Mu\(-Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        lmbl0 = 0;
        Xu0 = QsiDifUu(ubTilU0,ubTilpU0);
        y20 = y2(Xu0);
        yp20 = y2p(Xu0);
        y2p20 = y2pp(Xu0);
        y3p20 = y2ppp(Xu0);
        y4p20 = y2pppp(Xu0);
        
        ubEF0 = zeros(Ng,1);
        ubpEF0 = zeros(Ng,1);
        uTilT0 = Atu.'*ubEF0;
        uTilpT0 = Atu.'*ubpEF0;
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        xTil0 = [ubEF0; mid0; fiBdT0];
        mipd0 = fmipd(t0,mid0);
        xTilp0 = [ubpEF0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        ubppEF0 = zeros(Ng,1);
        ubppEF0 = M1gEF\(BefTx*Tx0 - CgEF*ubpEF0 - KgEF*ubEF0 -...
            fcEFb(ubEF0,ubpEF0,ubppEF0));
        
        mippd0 = fmippd(t0,mipd0);
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTilpp0 = [ubppEF0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        y1dTil0 = y1dTil(t0);
        lmb0 = hTilMF(t0);
        
        ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmb0));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmb0.'];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        lmbln = lmbl0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            lmblnM10 = lmbln;
            zlTilnM10 = [ubTilUnM10; lmblnM10];
            
            %Integração:
            epsilon = argumentos.epsNewtRaph;
            argumentos.epsNewtRaph = epsilon*10^-12;
            EtaU0 = zlTilnM10;
            PsiU = @(eta)(KuzNmk(dtnM1)*eta + FnzNmk(tnM1,eta) +...
                fiuNmkzn(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuzNmk(dtnM1) + dFnzdzT(eta));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            zTilnM1 = EtaConvU;
            ubTilUnM1 = Aulz.'*zTilnM1;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            lmblnM1 = Almblz.'*zTilnM1;
            argumentos.epsNewtRaph = epsilon;
            
            XunM1 = QsiDifUu(ubTilUnM1,ubTilpUnM1);
            y2nM1 = y2(XunM1);
            yp2nM1 = y2p(XunM1);
            y2p2nM1 = y2pp(XunM1);
            y3p2nM1 = y2ppp(XunM1);
            y4p2nM1 = y2pppp(XunM1);
            
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            zlTilnM10 = [xTilnM10; lmbnM10];
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            Eta0 = zlTilnM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Befz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Befz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                 dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            xTilnM1 = Axz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            lmbnM1 = Almbz.'*znM1;
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmbnM1.'];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
        end
    case 27
        uTilTx0 = zeros(Ng-1,1);
        uTilpTx0 = zeros(Ng-1,1);
        uTilppTx0 = zeros(Ng-1,1);
        lmb0 = zeros(rEstr+1,1);
        ZZ = zeros(nT,3*(Ng-1)+rEstr+1+3*1);
        ZZ(1,:) = [uTilTx0.', uTilpTx0.', uTilppTx0.', lmb0.',fTetaX1(t0),fTetapX1(t0),fTetappX1(t0)];
        uTilTxn = uTilTx0;
        uTilpTxn = uTilpTx0;
        uTilppTxn = uTilppTx0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            uTilTxnM10 = uTilTxn + dtnM1*uTilpTxn + 0.5*(dtnM1*uTilppTxn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [uTilTxnM10; lmbnM10];
            %Integração:
            Eta0 = znM10;
            Psi = @(eta)(Knmkz(dtnM1)*eta +...
                fcTilznM1(dtnM1,tnM1,eta,fTetaX1(tnM1),fTetapX1(tnM1),fTetappX1(tnM1),uTilTxn,uTilpTxn,uTilppTxn) +...
                KnmkTxz(fTetaX1(tnM1),fTetapX1(tnM1),fTetappX1(tnM1)) +...
                fiNmkzn(dtnM1,uTilTxn,uTilpTxn,uTilppTxn));
            dPsidEtaT = @(eta)(Knmkz(dtnM1) +...
                dfcTilzdzT(dtnM1,eta,fTetaX1(tnM1),fTetapX1(tnM1),fTetappX1(tnM1),uTilTxn,uTilpTxn,uTilppTxn));
            %{
            ad1 = Knmkz(dtnM1);
            ad2 = dfcTilzdzT(dtnM1,Eta0,fTetaX1(tnM1),fTetapX1(tnM1),fTetappX1(tnM1),uTilTxn,uTilpTxn,uTilppTxn);
            %}
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            uTilTxnM1 = Aurz.'*znM1;
            lmbnM1 = Almbz.'*znM1;
            uTilpTxnM1 = fupnM1(dtnM1,uTilTxnM1,uTilTxn,uTilpTxn,uTilppTxn);
            uTilppTxnM1 = fuppnM1(dtnM1,uTilTxnM1,uTilTxn,uTilpTxn,uTilppTxn);
            %Limitação de variáveis cíclicas:
            if ind_LVC
                %Limitação das variáveis cíclicas:
                tetaX1 = IdNg(6*1,:)*ubTilnM1;
                ubTilnM1 = ubTilnM1 - 2*pi*IdLVC*double(tetaX1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [uTilTxnM1.', uTilpTxnM1.', uTilppTxnM1.', lmbnM1.', fTetaX1(tnM1), fTetapX1(tnM1), fTetappX1(tnM1)];
            %Atualização:
            uTilTxn = uTilTxnM1;
            uTilpTxn = uTilpTxnM1;
            uTilppTxn = uTilppTxnM1;
        end
    case 28
        uTilTx0 = zeros(Ng-1,1);
        uTilpTx0 = zeros(Ng-1,1);
        uTilppTx0 = zeros(Ng-1,1);
        lmb0 = zeros(rEstr+1,1);
        ZZ = zeros(nT,3*(Ng-1)+rEstr+1+3*1);
        ZZ(1,:) = [uTilTx0.', uTilpTx0.', uTilppTx0.', lmb0.',fTetaX1(t0),fTetapX1(t0),fTetappX1(t0)];
        uTilTxn = uTilTx0;
        uTilpTxn = uTilpTx0;
        uTilppTxn = uTilppTx0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            uTilTxnM10 = uTilTxn + dtnM1*uTilpTxn + 0.5*(dtnM1*uTilppTxn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [uTilTxnM10; lmbnM10];
            %Integração:
            Eta0 = znM10;
            Psi = @(eta)(Knmkz(dtnM1)*eta +...
                fcTilznM1(dtnM1,tnM1,eta,fTetaX1(tnM1),fTetapX1(tnM1),fTetappX1(tnM1),uTilTxn,uTilpTxn,uTilppTxn) +...
                KnmkTxz(fTetaX1(tnM1),fTetapX1(tnM1),fTetappX1(tnM1)) +...
                fiNmkzn(dtnM1,uTilTxn,uTilpTxn,uTilppTxn));
            dPsidEtaT = @(eta)(Knmkz(dtnM1) +...
                dfcTilzdzT(dtnM1,eta,fTetaX1(tnM1),fTetapX1(tnM1),fTetappX1(tnM1),uTilTxn,uTilpTxn,uTilppTxn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            uTilTxnM1 = Aurz.'*znM1;
            lmbnM1 = Almbz.'*znM1;
            uTilpTxnM1 = fupnM1(dtnM1,uTilTxnM1,uTilTxn,uTilpTxn,uTilppTxn);
            uTilppTxnM1 = fuppnM1(dtnM1,uTilTxnM1,uTilTxn,uTilpTxn,uTilppTxn);
            %Limitação de variáveis cíclicas:
            if ind_LVC
                %Limitação das variáveis cíclicas:
                tetaX1 = IdNg(6*1,:)*ubTilnM1;
                ubTilnM1 = ubTilnM1 - 2*pi*IdLVC*double(tetaX1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [uTilTxnM1.', uTilpTxnM1.', uTilppTxnM1.', lmbnM1.',fTetaX1(tnM1),fTetapX1(tnM1),fTetappX1(tnM1)];
            %Atualização:
            uTilTxn = uTilTxnM1;
            uTilpTxn = uTilpTxnM1;
            uTilppTxn = uTilppTxnM1;
        end
    case 29
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        lmbl0 = 0;
        ubTilppU0 = Mu\(-Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0) - kNmk*GamaRu.'*lmbl0);
        
        ubAc0 = zeros(8,1);
        ubpAc0 = zeros(8,1);
        ubppAc0 = zeros(8,1);
        lmbAc0 = zeros(2,1);
        ubppAc0 = M1ac\(-Cac*ubpAc0 - Kac*ubAc0 -...
            fTilc(ubAc0,ubpAc0,ubppAc0) - kNmk*GamaAc.'*lmbAc0);
        
        ubEF0 = zeros(Ng,1);
        ubpEF0 = zeros(Ng,1);
        xTil0 = ubEF0;
        xTilp0 = ubpEF0;
        lmb0 = hTilMF(t0);
        ubppEF0 = zeros(Ng,1);
        ubppEF0 = M1gEF\(-CgEF*ubpEF0 - KgEF*ubEF0 -...
            fcEFb(ubEF0,ubpEF0,ubppEF0) - kNmk*GamaMF.'*lmb0);
        xTilpp0 = ubppEF0;
        
        ZZ = zeros(nT,3*length(ubTilU0) + length(lmbl0) + 3*length(ubAc0) + length(lmbAc0) + 3*length(xTil0) + length(lmb0));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', lmbl0, ubAc0.', ubpAc0.', ubppAc0.', lmbAc0.', xTil0.', xTilp0.', xTilpp0.', lmb0.'];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        lmbln = lmbl0;
        
        ubAcn = ubAc0;
        ubpAcn = ubpAc0;
        ubppAcn = ubppAc0;
        lmbAcn = lmbAc0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            
            %Modelo PC longitudinal:
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            lmblnM10 = lmbln;
            zlTilnM10 = [ubTilUnM10; lmblnM10];
            
            epsilon = argumentos.epsNewtRaph;
            argumentos.epsNewtRaph = epsilon*10^0;
            EtaU0 = zlTilnM10;
            PsiU = @(eta)(KuzNmk(dtnM1)*eta + FnzNmk(tnM1,eta) +...
                fiuNmkzn(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuzNmk(dtnM1) + dFnzdzT(eta));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            zTilnM1 = EtaConvU;
            ubTilUnM1 = Aulz.'*zTilnM1;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            lmblnM1 = Almblz.'*zTilnM1;
            argumentos.epsNewtRaph = epsilon;
            
            %Modelo PC com acoplamento:
            ubAcnM10 = ubAcn + dtnM1*ubpAcn + 0.5*(dtnM1*ubppAcn)*dtnM1;
            lmbAcnM10 = lmbAcn;
            zAcnM10 = [ubAcnM10; lmbAcnM10];
            
            Eta0 = zAcnM10;
            Psi = @(eta)(KacNmkz(dtnM1)*eta +...
                fTilcuAcznM1(dtnM1,tnM1,eta,ubAcn,ubpAcn,ubppAcn) +...
                fiNmkuAczn(dtnM1,ubAcn,ubpAcn,ubppAcn));
            dPsidEtaT = @(eta)(KacNmkz(dtnM1) +...
                dfTilcAczdzT(dtnM1,eta,ubAcn,ubpAcn,ubppAcn));
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            zAcnM1 = EtaConv;
            ubAcnM1 = AuAcz.'*zAcnM1;
            ubpAcnM1 = fupnM1(dtnM1,ubAcnM1,ubAcn,ubpAcn,ubppAcn);
            ubppAcnM1 = fuppnM1(dtnM1,ubAcnM1,ubAcn,ubpAcn,ubppAcn);
            lmbAcnM1 = AlmbAcz.'*zAcnM1;
            
            %Modelo EF:
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [xTilnM10; lmbnM10];
            
            Eta0 = znM10;
            Psi = @(eta)(Knmkxz(dtnM1)*eta +...
                fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) +...
                fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn));
            dPsidEtaT = @(eta)(Knmkxz(dtnM1) +...
                dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn));
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            xTilnM1 = Axz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            lmbnM1 = Almbz.'*znM1;
            
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', lmblnM1, ubAcnM1.', ubpAcnM1.', ubppAcnM1.', lmbAcnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', lmbnM1.'];
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            lmbln = lmblnM1;
            ubAcn = ubAcnM1;
            ubpAcn = ubpAcnM1;
            ubppAcn = ubppAcnM1;
            lmbAcn = lmbAcnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
        end
    case 30
        ubAc0 = zeros(8,1);
        ubpAc0 = zeros(8,1);
        ubppAc0 = zeros(8,1);
        lmbAc0 = hTilAc(t0);
        ubppAc0 = M1ac\(-Cac*ubpAc0 - Kac*ubAc0 -...
            fTilc(ubAc0,ubpAc0,ubppAc0) - kNmk*GamaAc.'*lmbAc0);
        
        ubEF0 = zeros(Ng,1);
        ubpEF0 = zeros(Ng,1);
        lmb0 = hTilMF(t0);
        ubppEF0 = zeros(Ng,1);
        ubppEF0 = M1gEF\(-CgEF*ubpEF0 - KgEF*ubEF0 -...
            fcEFb(ubEF0,ubpEF0,ubppEF0) - kNmk*GamaMF.'*lmb0);
        
        ubrk0 = zeros(cPk,1);
        ubprk0 = zeros(cPk,1);
        Tx0 = fTx(t0);
        ub10 = fub1IN(t0);
        ubp10 = fubp1IN(t0);
        ubpp10 = fubpp1IN(t0);
        ubpprk0 = M1s\(Bs*Tx0 - Cs*ubprk0 - fKs(ub10,ubrk0)*ubrk0 -...
            M1s1*ubpp10 - Cs1*ubp10 - fKs1(ub10,ubrk0)*ub10);
        
        ZZ = zeros(nT,3*length(ubAc0) + length(lmbAc0) + 3*length(ubEF0) + length(lmb0) + 3*length(ubrk0) + 3*length(ub10) + length(Tx0));
        ZZ(1,:) = [ubAc0.', ubpAc0.', ubppAc0.', lmbAc0.', ubEF0.', ubpEF0.', ubppEF0.', lmb0.', ubrk0.', ubprk0.', ubpprk0.', ub10, ubp10, ubpp10, Tx0];
        
        ubAcn = ubAc0;
        ubpAcn = ubpAc0;
        ubppAcn = ubppAc0;
        lmbAcn = lmbAc0;
        
        ubEFn = ubEF0;
        ubpEFn = ubpEF0;
        ubppEFn = ubppEF0;
        lmbn = lmb0;
        
        ubrkn = ubrk0;
        ubprkn = ubprk0;
        ubpprkn = ubpprk0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            
            %Modelo PC com acoplamento:
            ubAcnM10 = ubAcn + dtnM1*ubpAcn + 0.5*(dtnM1*ubppAcn)*dtnM1;
            lmbAcnM10 = lmbAcn;
            zAcnM10 = [ubAcnM10; lmbAcnM10];
            
            Eta0 = zAcnM10;
            Psi = @(eta)(KacNmkz(dtnM1)*eta +...
                fTilcuAcznM1(dtnM1,tnM1,eta,ubAcn,ubpAcn,ubppAcn) +...
                fiNmkuAczn(dtnM1,ubAcn,ubpAcn,ubppAcn));
            dPsidEtaT = @(eta)(KacNmkz(dtnM1) +...
                dfTilcAczdzT(dtnM1,eta,ubAcn,ubpAcn,ubppAcn));
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            zAcnM1 = EtaConv;
            ubAcnM1 = AuAcz.'*zAcnM1;
            ubpAcnM1 = fupnM1(dtnM1,ubAcnM1,ubAcn,ubpAcn,ubppAcn);
            ubppAcnM1 = fuppnM1(dtnM1,ubAcnM1,ubAcn,ubpAcn,ubppAcn);
            lmbAcnM1 = AlmbAcz.'*zAcnM1;
            
            %Modelo EF:
            ubEFnM10 = ubEFn + dtnM1*ubpEFn + 0.5*(dtnM1*ubppEFn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [ubEFnM10; lmbnM10];
            
            Eta0 = znM10;
            Psi = @(eta)(KefNmkz(dtnM1)*eta +...
                fTilcEFznM1(dtnM1,tnM1,eta,ubEFn,ubpEFn,ubppEFn) +...
                fiNmkuzn(dtnM1,ubEFn,ubpEFn,ubppEFn));
            dPsidEtaT = @(eta)(KefNmkz(dtnM1) +...
                dfTilcEFzdzT(dtnM1,eta,ubEFn,ubpEFn,ubppEFn));
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            ubEFnM1 = AuEFz.'*znM1;
            ubpEFnM1 = fupnM1(dtnM1,ubEFnM1,ubEFn,ubpEFn,ubppEFn);
            ubppEFnM1 = fuppnM1(dtnM1,ubEFnM1,ubEFn,ubpEFn,ubppEFn);
            lmbnM1 = Almbz.'*znM1;
            
            %Modelo reduzido:
            ubrknM10 = ubrkn + dtnM1*ubprkn + 0.5*(dtnM1*ubpprkn)*dtnM1;
            TxnM1 = fTx(tnM1);
            ub1nM1 = fub1IN(tnM1);
            ubp1nM1 = fubp1IN(tnM1);
            ubpp1nM1 = fubpp1IN(tnM1);
            
            Eta0 = ubrknM10;
            Psi = @(eta)(Khats(dtnM1,tnM1,eta)*eta +...
                fiTilskn(dtnM1,ubrkn,ubprkn,ubpprkn) - Bs*TxnM1 +...
                M1s1*ubpp1nM1 + Cs1*ubp1nM1 + fKs1(ub10,eta)*ub1nM1);
            dPsidEtaT = @(eta)(Khats(dtnM1,tnM1,eta));
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            ubrknM1 = EtaConv;
            ubprknM1 = fupnM1(dtnM1,ubrknM1,ubrkn,ubprkn,ubpprkn);
            ubpprknM1 = fuppnM1(dtnM1,ubrknM1,ubrkn,ubprkn,ubpprkn);
            
            %Armazenamento:
            ZZ(i,:) = [ubAcnM1.', ubpAcnM1.', ubppAcnM1.', lmbAcnM1.', ubEFnM1.', ubpEFnM1.', ubppEFnM1.', lmbnM1.', ubrknM1.', ubprknM1.', ubpprknM1.', ub1nM1, ubp1nM1, ubpp1nM1, TxnM1];
            %Atualização:
            ubAcn = ubAcnM1;
            ubpAcn = ubpAcnM1;
            ubppAcn = ubppAcnM1;
            lmbAcn = lmbAcnM1;
            ubEFn = ubEFnM1;
            ubpEFn = ubpEFnM1;
            ubppEFn = ubppEFnM1;
            lmbn = lmbnM1;
            ubrkn = ubrknM1;
            ubprkn = ubprknM1;
            ubpprkn = ubpprknM1;
        end
    case 31
        ubrk0 = zeros(cPk,1);
        ubprk0 = zeros(cPk,1);
        Tx0 = fTx(t0);
        ubpprk0 = BtX*Tx0 - CtX*ubprk0 - KtX*ubrk0;
        lmb0 = 0;
        ZZ = zeros(nT,3*length(ubrk0) + length(lmb0) + length(Tx0));
        ZZ(1,:) = [ubrk0.', ubprk0.', ubpprk0.', lmb0, Tx0];
        
        ubrkn = ubrk0;
        ubprkn = ubprk0;
        ubpprkn = ubpprk0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            
            ubrknM10 = ubrkn + dtnM1*ubprkn + 0.5*(dtnM1*ubpprkn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [ubrknM10; lmbnM10];
            TxnM1 = fTx(tnM1);
            
            Eta0 = znM10;
            Psi = @(eta)(Knmk(dtnM1,eta)*eta +...
                fiTilrkzn(dtnM1,tnM1,ubrkn,ubprkn,ubpprkn) - Bkz*TxnM1);
            dPsidEtaT = @(eta)(Knmk(dtnM1,eta));
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            ubrknM1 = Auz.'*znM1;
            ubprknM1 = fupnM1(dtnM1,ubrknM1,ubrkn,ubprkn,ubpprkn);
            ubpprknM1 = fuppnM1(dtnM1,ubrknM1,ubrkn,ubprkn,ubpprkn);
            lmbnM1 = Almbz.'*znM1;
            
            %Armazenamento:
            ZZ(i,:) = [ubrknM1.', ubprknM1.', ubpprknM1.', lmbnM1, TxnM1];
            %Atualização:
            ubrkn = ubrknM1;
            ubprkn = ubprknM1;
            ubpprkn = ubpprknM1;
            lmbn = lmbnM1;
        end
    case 32
        ubrk0 = zeros(cPk,1);
        ubprk0 = zeros(cPk,1);
        Tx0 = fTx(t0);
        ub10 = fub1IN(t0);
        ubp10 = fubp1IN(t0);
        ubpp10 = fubpp1IN(t0);
        ubpprk0 = M1s\(Bs*Tx0 - Cs*ubprk0 - fKs(ub10,ubrk0)*ubrk0 -...
            M1s1*ubpp10 - Cs1*ubp10 - fKs1(ub10,ubrk0)*ub10);
        ZZ = zeros(nT,3*length(ubrk0) + 3*length(ub10) + length(Tx0));
        ZZ(1,:) = [ubrk0.', ubprk0.', ubpprk0.', ub10, ubp10, ubpp10, Tx0];
        
        ubrkn = ubrk0;
        ubprkn = ubprk0;
        ubpprkn = ubpprk0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            
            ubrknM10 = ubrkn + dtnM1*ubprkn + 0.5*(dtnM1*ubpprkn)*dtnM1;
            TxnM1 = fTx(tnM1);
            ub1nM1 = fub1IN(tnM1);
            ubp1nM1 = fubp1IN(tnM1);
            ubpp1nM1 = fubpp1IN(tnM1);
            
            Eta0 = ubrknM10;
            Psi = @(eta)(Khats(dtnM1,tnM1,eta)*eta +...
                fiTilskn(dtnM1,ubrkn,ubprkn,ubpprkn) - Bs*TxnM1 +...
                M1s1*ubpp1nM1 + Cs1*ubp1nM1 + fKs1(ub10,eta)*ub1nM1);
            dPsidEtaT = @(eta)(Khats(dtnM1,tnM1,eta));
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            ubrknM1 = EtaConv;
            ubprknM1 = fupnM1(dtnM1,ubrknM1,ubrkn,ubprkn,ubpprkn);
            ubpprknM1 = fuppnM1(dtnM1,ubrknM1,ubrkn,ubprkn,ubpprkn);
            
            %Armazenamento:
            ZZ(i,:) = [ubrknM1.', ubprknM1.', ubpprknM1.', ub1nM1, ubp1nM1, ubpp1nM1, TxnM1];
            %Atualização:
            ubrkn = ubrknM1;
            ubprkn = ubprknM1;
            ubpprkn = ubpprknM1;
        end
    case 33
        fZp = @(t,Z)(fZpModAlt( t, Z, argumentos, repositorio ));
        tspan = T;
        Z0 = zeros(6,1);
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [T,ZZ] = ode45(fZp,tspan,Z0,options);
    case 34
        fZp = @(t,Z)(fZpModAlt23( t, Z, argumentos, repositorio ));
        tspan = T;
        Z0 = zeros(4,1);
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [T,ZZ] = ode45(fZp,tspan,Z0,options);
    case 35
        ubEF0 = zeros(Ng,1);
        ubpEF0 = zeros(Ng,1);
        y20 = FnbNM1(IdNg(Ng-5,:)*ubEF0);
        yp20 = 0;
        y2p20 = 0;
        y3p20 = 0;
        y4p20 = 0;
        
        uTilT0 = Atu.'*ubEF0;
        uTilpT0 = Atu.'*ubpEF0;
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        xTil0 = [ubEF0; mid0; fiBdT0];
        mipd0 = fmipd(t0,mid0);
        xTilp0 = [ubpEF0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        ubppEF0 = zeros(Ng,1);
        ubppEF0 = M1gEF\(BefTx*Tx0 - CgEF*ubpEF0 - KgEF*ubEF0 -...
            fcEFb(ubEF0,ubpEF0,ubppEF0));
        
        mippd0 = fmippd(t0,mipd0);
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTilpp0 = [ubppEF0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        y1dTil0 = y1dTil(t0);
        lmb0 = hTilMF(t0);
        
        ZZ = zeros(nT,3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmb0));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmb0.'];
        
        y2n = y20;
        yp2n = yp20;
        y2p2n = y2p20;
        y3p2n = y3p20;
        y4p2n = y4p20;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            zlTilnM10 = [xTilnM10; lmbnM10];
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            ubEFnM10 = fubTilEF(xTilnM10);
            ubNM1nM10 = IdNg(Ng-5,:)*ubEFnM10;
            y2nM10 = kNM1*ubNM1nM10*double(ubNM1nM10 >= 0);
            xTilpnM10 = fupnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            ubpEFnM10 = fubTilEF(xTilpnM10);
            ubpNM1nM10 = IdNg(Ng-5,:)*ubpEFnM10;
            yp2nM10 = kNM1*ubpNM1nM10*double(ubNM1nM10 >= 0);
            xTilppnM10 = fuppnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            ubppEFnM10 = fubTilEF(xTilppnM10);
            ubppNM1nM10 = IdNg(Ng-5,:)*ubppEFnM10;
            y2p2nM10 = kNM1*ubppNM1nM10*double(ubNM1nM10 >= 0);
            y3p2nM1 = y3p2n;
            y4p2nM1 = y4p2n;
            
            Eta0 = zlTilnM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Befz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM10,yp2nM10,y2p2nM10,y3p2nM1,y4p2nM1,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM10,yp2nM10,y2p2nM10,y3p2nM1,ZnM10);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM10,yp2nM10,y2p2nM10,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Befz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM10,yp2nM10,y2p2nM10,y3p2nM1,y4p2nM1,ZnM10);
                 dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM10,yp2nM10,y2p2nM10,y3p2nM1,ZnM10);
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM10,yp2nM10,y2p2nM10,y3p2nM1,y4p2nM1,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            xTilnM1 = Axz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            lmbnM1 = Almbz.'*znM1;
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            ubEFnM1 = fubTilEF(xTilnM1);
            ubNM1nM1 = IdNg(Ng-5,:)*ubEFnM1;
            y2nM1 = kNM1*ubNM1nM1*double(ubNM1nM1 >= 0);
            ubpEFnM1 = fubTilEF(xTilpnM1);
            ubpNM1nM1 = IdNg(Ng-5,:)*ubpEFnM1;
            yp2nM1 = kNM1*ubpNM1nM1*double(ubNM1nM1 >= 0);
            ubppEFnM1 = fubTilEF(xTilppnM1);
            ubppNM1nM1 = IdNg(Ng-5,:)*ubppEFnM1;
            y2p2nM1 = kNM1*ubppNM1nM1*double(ubNM1nM1 >= 0);
            y3p2nM1 = 0;
            y4p2nM1 = 0;
            
            %Reintegração com valores atualizados de y2 e suas derivadas:
            %Obs.: 
            %(1) Número pré-determinado de iterações
            %(2) Objetivo: melhorar aproximação de y2 e suas derivadas
            Ny = 1;
            ky = 0;
            while ky < Ny
                Eta0 = zlTilnM10;
                Psi = @(eta)(...
                    [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Befz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                     fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                     GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
                dPsidEtaT = @(eta)(...
                    [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Befz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                     dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                     dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
                [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
                znM1 = EtaConv;
                xTilnM1 = Axz.'*znM1;
                xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
                xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
                lmbnM1 = Almbz.'*znM1;
                %Maior precisão na estimação de fippBdnM1:
                midnM1 = Amix.'*xTilnM1;
                mipdnM1 = fmipd(tnM1,midnM1);
                mippdnM1 = fmippd(tnM1,mipdnM1);
                fipBdTnM1 = Afix.'*xTilpnM1;
                fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
                xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];

                ubEFnM1 = fubTilEF(xTilnM1);
                ubNM1nM1 = IdNg(Ng-5,:)*ubEFnM1;
                y2nM1 = kNM1*ubNM1nM1*double(ubNM1nM1 >= 0);
                ubpEFnM1 = fubTilEF(xTilpnM1);
                ubpNM1nM1 = IdNg(Ng-5,:)*ubpEFnM1;
                yp2nM1 = kNM1*ubpNM1nM1*double(ubNM1nM1 >= 0);
                ubppEFnM1 = fubTilEF(xTilppnM1);
                ubppNM1nM1 = IdNg(Ng-5,:)*ubppEFnM1;
                y2p2nM1 = kNM1*ubppNM1nM1*double(ubNM1nM1 >= 0);
                y3p2nM1 = 0;
                y4p2nM1 = 0;
                ky = ky + 1;
            end
            
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmbnM1.'];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
            y2n = y2nM1;
            yp2n = yp2nM1;
            y2p2n = y2p2nM1;
            y3p2n = y3p2nM1;
            y4p2n = y4p2nM1;
        end
    case 36
        ubEF0 = zeros(Ng,1);
        ubpEF0 = zeros(Ng,1);
        y20 = FnbNM1(IdNg(Ng-5,:)*ubEF0);
        yp20 = 0;
        y2p20 = 0;
        y3p20 = 0;
        y4p20 = 0;
        
        uTilT0 = Atu.'*ubEF0;
        uTilpT0 = Atu.'*ubpEF0;
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        xTil0 = [ubEF0; mid0; fiBdT0];
        mipd0 = fmipd(t0,mid0);
        xTilp0 = [ubpEF0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        ubppEF0 = zeros(Ng,1);
        ubppEF0 = M1gEF\(BefTx*Tx0 - CgEF*ubpEF0 - KgEF*ubEF0 -...
            fcEFb(ubEF0,ubpEF0,ubppEF0));
        
        mippd0 = fmippd(t0,mipd0);
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTilpp0 = [ubppEF0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        y1dTil0 = y1dTil(t0);
        lmb0 = hTilMF(t0);
        
        %Para aproximação de y2 e suas derivadas:
        xTil3p0 = zeros(Ng+2,1);
        xTil4p0 = zeros(Ng+2,1);
        KY = zeros(nT,1);
        DY2 = zeros(nT,1);
        
        %Vetor de armazenamento de dados:
        ZZ = zeros(nT,3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmb0));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmb0.'];
        
        y2n = y20;
        yp2n = yp20;
        y2p2n = y2p20;
        y3p2n = y3p20;
        y4p2n = y4p20;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        xTil3pn = xTil3p0;
        xTil4pn = xTil4p0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            zlTilnM10 = [xTilnM10; lmbnM10];
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            ubEFnM10 = fubTilEF(xTilnM10);
            ubNM1nM10 = IdNg(Ng-5,:)*ubEFnM10;
            y2nM10 = kNM1*ubNM1nM10*double(ubNM1nM10 >= 0);
            xTilpnM10 = fupnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            ubpEFnM10 = fubTilEF(xTilpnM10);
            ubpNM1nM10 = IdNg(Ng-5,:)*ubpEFnM10;
            yp2nM10 = kNM1*ubpNM1nM10*double(ubNM1nM10 >= 0);
            xTilppnM10 = fuppnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            ubppEFnM10 = fubTilEF(xTilppnM10);
            ubppNM1nM10 = IdNg(Ng-5,:)*ubppEFnM10;
            y2p2nM10 = kNM1*ubppNM1nM10*double(ubNM1nM10 >= 0);
            xTil3pnM10 = fupnM1(dtnM1,xTilppnM10,xTilppn,xTil3pn,xTil4pn);
            ub3pEFnM10 = fubTilEF(xTil3pnM10);
            ub3pNM1nM10 = IdNg(Ng-5,:)*ub3pEFnM10;
            y3p2nM10 = kNM1*ub3pNM1nM10*double(ubNM1nM10 >= 0);
            xTil4pnM10 = fuppnM1(dtnM1,xTilppnM10,xTilppn,xTil3pn,xTil4pn);
            ub4pEFnM10 = fubTilEF(xTil4pnM10);
            ub4pNM1nM10 = IdNg(Ng-5,:)*ub4pEFnM10;
            y4p2nM10 = kNM1*ub4pNM1nM10*double(ubNM1nM10 >= 0);
            
            Eta0 = zlTilnM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Befz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM10,yp2nM10,y2p2nM10,y3p2nM10,y4p2nM10,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM10,yp2nM10,y2p2nM10,y3p2nM10,ZnM10);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM10,yp2nM10,y2p2nM10,y3p2nM10,y4p2nM10,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Befz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM10,yp2nM10,y2p2nM10,y3p2nM10,y4p2nM10,ZnM10);
                 dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM10,yp2nM10,y2p2nM10,y3p2nM10,ZnM10);
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM10,yp2nM10,y2p2nM10,y3p2nM10,y4p2nM10,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            xTilnM1 = Axz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            lmbnM1 = Almbz.'*znM1;
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            ubEFnM1 = fubTilEF(xTilnM1);
            ubNM1nM1 = IdNg(Ng-5,:)*ubEFnM1;
            y2nM1 = kNM1*ubNM1nM1*double(ubNM1nM1 >= 0);
            ubpEFnM1 = fubTilEF(xTilpnM1);
            ubpNM1nM1 = IdNg(Ng-5,:)*ubpEFnM1;
            yp2nM1 = kNM1*ubpNM1nM1*double(ubNM1nM1 >= 0);
            ubppEFnM1 = fubTilEF(xTilppnM1);
            ubppNM1nM1 = IdNg(Ng-5,:)*ubppEFnM1;
            y2p2nM1 = kNM1*ubppNM1nM1*double(ubNM1nM1 >= 0);
            xTil3pnM1 = fupnM1(dtnM1,xTilppnM1,xTilppn,xTil3pn,xTil4pn);
            ub3pEFnM1 = fubTilEF(xTil3pnM1);
            ub3pNM1nM1 = IdNg(Ng-5,:)*ub3pEFnM1;
            y3p2nM1 = kNM1*ub3pNM1nM1*double(ubNM1nM1 >= 0);
            xTil4pnM1 = fuppnM1(dtnM1,xTilppnM1,xTilppn,xTil3pn,xTil4pn);
            ub4pEFnM1 = fubTilEF(xTil4pnM1);
            ub4pNM1nM1 = IdNg(Ng-5,:)*ub4pEFnM1;
            y4p2nM1 = kNM1*ub4pNM1nM1*double(ubNM1nM1 >= 0);
            
            %Reintegração com valores atualizados de y2 e suas derivadas:
            %Obs.: 
            %(1) Número pré-determinado máximo de iterações
            %(2) Objetivo: melhorar aproximação de y2 e suas derivadas
            NyMax = 15;
            ky = 1;
            epsDy2 = 10^-15;
            Y20 = [y2nM10, yp2nM10, y2p2nM10, y3p2nM1, y4p2nM1];
            Y2 = [y2nM1, yp2nM1, y2p2nM1, y3p2nM1, y4p2nM1];
            dy2 = sqrt((Y2 - Y20)*(Y2 - Y20).');
            dY2 = zeros(1,1);
            while (dy2 > epsDy2) && (ky <= NyMax)
                Y20 = Y2;
                Eta0 = zlTilnM10;
                Psi = @(eta)(...
                    [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Befz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                     fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                     GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
                dPsidEtaT = @(eta)(...
                    [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Befz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                     dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                     dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
                [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
                znM1 = EtaConv;
                xTilnM1 = Axz.'*znM1;
                xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
                xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
                lmbnM1 = Almbz.'*znM1;
                %Maior precisão na estimação de fippBdnM1:
                midnM1 = Amix.'*xTilnM1;
                mipdnM1 = fmipd(tnM1,midnM1);
                mippdnM1 = fmippd(tnM1,mipdnM1);
                fipBdTnM1 = Afix.'*xTilpnM1;
                fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
                xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];

                ubEFnM1 = fubTilEF(xTilnM1);
                ubNM1nM1 = IdNg(Ng-5,:)*ubEFnM1;
                y2nM1 = kNM1*ubNM1nM1*double(ubNM1nM1 >= 0);
                ubpEFnM1 = fubTilEF(xTilpnM1);
                ubpNM1nM1 = IdNg(Ng-5,:)*ubpEFnM1;
                yp2nM1 = kNM1*ubpNM1nM1*double(ubNM1nM1 >= 0);
                ubppEFnM1 = fubTilEF(xTilppnM1);
                ubppNM1nM1 = IdNg(Ng-5,:)*ubppEFnM1;
                y2p2nM1 = kNM1*ubppNM1nM1*double(ubNM1nM1 >= 0);
                xTil3pnM1 = fupnM1(dtnM1,xTilppnM1,xTilppn,xTil3pn,xTil4pn);
                ub3pEFnM1 = fubTilEF(xTil3pnM1);
                ub3pNM1nM1 = IdNg(Ng-5,:)*ub3pEFnM1;
                y3p2nM1 = kNM1*ub3pNM1nM1*double(ubNM1nM1 >= 0);
                xTil4pnM1 = fuppnM1(dtnM1,xTilppnM1,xTilppn,xTil3pn,xTil4pn);
                ub4pEFnM1 = fubTilEF(xTil4pnM1);
                ub4pNM1nM1 = IdNg(Ng-5,:)*ub4pEFnM1;
                y4p2nM1 = kNM1*ub4pNM1nM1*double(ubNM1nM1 >= 0);
                
                Y2 = [y2nM1, yp2nM1, y2p2nM1, y3p2nM1, y4p2nM1];
                dy2 = sqrt((Y2 - Y20)*(Y2 - Y20).');
                dY2(ky,1) = dy2;
                ky = ky + 1
            end
            KY(i,1) = ky;
            DY2(i,1) = dy2;
            
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmbnM1.'];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            xTil3pn = xTil3pnM1;
            xTil4pn = xTil4pnM1;
            lmbn = lmbnM1;
            y2n = y2nM1;
            yp2n = yp2nM1;
            y2p2n = y2p2nM1;
            y3p2n = y3p2nM1;
            y4p2n = y4p2nM1;
        end
    case 37
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        ubTilppU0 = Mu\(-Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        lmbl0 = hlU(t0);
        Xu0 = QsiDifUu(ubTilU0,ubTilpU0);
        y20 = y2(Xu0);
        yp20 = y2p(Xu0);
        y2p20 = y2pp(Xu0);
        y3p20 = y2ppp(Xu0);
        y4p20 = y2pppp(Xu0);
        
        ubAc0 = zeros(8,1);
        ubpAc0 = zeros(8,1);
        uTilT0 = Atu.'*ubAc0;
        uTilpT0 = Atu.'*ubpAc0;
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        xTil0 = [ubAc0; mid0; fiBdT0];
        mipd0 = fmipd(t0,mid0);
        xTilp0 = [ubpAc0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        ubppAc0 = zeros(8,1);
        ubppAc0 = M1ac\(Bac(:,1)*Tx0 - Cac*ubpAc0 - Kac*ubAc0 -...
            fTilc(ubAc0,ubpAc0,ubppAc0));
        
        mippd0 = fmippd(t0,mipd0);
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTilpp0 = [ubppAc0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        y1dTil0 = y1dTil(t0);
        lmb0 = hlU(t0);
        
        ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmbl0) + length(lmb0));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmbl0, lmb0];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        lmbln = lmbl0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            lmblnM10 = lmbln;
            zlTilnM10 = [ubTilUnM10; lmblnM10];
            
            %Integração:
            epsilon = argumentos.epsNewtRaph;
            argumentos.epsNewtRaph = epsilon*10^-12;
            EtaU0 = zlTilnM10;
            PsiU = @(eta)(KuzNmk(dtnM1)*eta + FnzNmk(tnM1,eta) +...
                fiuNmkzn(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuzNmk(dtnM1) + dFnzdzT(eta));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            zTilnM1 = EtaConvU;
            ubTilUnM1 = Aulz.'*zTilnM1;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            lmblnM1 = Almblz.'*zTilnM1;
            argumentos.epsNewtRaph = epsilon;
            
            XunM1 = QsiDifUu(ubTilUnM1,ubTilpUnM1);
            y2nM1 = y2(XunM1);
            yp2nM1 = y2p(XunM1);
            y2p2nM1 = y2pp(XunM1);
            y3p2nM1 = y2ppp(XunM1);
            y4p2nM1 = y2pppp(XunM1);
            
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [xTilnM10; lmbnM10];
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            Eta0 = znM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Bacz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bacz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                 dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            xTilnM1 = Axz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            lmbnM1 = Almbz.'*znM1;
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([8+1,8+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmblnM1, lmbnM1];
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            lmbln = lmblnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
        end
    case 38
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        ubTilppU0 = Mu\(-Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        lmbl0 = hU(t0);
        Xu0 = QsiDifUu(ubTilU0,ubTilpU0);
        y20 = y2(Xu0);
        yp20 = y2p(Xu0);
        y2p20 = y2pp(Xu0);
        y3p20 = y2ppp(Xu0);
        y4p20 = y2pppp(Xu0);
        
        ubAc0 = zeros(Ng,1);
        ubpAc0 = zeros(Ng,1);
        uTilT0 = Atu.'*ubAc0;
        uTilpT0 = Atu.'*ubpAc0;
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        xTil0 = [ubAc0; mid0; fiBdT0];
        mipd0 = fmipd(t0,mid0);
        xTilp0 = [ubpAc0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        ubppAc0 = zeros(Ng,1);
        ubppAc0 = M1gEF\(BefTx*Tx0 - CgEF*ubpAc0 - KgEF*ubAc0 -...
            fcEFb(ubAc0,ubpAc0,ubppAc0));
        
        mippd0 = fmippd(t0,mipd0);
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTilpp0 = [ubppAc0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        y1dTil0 = y1dTil(t0);
        lmb0 = hTilAc(t0);
        
        ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmbl0) + length(lmb0));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmbl0, lmb0.'];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        lmbln = lmbl0;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            lmblnM10 = lmbln;
            zlTilnM10 = [ubTilUnM10; lmblnM10];
            
            %Integração:
            epsilon = argumentos.epsNewtRaph;
            argumentos.epsNewtRaph = epsilon*10^-12;
            EtaU0 = zlTilnM10;
            PsiU = @(eta)(KuzNmk(dtnM1)*eta + FnzNmk(tnM1,eta) +...
                fiuNmkzn(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuzNmk(dtnM1) + dFnzdzT(eta));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            zTilnM1 = EtaConvU;
            ubTilUnM1 = Aulz.'*zTilnM1;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            lmblnM1 = Almblz.'*zTilnM1;
            argumentos.epsNewtRaph = epsilon;
            
            XunM1 = QsiDifUu(ubTilUnM1,ubTilpUnM1);
            y2nM1 = y2(XunM1);
            yp2nM1 = y2p(XunM1);
            y2p2nM1 = y2pp(XunM1);
            y3p2nM1 = y2ppp(XunM1);
            y4p2nM1 = y2pppp(XunM1);
            
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [xTilnM10; lmbnM10];
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            Eta0 = znM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Bacz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bacz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                 dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            xTilnM1 = Axz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            lmbnM1 = Almbz.'*znM1;
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmblnM1, lmbnM1.'];
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            lmbln = lmblnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
        end
    case 39
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        ubTilppU0 = Mu\(-Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        lmbl0 = hU(t0);
        Xu0 = QsiDifUu(ubTilU0,ubTilpU0);
        y20 = y2(Xu0);
        yp20 = y2p(Xu0);
        y2p20 = y2pp(Xu0);
        y3p20 = y2ppp(Xu0);
        y4p20 = y2pppp(Xu0);
        
        ubAc0 = zeros(Ng,1);
        ubpAc0 = zeros(Ng,1);
        uTilT0 = Atu.'*ubAc0;
        uTilpT0 = Atu.'*ubpAc0;
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        xTil0 = [ubAc0; mid0; fiBdT0];
        mipd0 = fmipd(t0,mid0);
        xTilp0 = [ubpAc0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        ubppAc0 = zeros(Ng,1);
        ubppAc0 = M1gEF\(BefTx*Tx0 - CgEF*ubpAc0 - KgEF*ubAc0 -...
            fcEFb(ubAc0,ubpAc0,ubppAc0));
        
        mippd0 = fmippd(t0,mipd0);
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTilpp0 = [ubppAc0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        y1dTil0 = y1dTil(t0);
        lmb0 = hTilAc(t0);
        %y10 = beta0*Z0;
        
        TetaTil0 = AtetaEux.'*xTil0;
        TetapTil0 = AtetaEux.'*xTilp0;
        TetappTil0 = AtetaEux.'*xTilpp0;
        FiTil0 = zeros(3,1);
        FiTil0(2,1) = -10^-5; %Se for nulo, ocorre inversão de matriz sing.
        FipTil0 = fFipTil(FiTil0,TetaTil0,TetapTil0);
        FippTil0 = fFippTil(FiTil0,FipTil0,TetaTil0,TetapTil0,TetappTil0);
        
        %ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmbl0) + length(lmb0) + 3*length(FiTil0) + length(y10));
        %ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmbl0, lmb0.', FiTil0.', FipTil0.', FippTil0.', y10];
        ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmbl0) + length(lmb0) + 3*length(FiTil0));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmbl0, lmb0.', FiTil0.', FipTil0.', FippTil0.'];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        lmbln = lmbl0;
        %y1n = y10;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        
        FiTiln = FiTil0;
        FipTiln = FipTil0;
        FippTiln = FippTil0;
        
        %Y1 = zeros(NJ+1,1);
        %Yb1 = zeros(NJ+1,1);
        %y1nM1 = y1n;
        %yb1n = y1n;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            lmblnM10 = lmbln;
            zlTilnM10 = [ubTilUnM10; lmblnM10];
            %if i-1<=NJ-1
            %    yb1nM1 = ((i-1)*yb1n + y1n)/(i);
            %    yb1n = yb1nM1;
            %    Y1(i,1) = y1nM1;
            %    Yb1(i,1) = yb1nM1;
            %else
            %    yb1nM1 = yb1n + (y1nM1 - Y1(1,1))/NJ;
            %    yb1n = yb1nM1;
            %    Y1(1:NJ,1) = Y1(2:NJ+1,1);
            %    Y1(NJ+1,1) = y1nM1;
            %    Yb1(1:NJ,1) = Yb1(2:NJ+1,1);
            %    Yb1(NJ+1,1) = yb1nM1;
            %end
            
            %Integração:
            epsilon = argumentos.epsNewtRaph;
            argumentos.epsNewtRaph = epsilon*10^-12;
            EtaU0 = zlTilnM10;
            PsiU = @(eta)(KuzNmk(dtnM1)*eta + FnzNmk(tnM1,eta) +...
                fiuNmkzn(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuzNmk(dtnM1) + dFnzdzT(eta));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            zTilnM1 = EtaConvU;
            ubTilUnM1 = Aulz.'*zTilnM1;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            lmblnM1 = Almblz.'*zTilnM1;
            argumentos.epsNewtRaph = epsilon;
            
            XunM1 = QsiDifUu(ubTilUnM1,ubTilpUnM1);
            y2nM1 = y2(XunM1);
            yp2nM1 = y2p(XunM1);
            y2p2nM1 = y2pp(XunM1);
            y3p2nM1 = y2ppp(XunM1);
            y4p2nM1 = y2pppp(XunM1);
            
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [xTilnM10; lmbnM10];
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            Eta0 = znM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Bacz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bacz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                 dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            xTilnM1 = Axz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            lmbnM1 = Almbz.'*znM1;
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            %ZnM1 = fZxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %y1nM1 = beta0*ZnM1;
            
            %Vetor velocidade angular (rotor intermediário):
            TetaTilnM1 = AtetaEux.'*xTilnM1;
            TetapTilnM1 = AtetaEux.'*xTilpnM1;
            TetappTilnM1 = AtetaEux.'*xTilppnM1;
            FiTilnM10 = FiTiln + dtnM1*FipTiln + 0.5*(dtnM1*FippTiln)*dtnM1;
            EtaEu0 = FiTilnM10;
            PsiEu = @(eta)(fWce0nM1(dtnM1,eta,FiTiln,FipTiln,FippTiln,TetaTilnM1,TetapTilnM1));
            dPsiEudEtaT = @(eta)(dfWce0dFinM1T(dtnM1,eta,FiTiln,FipTiln,FippTiln));
            [ EtaEuConv ] = f_Newton_Raphson( EtaEu0, PsiEu, dPsiEudEtaT, argumentos );
            FiTilnM1 = EtaEuConv;
            %FipTilnM1 = fupnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln);
            FipTilnM1 = fFipTil(FiTilnM1,TetaTilnM1,TetapTilnM1);
            %Obs.: Estimativa mais precisa de FipTilnM1
            %FippTilnM1 = fuppnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln);
            FippTilnM1 = fFippTil(FiTilnM1,FipTilnM1,TetaTilnM1,TetapTilnM1,TetappTilnM1);
            %Obs.: Estimativa mais precisa de FippTilnM1
            
            %Armazenamento:
            %ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmblnM1, lmbnM1.', FiTilnM1.', FipTilnM1.', FippTilnM1.', yb1nM1];
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmblnM1, lmbnM1.', FiTilnM1.', FipTilnM1.', FippTilnM1.'];
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            lmbln = lmblnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
            %y1n = y1nM1;
            FiTiln = FiTilnM1;
            FipTiln = FipTilnM1;
            FippTiln = FippTilnM1;
        end
    case 40
        ubTilU0 = zeros(3,1);
        ubTilpU0 = zeros(3,1);
        ubTilppU0 = Mu\(-Cu*ubTilpU0 - Ku*ubTilU0 - fNuTil(ubTilU0));
        lmbl0 = 0;
        Xu0 = QsiDifUu(ubTilU0,ubTilpU0);
        y20 = y2(Xu0);
        yp20 = y2p(Xu0);
        y2p20 = y2pp(Xu0);
        y3p20 = y2ppp(Xu0);
        y4p20 = y2pppp(Xu0);
        
        ubAc0 = zeros(Ng,1);
        ubpAc0 = zeros(Ng,1);
        uTilT0 = Atu.'*ubAc0;
        uTilpT0 = Atu.'*ubpAc0;
        Z0 = [uTilT0; uTilpT0];
        mid0 = AmiT*Z0;
        fiBdT0 = KfiT\KlcTd(t0,mid0);
        fipBdT0 = 0;
        xTil0 = [ubAc0; mid0; fiBdT0];
        mipd0 = fmipd(t0,mid0);
        xTilp0 = [ubpAc0; mipd0; fipBdT0];
        Tx0 = UcTx(t0,xTil0,xTilp0,y20,yp20,y2p20,y3p20,y4p20,Z0);
        ubppAc0 = zeros(Ng,1);
        ubppAc0 = M1gEF\(BefTx*Tx0 - CgEF*ubpAc0 - KgEF*ubAc0 -...
            fcEFb(ubAc0,ubpAc0,ubppAc0));
        
        mippd0 = fmippd(t0,mipd0);
        fippBdT0 = ffippBdT(t0,mid0,mipd0,fipBdT0);
        
        xTilpp0 = [ubppAc0; mippd0; fippBdT0];
        Zd0 = fZtdx(t0,xTil0,y20,yp20,y2p20,y3p20,Z0);
        s10 = s1uTil(t0,uTilT0,uTilpT0,y20,yp20,y2p20,y3p20);
        y1dTil0 = y1dTil(t0);
        lmb0 = zeros(rEstrT,1);
        y10 = beta0*Z0;
        
        TetaTil0 = AtetaEux.'*xTil0;
        TetapTil0 = AtetaEux.'*xTilp0;
        TetappTil0 = AtetaEux.'*xTilpp0;
        FiTil0 = zeros(3,1);
        FiTil0(2,1) = -10^-5; %Se for nulo, ocorre inversão de matriz sing.
        FipTil0 = fFipTil(FiTil0,TetaTil0,TetapTil0);
        FippTil0 = fFippTil(FiTil0,FipTil0,TetaTil0,TetapTil0,TetappTil0);
        
        ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmbl0) + length(lmb0) + 3*length(FiTil0) + length(y10));
        ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmbl0, lmb0.', FiTil0.', FipTil0.', FippTil0.', y10];
        
        ubTilUn = ubTilU0;
        ubTilpUn = ubTilpU0;
        ubTilppUn = ubTilppU0;
        lmbln = lmbl0;
        y1n = y10;
        
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lmbn = lmb0;
        
        FiTiln = FiTil0;
        FipTiln = FipTil0;
        FippTiln = FippTil0;
        
        Y1 = zeros(NJ+1,1);
        Yb1 = zeros(NJ+1,1);
        y1nM1 = y1n;
        yb1n = y1n;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            ubTilUnM10 = ubTilUn + dtnM1*ubTilpUn + 0.5*(dtnM1*ubTilppUn)*dtnM1;
            lmblnM10 = lmbln;
            zlTilnM10 = [ubTilUnM10; lmblnM10];
            if i-1<=NJ-1
                yb1nM1 = ((i-1)*yb1n + y1n)/(i);
                yb1n = yb1nM1;
                Y1(i,1) = y1nM1;
                Yb1(i,1) = yb1nM1;
            else
                yb1nM1 = yb1n + (y1nM1 - Y1(1,1))/NJ;
                yb1n = yb1nM1;
                Y1(1:NJ,1) = Y1(2:NJ+1,1);
                Y1(NJ+1,1) = y1nM1;
                Yb1(1:NJ,1) = Yb1(2:NJ+1,1);
                Yb1(NJ+1,1) = yb1nM1;
            end
            
            %Integração:
            epsilon = argumentos.epsNewtRaph;
            argumentos.epsNewtRaph = epsilon*10^-12;
            EtaU0 = zlTilnM10;
            PsiU = @(eta)(KuzNmk(dtnM1)*eta + FnzNmk(tnM1,eta,yb1nM1) +...
                fiuNmkzn(dtnM1,ubTilUn,ubTilpUn,ubTilppUn));
            dPsiUdEtaUT = @(eta)(KuzNmk(dtnM1) + dFnzdzT(eta));
            [ EtaConvU ] = f_Newton_Raphson( EtaU0, PsiU, dPsiUdEtaUT, argumentos );
            zTilnM1 = EtaConvU;
            ubTilUnM1 = Aulz.'*zTilnM1;
            ubTilpUnM1 = fupnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            ubTilppUnM1 = fuppnM1(dtnM1,ubTilUnM1,ubTilUn,ubTilpUn,ubTilppUn);
            lmblnM1 = Almblz.'*zTilnM1;
            argumentos.epsNewtRaph = epsilon;
            
            XunM1 = QsiDifUu(ubTilUnM1,ubTilpUnM1);
            y2nM1 = y2(XunM1);
            yp2nM1 = y2p(XunM1);
            y2p2nM1 = y2pp(XunM1);
            y3p2nM1 = y2ppp(XunM1);
            y4p2nM1 = y2pppp(XunM1);
            
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lmbnM10 = lmbn;
            znM10 = [xTilnM10; lmbnM10];
            ZnM10 = fZxnM1(dtnM1,xTilnM10,xTiln,xTilpn,xTilppn);
            
            Eta0 = znM10;
            Psi = @(eta)(...
                [Knmkxz(dtnM1)*eta + fTilcxznM1(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn) - Bacz*UcTxznM1(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10) + fiNmkxzn(dtnM1,xTiln,xTilpn,xTilppn);
                 fmipdznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) - fmipdAtznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 GamaFiBdTxznM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz*eta - KlcTdxznM1(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            dPsidEtaT = @(eta)(...
                [Knmkxz(dtnM1) + dfTilczdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) - Bacz*dUcTxzdzT(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
                 dmipddzT(dtnM1) - fdmipdAtdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
                 dGamaFiTxzdzT(dtnM1,eta,xTiln,xTilpn,xTilppn) + KfiTxz - dKlcTdxzdzT(tnM1,eta,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10)]);
            [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = EtaConv;
            xTilnM1 = Axz.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            lmbnM1 = Almbz.'*znM1;
            %
            %Maior precisão na estimação de fippBdnM1:
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fipBdTnM1 = Afix.'*xTilpnM1;
            fippBdTnM1 = ffippBdT(tnM1,midnM1,mipdnM1,fipBdTnM1);
            xTilppnM1([Ng+1,Ng+2],1) = [mippdnM1; fippBdTnM1];
            %}
            if ind_LVC
                uTilTnM1 = Atx.'*xTilnM1;
                teta1nM1 = Id3(1,:)*uTilTnM1;
                xTilnM1 = xTilnM1 - 2*pi*[0;1;0;0;0;1;0;1;1;0]*double(teta1nM1>=2*pi);
            end
            TxnM1 = UcTxnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,y4p2nM1,ZnM10);
            y1dTilnM1 = y1dTil(tnM1);
            ZdnM1 = fZtdx(tnM1,xTilnM1,y2nM1,yp2nM1,y2p2nM1,y3p2nM1,ZnM10);
            s1nM1 = s1xnM1(tnM1,dtnM1,xTilnM1,xTiln,xTilpn,xTilppn,y2nM1,yp2nM1,y2p2nM1,y3p2nM1);
            ZnM1 = fZxnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            y1nM1 = beta0*ZnM1;
            
            %Vetor velocidade angular (rotor intermediário):
            TetaTilnM1 = AtetaEux.'*xTilnM1;
            TetapTilnM1 = AtetaEux.'*xTilpnM1;
            TetappTilnM1 = AtetaEux.'*xTilppnM1;
            FiTilnM10 = FiTiln + dtnM1*FipTiln + 0.5*(dtnM1*FippTiln)*dtnM1;
            EtaEu0 = FiTilnM10;
            PsiEu = @(eta)(fWce0nM1(dtnM1,eta,FiTiln,FipTiln,FippTiln,TetaTilnM1,TetapTilnM1));
            dPsiEudEtaT = @(eta)(dfWce0dFinM1T(dtnM1,eta,FiTiln,FipTiln,FippTiln));
            [ EtaEuConv ] = f_Newton_Raphson( EtaEu0, PsiEu, dPsiEudEtaT, argumentos );
            FiTilnM1 = EtaEuConv;
            %FipTilnM1 = fupnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln);
            FipTilnM1 = fFipTil(FiTilnM1,TetaTilnM1,TetapTilnM1);
            %Obs.: Estimativa mais precisa de FipTilnM1
            %FippTilnM1 = fuppnM1(dtnM1,FiTilnM1,FiTiln,FipTiln,FippTiln);
            FippTilnM1 = fFippTil(FiTilnM1,FipTilnM1,TetaTilnM1,TetapTilnM1,TetappTilnM1);
            %Obs.: Estimativa mais precisa de FippTilnM1
            
            %Armazenamento:
            ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmblnM1, lmbnM1.', FiTilnM1.', FipTilnM1.', FippTilnM1.', yb1nM1];
            %Atualização:
            ubTilUn = ubTilUnM1;
            ubTilpUn = ubTilpUnM1;
            ubTilppUn = ubTilppUnM1;
            lmbln = lmblnM1;
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lmbn = lmbnM1;
            y1n = y1nM1;
            FiTiln = FiTilnM1;
            FipTiln = FipTilnM1;
            FippTiln = FippTilnM1;
        end
    otherwise
        disp('other value')
end

%%
%Saídas:

switch simulacao
    case 1
        %Newmark:
        uTilT = ZZ(:,1:3);
        uTilpT = ZZ(:,4:6);
        uTilppT = ZZ(:,7:9);
        Tx = ZZ(:,10);
        y1 = ZZ(:,1:6)*C1.';
        teta1 = uTilT(:,1);
        teta2 = uTilT(:,2);
        teta3 = uTilT(:,3);
        tetap1 = uTilpT(:,1);
        tetap2 = uTilpT(:,2);
        tetap3 = uTilpT(:,3);
        tetapp1 = uTilppT(:,1);
        tetapp2 = uTilppT(:,2);
        tetapp3 = uTilppT(:,3);
        %Entrada:
        figure
        plot(T,Tx,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variável de saída:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = tetap3')
        grid
        %Posições angulares:
        figure
        plot(T,teta1,'k',T,teta2,'r',T,teta3,'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Posições angulares (teta)')
        legend('teta1','teta2','teta3')
        grid
        %Velocidades angulares:
        figure
        plot(T,tetap1*(30/pi),'k',T,tetap2*(30/pi),'r',T,tetap3*(30/pi),'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidades angulares (tetap)')
        legend('tetap1','tetap2','tetap3')
        grid
        %Acelerações angulares:
        figure
        plot(T,tetapp1,'k',T,tetapp2,'r',T,tetapp3,'b')
        xlabel('Tempo (s)')
        ylabel('tetapp (rad/s2)')
        title('Acelerações angulares (tetapp)')
        legend('tetapp1','tetapp2','tetapp3')
        grid
        
        %%
        %ODE45:
        if ODEativo
            uTilT = ZZode(:,1:3);
            uTilpT = ZZode(:,4:6);
            y1ode = ZZode*C1.';
            teta1od = uTilT(:,1);
            teta2od = uTilT(:,2);
            teta3od = uTilT(:,3);
            tetap1od = uTilpT(:,1);
            tetap2od = uTilpT(:,2);
            tetap3od = uTilpT(:,3);
            %Variável de saída:
            figure
            plot(T,y1ode*30/pi,'k')
            xlabel('Tempo (s)')
            ylabel('y1 (RPM)')
            title('Variável de saída: y1 = tetap3 [ODE]')
            grid
            %Posições angulares:
            figure
            plot(T,teta1od,'k',T,teta2od,'r',T,teta3od,'b')
            xlabel('Tempo (s)')
            ylabel('teta (rad)')
            title('Posições angulares (teta) [ODE]')
            legend('teta1','teta2','teta3')
            grid
            %Velocidades angulares:
            figure
            plot(T,tetap1od*(30/pi),'k',T,tetap2od*(30/pi),'r',T,tetap3od*(30/pi),'b')
            xlabel('Tempo (s)')
            ylabel('tetap (RPM)')
            title('Velocidades angulares (tetap) [ODE]')
            legend('tetap1','tetap2','tetap3')
            grid
        end
    case 2
        %Newmark:
        ubTilU = ZZ(:,1:3);
        ubTilpU = ZZ(:,4:6);
        ubTilppU = ZZ(:,7:9);
        Fbu = ZZ(:,10);
        u1 = ubTilU(:,1);
        u2 = ubTilU(:,2);
        u3 = ubTilU(:,3);
        up1 = ubTilpU(:,1);
        up2 = ubTilpU(:,2);
        up3 = ubTilpU(:,3);
        upp1 = ubTilppU(:,1);
        upp2 = ubTilppU(:,2);
        upp3 = ubTilppU(:,3);
        y2 = fNub3(u3);
        %Entrada:
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        %Variável de saída:
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,u1*10^3,'k',T,u2*10^3,'r',T,u3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,up1*10^3,'k',T,up2*10^3,'r',T,up3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,upp1,'k',T,upp2,'r',T,upp3,'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        if MEANativo
            Tpd = T(nPD:nT,1);
            ubTilpd = ZZpd(:,1:3);
            ubTilpPD = ZZpd(:,4:6);
            ubTilppPD = ZZpd(:,7:9);
            FbuPD = ZZpd(:,10);
            ub1PD = ubTilpd(:,1);
            ub2PD = ubTilpd(:,2);
            ub3PD = ubTilpd(:,3);
            ubp1PD = ubTilpPD(:,1);
            ubp2PD = ubTilpPD(:,2);
            ubp3PD = ubTilpPD(:,3);
            ubpp1PD = ubTilppPD(:,1);
            ubpp2PD = ubTilppPD(:,2);
            ubpp3PD = ubTilppPD(:,3);
            y2PD = fNub3(ub3PD);
            %Entrada:
            figure
            plot(Tpd,FbuPD,'k')
            xlabel('Tempo (s)')
            ylabel('Fbu (N)')
            title('Força de atuação longitudinal [PD]')
            grid
            %Variável de saída:
            figure
            plot(Tpd,y2PD,'k')
            xlabel('Tempo (s)')
            ylabel('y2 (N)')
            title('Variável de saída: y2 = Nr(ub3) [PD]')
            grid
            %Posições longitudinais:
            figure
            plot(Tpd,ub1PD*10^3,'k',Tpd,ub2PD*10^3,'r',Tpd,ub3PD*10^3,'b')
            xlabel('Tempo (s)')
            ylabel('u (mm)')
            title('Posições longitudinais (u) [PD]')
            legend('ub1','ub2','ub3')
            grid
            %Velocidades longitudinais:
            figure
            plot(Tpd,ubp1PD*10^3,'k',Tpd,ubp2PD*10^3,'r',Tpd,ubp3PD*10^3,'b')
            xlabel('Tempo (s)')
            ylabel('up (mm/s)')
            title('Velocidades longitudinais (up) [PD]')
            legend('up1','up2','up3')
            grid
            %Acelerações longitudinais:
            figure
            plot(Tpd,ubpp1PD,'k',Tpd,ubpp2PD,'r',Tpd,ubpp3PD,'b')
            xlabel('Tempo (s)')
            ylabel('upp (mm/s2)')
            title('Acelerações longitudinais (upp) [PD]')
            legend('upp1','upp2','upp3')
            grid
        end
        %%
        %ODE45:
        if ODEativo
            ubTilU = ZZode(:,1:3);
            ubTilpU = ZZode(:,4:6);
            ub1od = ubTilU(:,1);
            ub2od = ubTilU(:,2);
            ub3od = ubTilU(:,3);
            ubp1od = ubTilpU(:,1);
            ubp2od = ubTilpU(:,2);
            ubp3od = ubTilpU(:,3);
            y2ode = fNub3(ub3od);
            %Variável de saída:
            figure
            plot(T,y2ode,'k')
            xlabel('Tempo (s)')
            ylabel('y2 (N)')
            title('Variável de saída: y2 = Nr(ub3) [ODE]')
            grid
            %Posições longitudinais:
            figure
            plot(T,ub1od*10^3,'k',T,ub2od*10^3,'r',T,ub3od*10^3,'b')
            xlabel('Tempo (s)')
            ylabel('ub (mm)')
            title('Posições lineares (ub) [ODE]')
            legend('ub1','ub2','ub3')
            grid
            %Velocidades angulares:
            figure
            plot(T,ubp1od*10^3,'k',T,ubp2od*10^3,'r',T,ubp3od*10^3,'b')
            xlabel('Tempo (s)')
            ylabel('up (mm/s)')
            title('Velocidades longitudinais (up) [ODE]')
            legend('up1','up2','up3')
            grid
        end
    case 3
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        Y1d = ZZ(:,3*nx+1:3*nx+ny);
        Zd = ZZ(:,3*nx+ny+1:3*nx+ny+nz);
        Mid = xTil(:,1);
        FiBd = xTil(:,2);
        Mipd = xTilp(:,1);
        FipBd = xTilp(:,2);
        Mippd = xTilpp(:,1);
        FippBd = xTilpp(:,2);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        %Camada limites:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k')
        xlabel('Tempo (s)')
        ylabel('FiBd')
        title('Camada limite - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
    case 4
        nFi = length(fiBdU0);
        ny = length(y2dTil0);
        nz = length(Zd0);
        FiBd = ZZ(:,1:nFi);
        FipBd = ZZ(:,nFi+1:2*nFi);
        FippBd = ZZ(:,2*nFi+1:3*nFi);
        Y2d = ZZ(:,3*nFi+1:3*nFi+ny);
        Zd = ZZ(:,3*nFi+ny+1:3*nFi+ny+nz);
        ub1d = Zd(:,1);
        ub2d = Zd(:,2);
        ub3d = Zd(:,3);
        ubp1d = Zd(:,4);
        ubp2d = Zd(:,5);
        ubp3d = Zd(:,6);
        %Camada limites:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k')
        xlabel('Tempo (s)')
        ylabel('FiBd')
        title('Camada limite - Movimento longitudinal')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Longitudinal')
        legend('FiBd','FipBd','FippBd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,ub1d*10^3,'k',T,ub2d*10^3,'r',T,ub3d*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubd (mm)')
        title('Posição desejada - Mov. Longitudinal - Zd')
        legend('ub1d','ub2d','ub3d')
        grid
        figure
        plot(T,ubp1d*10^3,'k',T,ubp2d*10^3,'r',T,ubp3d*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpd (mm/s)')
        title('Velocidade desejada - Mov. Longitudinal - Zd')
        legend('ubp1d','ubp2d','ubp3d')
        grid
        figure
        plot(T,Y2d(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('y2d (N)')
        title('Saída desejada - Y2d')
        grid
        figure
        %plot(T,Y2d(:,1),'k',T,Y2d(:,2),'r',T,Y2d(:,3),'b',T,Y2d(:,4),'b',T,Y2d(:,5),'b',T,Y2d(:,6),'b')
        plot(T,Y2d(:,1),'k',T,Y2d(:,2),'r',T,Y2d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y2d')
        title('Saída desejada e derivadas - Y2d')
        %legend('y2d','yp2d','ypp2d','yppp2d','ypppp2d','yppppp2d')
        legend('y2d (N)','yp2d (N/s)','ypp2d (N/s2)')
        grid
    case 5
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nu = length(uTilT0);
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        Y1d = ZZ(:,3*nx+1:3*nx+ny);
        Zd = ZZ(:,3*nx+ny+1:3*nx+ny+nz);
        Uc = ZZ(:,3*nx+ny+nz+1);
        S1 = ZZ(:,3*nx+ny+nz+2);
        
        UtilT = xTil(:,1:nu);
        Mid = xTil(:,nu+1);
        FiBd = xTil(:,nu+2);
        UtilpT = xTilp(:,1:nu);
        Mipd = xTilp(:,nu+1);
        FipBd = xTilp(:,nu+2);
        UtilppT = xTilpp(:,1:nu);
        Mippd = xTilpp(:,nu+1);
        FippBd = xTilpp(:,nu+2);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        %Vetor deslocamento e derivadas:
        figure
        plot(T,UtilT(:,1),'k',T,UtilT(:,2),'r',T,UtilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        figure
        plot(T,UtilpT(:,1)*30/pi,'k',T,UtilpT(:,2)*30/pi,'r',T,UtilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        figure
        plot(T,UtilppT(:,1),'k',T,UtilppT(:,2),'r',T,UtilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Variável de controle:
        figure
        plot(T,Uc,'k')
        xlabel('Tempo (s)')
        ylabel('Uc (Nm)')
        title('Variável de controle - Tx')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
    case 6
        nx = length(xTil0);
        ny = length(y2dTil0);
        nz = length(Zd0);
        nu = length(uTilU0);
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        Y2d = ZZ(:,3*nx+1:3*nx+ny);
        Zd = ZZ(:,3*nx+ny+1:3*nx+ny+nz);
        Uc = ZZ(:,3*nx+ny+nz+1);
        S2 = ZZ(:,3*nx+ny+nz+2);
        
        UtilU = xTil(:,1:nu);
        FiBd = xTil(:,nu+1);
        UtilpU = xTilp(:,1:nu);
        FipBd = xTilp(:,nu+1);
        UtilppU = xTilpp(:,1:nu);
        FippBd = xTilpp(:,nu+1);
        ub1d = Zd(:,1);
        ub2d = Zd(:,2);
        ub3d = Zd(:,3);
        ubp1d = Zd(:,4);
        ubp2d = Zd(:,5);
        ubp3d = Zd(:,6);
        %Vetor deslocamento e derivadas:
        figure
        plot(T,UtilU(:,1)*10^3,'k',T,UtilU(:,2)*10^3,'r',T,UtilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        title('Deslocamento linear - Longitudinal')
        legend('ub1','ub2','ub3')
        grid
        figure
        plot(T,UtilpU(:,1)*10^3,'k',T,UtilpU(:,2)*10^3,'r',T,UtilpU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubp (mm)')
        title('Velocidade linear - Longitudinal')
        legend('ubp1','ubp2','ubp3')
        grid
        figure
        plot(T,UtilppU(:,1)*10^3,'k',T,UtilppU(:,2)*10^3,'r',T,UtilppU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        title('Aceleração linear - Longitudinal')
        legend('ubpp1','ubpp2','ubpp3')
        grid
        %Variável de controle:
        figure
        plot(T,Uc,'k')
        xlabel('Tempo (s)')
        ylabel('Uc (N)')
        title('Variável de controle - Fu')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S2,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s2')
        title('Camada limite e variável de escorregamento - Longitudinal')
        legend('FiBdU','s2')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Longitudinal')
        legend('FiBd','FipBd','FippBd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,ub1d*10^3,'k',T,ub2d*10^3,'r',T,ub3d*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubd (mm)')
        title('Posição desejada - Longitudinal - Zd')
        legend('ub1d','ub2d','ub3d')
        grid
        figure
        plot(T,ubp1d*10^3,'k',T,ubp2d*10^3,'r',T,ubp3d*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpd (mm/s)')
        title('Velocidade desejada - Longitudinal - Zd')
        legend('ubp1d','ubp2d','ubp3d')
        grid
        figure
        plot(T,Y2d(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('y2d (N)')
        title('Saída desejada - Y2d')
        grid
        figure
        %plot(T,Y2d(:,1),'k',T,Y2d(:,2),'r',T,Y2d(:,3),'b',T,Y2d(:,4),'b',T,Y2d(:,5),'b',T,Y2d(:,6),'b')
        plot(T,Y2d(:,1),'k',T,Y2d(:,2),'r',T,Y2d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y2d')
        title('Saída desejada e derivadas - Y2d')
        %legend('y2d','yp2d','ypp2d','yppp2d','ypppp2d')
        legend('y2d','yp2d','ypp2d')
        grid
    case 7
        nx = length(xTil0);
        ny = length(y2dTil0);
        nz = length(Zd0);
        nu = length(uTilU0);
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        Y2d = ZZ(:,3*nx+1:3*nx+ny);
        Zd = ZZ(:,3*nx+ny+1:3*nx+ny+nz);
        Uc = ZZ(:,3*nx+ny+nz+1);
        S2 = ZZ(:,3*nx+ny+nz+2);
        
        UtilU = xTil(:,1:nu);
        FiBd = xTil(:,nu+1);
        UtilpU = xTilp(:,1:nu);
        FipBd = xTilp(:,nu+1);
        UtilppU = xTilpp(:,1:nu);
        FippBd = xTilpp(:,nu+1);
        ub1d = Zd(:,1);
        ub2d = Zd(:,2);
        ub3d = Zd(:,3);
        ubp1d = Zd(:,4);
        ubp2d = Zd(:,5);
        ubp3d = Zd(:,6);
        %Vetor deslocamento e derivadas:
        figure
        plot(T,UtilU(:,1)*10^3,'k',T,UtilU(:,2)*10^3,'r',T,UtilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        title('Deslocamento linear - Longitudinal')
        legend('ub1','ub2','ub3')
        grid
        figure
        plot(T,UtilpU(:,1)*10^3,'k',T,UtilpU(:,2)*10^3,'r',T,UtilpU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubp (mm)')
        title('Velocidade linear - Longitudinal')
        legend('ubp1','ubp2','ubp3')
        grid
        figure
        plot(T,UtilppU(:,1)*10^3,'k',T,UtilppU(:,2)*10^3,'r',T,UtilppU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        title('Aceleração linear - Longitudinal')
        legend('ubpp1','ubpp2','ubpp3')
        grid
        %Variável de controle:
        figure
        plot(T,Uc,'k')
        xlabel('Tempo (s)')
        ylabel('Uc (N)')
        title('Variável de controle - Fu')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S2,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s2')
        title('Camada limite e variável de escorregamento - Longitudinal')
        legend('FiBdU','s2')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Longitudinal')
        legend('FiBd','FipBd','FippBd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,ub1d*10^3,'k',T,ub2d*10^3,'r',T,ub3d*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubd (mm)')
        title('Posição desejada - Longitudinal - Zd')
        legend('ub1d','ub2d','ub3d')
        grid
        figure
        plot(T,ubp1d*10^3,'k',T,ubp2d*10^3,'r',T,ubp3d*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpd (mm/s)')
        title('Velocidade desejada - Longitudinal - Zd')
        legend('ubp1d','ubp2d','ubp3d')
        grid
        figure
        plot(T,Y2d(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('y2d (N)')
        title('Saída desejada - Y2d')
        grid
        figure
        %plot(T,Y2d(:,1),'k',T,Y2d(:,2),'r',T,Y2d(:,3),'b',T,Y2d(:,4),'b',T,Y2d(:,5),'b',T,Y2d(:,6),'b')
        plot(T,Y2d(:,1),'k',T,Y2d(:,2),'r',T,Y2d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y2d')
        title('Saída desejada e derivadas - Y2d')
        %legend('y2d','yp2d','ypp2d','yppp2d','ypppp2d')
        legend('y2d','yp2d','ypp2d')
        grid
    case 8
        nx = length(xTil0);
        ny = length(y2dTil0);
        nz = length(Zd0);
        nu = length(uTilU0);
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        Y2d = ZZ(:,3*nx+1:3*nx+ny);
        Zd = ZZ(:,3*nx+ny+1:3*nx+ny+nz);
        Uc = ZZ(:,3*nx+ny+nz+1);
        S2 = ZZ(:,3*nx+ny+nz+2);
        
        UtilU = xTil(:,1:nu);
        FiBd = xTil(:,nu+1);
        UtilpU = xTilp(:,1:nu);
        FipBd = xTilp(:,nu+1);
        UtilppU = xTilpp(:,1:nu);
        FippBd = xTilpp(:,nu+1);
        ub1d = Zd(:,1);
        ub2d = Zd(:,2);
        ub3d = Zd(:,3);
        ubp1d = Zd(:,4);
        ubp2d = Zd(:,5);
        ubp3d = Zd(:,6);
        %Vetor deslocamento e derivadas:
        figure
        plot(T,UtilU(:,1)*10^3,'k',T,UtilU(:,2)*10^3,'r',T,UtilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        title('Deslocamento linear - Longitudinal')
        legend('ub1','ub2','ub3')
        grid
        figure
        plot(T,UtilpU(:,1)*10^3,'k',T,UtilpU(:,2)*10^3,'r',T,UtilpU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubp (mm)')
        title('Velocidade linear - Longitudinal')
        legend('ubp1','ubp2','ubp3')
        grid
        figure
        plot(T,UtilppU(:,1)*10^3,'k',T,UtilppU(:,2)*10^3,'r',T,UtilppU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        title('Aceleração linear - Longitudinal')
        legend('ubpp1','ubpp2','ubpp3')
        grid
        %Variável de controle:
        figure
        plot(T,Uc,'k')
        xlabel('Tempo (s)')
        ylabel('Uc (N)')
        title('Variável de controle - Fu')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S2,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s2')
        title('Camada limite e variável de escorregamento - Longitudinal')
        legend('FiBdU','s2')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Longitudinal')
        legend('FiBd','FipBd','FippBd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,ub1d*10^3,'k',T,ub2d*10^3,'r',T,ub3d*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubd (mm)')
        title('Posição desejada - Longitudinal - Zd')
        legend('ub1d','ub2d','ub3d')
        grid
        figure
        plot(T,ubp1d*10^3,'k',T,ubp2d*10^3,'r',T,ubp3d*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpd (mm/s)')
        title('Velocidade desejada - Longitudinal - Zd')
        legend('ubp1d','ubp2d','ubp3d')
        grid
        figure
        plot(T,Y2d(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('y2d (N)')
        title('Saída desejada - Y2d')
        grid
        figure
        %plot(T,Y2d(:,1),'k',T,Y2d(:,2),'r',T,Y2d(:,3),'b',T,Y2d(:,4),'b',T,Y2d(:,5),'b',T,Y2d(:,6),'b')
        plot(T,Y2d(:,1),'k',T,Y2d(:,2),'r',T,Y2d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y2d')
        title('Saída desejada e derivadas - Y2d')
        %legend('y2d','yp2d','ypp2d','yppp2d','ypppp2d')
        legend('y2d','yp2d','ypp2d')
        grid
    case 9
        %Newmark:
        ubTilU = ZZ(:,1:3);
        ubTilpU = ZZ(:,4:6);
        ubTilppU = ZZ(:,7:9);
        Fbu = ZZ(:,10);
        y2 = ZZ(:,1:6)*C2.'*kNM1;
        u1 = ubTilU(:,1);
        u2 = ubTilU(:,2);
        u3 = ubTilU(:,3);
        up1 = ubTilpU(:,1);
        up2 = ubTilpU(:,2);
        up3 = ubTilpU(:,3);
        upp1 = ubTilppU(:,1);
        upp2 = ubTilppU(:,2);
        upp3 = ubTilppU(:,3);
        %Entrada:
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        %Variável de saída:
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,u1*10^3,'k',T,u2*10^3,'r',T,u3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,up1*10^3,'k',T,up2*10^3,'r',T,up3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,upp1,'k',T,upp2,'r',T,upp3,'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        if MEANativo
            Tpd = T(nPD:nT,1);
            ubTilpd = ZZpd(:,1:3);
            ubTilpPD = ZZpd(:,4:6);
            ubTilppPD = ZZpd(:,7:9);
            FbuPD = ZZpd(:,10);
            y2PD = ZZpd*C2.'*kNM1;
            ub1PD = ubTilpd(:,1);
            ub2PD = ubTilpd(:,2);
            ub3PD = ubTilpd(:,3);
            ubp1PD = ubTilpPD(:,1);
            ubp2PD = ubTilpPD(:,2);
            ubp3PD = ubTilpPD(:,3);
            ubpp1PD = ubTilppPD(:,1);
            ubpp2PD = ubTilppPD(:,2);
            ubpp3PD = ubTilppPD(:,3);
            %Entrada:
            figure
            plot(Tpd,FbuPD,'k')
            xlabel('Tempo (s)')
            ylabel('Fbu (N)')
            title('Força de atuação longitudinal [PD]')
            grid
            %Variável de saída:
            figure
            plot(Tpd,y2PD,'k')
            xlabel('Tempo (s)')
            ylabel('y2 (N)')
            title('Variável de saída: y2 = Nr(ub3) [PD]')
            grid
            %Posições longitudinais:
            figure
            plot(Tpd,ub1PD*10^3,'k',Tpd,ub2PD*10^3,'r',Tpd,ub3PD*10^3,'b')
            xlabel('Tempo (s)')
            ylabel('u (mm)')
            title('Posições longitudinais (u) [PD]')
            legend('ub1','ub2','ub3')
            grid
            %Velocidades longitudinais:
            figure
            plot(Tpd,ubp1PD*10^3,'k',Tpd,ubp2PD*10^3,'r',Tpd,ubp3PD*10^3,'b')
            xlabel('Tempo (s)')
            ylabel('up (mm/s)')
            title('Velocidades longitudinais (up) [PD]')
            legend('up1','up2','up3')
            grid
            %Acelerações longitudinais:
            figure
            plot(Tpd,ubpp1PD,'k',Tpd,ubpp2PD,'r',Tpd,ubpp3PD,'b')
            xlabel('Tempo (s)')
            ylabel('upp (mm/s2)')
            title('Acelerações longitudinais (upp) [PD]')
            legend('upp1','upp2','upp3')
            grid
        end
        %%
        %ODE45:
        if ODEativo
            ubTilU = ZZode(:,1:3);
            ubTilpU = ZZode(:,4:6);
            y2ode = ZZode*C2.'*kNM1;
            ub1od = ubTilU(:,1);
            ub2od = ubTilU(:,2);
            ub3od = ubTilU(:,3);
            ubp1od = ubTilpU(:,1);
            ubp2od = ubTilpU(:,2);
            ubp3od = ubTilpU(:,3);
            %Variável de saída:
            figure
            plot(T,y2ode,'k')
            xlabel('Tempo (s)')
            ylabel('y2 (N)')
            title('Variável de saída: y2 = Nr(ub3) [ODE]')
            grid
            %Posições longitudinais:
            figure
            plot(T,ub1od*10^3,'k',T,ub2od*10^3,'r',T,ub3od*10^3,'b')
            xlabel('Tempo (s)')
            ylabel('ub (mm)')
            title('Posições lineares (ub) [ODE]')
            legend('ub1','ub2','ub3')
            grid
            %Velocidades angulares:
            figure
            plot(T,ubp1od*10^3,'k',T,ubp2od*10^3,'r',T,ubp3od*10^3,'b')
            xlabel('Tempo (s)')
            ylabel('up (mm/s)')
            title('Velocidades longitudinais (up) [ODE]')
            legend('up1','up2','up3')
            grid
        end
    case 10
        nu = length(ubTilU0);
        nf = length(Fbu0);
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        ns = length(s10);
        
        %Newmark:
        ubTilU = ZZ(:,1:nu);
        ubTilpU = ZZ(:,nu+1:2*nu);
        ubTilppU = ZZ(:,2*nu+1:3*nu);
        Fbu = ZZ(:,3*nu+1:3*nu+nf);
        xTil = ZZ(:,3*nu+nf+1:3*nu+nf+nx);
        xTilp = ZZ(:,3*nu+nf+nx+1:3*nu+nf+2*nx);
        xTilpp = ZZ(:,3*nu+nf+2*nx+1:3*nu+nf+3*nx);
        Y1d = ZZ(:,3*nu+nf+3*nx+1:3*nu+nf+3*nx+ny);
        Zd = ZZ(:,3*nu+nf+3*nx+ny+1:3*nu+nf+3*nx+ny+nz);
        Uc = ZZ(:,3*nu+nf+3*nx+ny+nz+1:3*nu+nf+3*nx+ny+nz+nt);
        S1 = ZZ(:,3*nu+nf+3*nx+ny+nz+nt+1:3*nu+nf+3*nx+ny+nz+nt+ns);
        
        y2 = ZZ(:,1:2*nu)*C2.'*kNM1;
        u1 = ubTilU(:,1);
        u2 = ubTilU(:,2);
        u3 = ubTilU(:,3);
        up1 = ubTilpU(:,1);
        up2 = ubTilpU(:,2);
        up3 = ubTilpU(:,3);
        upp1 = ubTilppU(:,1);
        upp2 = ubTilppU(:,2);
        upp3 = ubTilppU(:,3);
        
        UtilT = xTil(:,1:nu);
        Mid = xTil(:,nu+1);
        FiBd = xTil(:,nu+2);
        UtilpT = xTilp(:,1:nu);
        Mipd = xTilp(:,nu+1);
        FipBd = xTilp(:,nu+2);
        UtilppT = xTilpp(:,1:nu);
        Mippd = xTilpp(:,nu+1);
        FippBd = xTilpp(:,nu+2);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        %%
        %Entrada:
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        %Variável de saída:
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,u1*10^3,'k',T,u2*10^3,'r',T,u3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,up1*10^3,'k',T,up2*10^3,'r',T,up3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,upp1,'k',T,upp2,'r',T,upp3,'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Vetor deslocamento e derivadas:
        figure
        plot(T,UtilT(:,1),'k',T,UtilT(:,2),'r',T,UtilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        figure
        plot(T,UtilpT(:,1)*30/pi,'k',T,UtilpT(:,2)*30/pi,'r',T,UtilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        figure
        plot(T,UtilppT(:,1),'k',T,UtilppT(:,2),'r',T,UtilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Variável de controle:
        figure
        plot(T,Uc,'k')
        xlabel('Tempo (s)')
        ylabel('Uc (Nm)')
        title('Variável de controle - Tx')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %%
        figure
        plot(Y1d(:,1),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
    case 11
        nu = length(ubTilU0);
        nf = length(Fbu0);
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        ns = length(s10);
        
        %Newmark:
        ubTilU = ZZ(:,1:nu);
        ubTilpU = ZZ(:,nu+1:2*nu);
        ubTilppU = ZZ(:,2*nu+1:3*nu);
        Fbu = ZZ(:,3*nu+1:3*nu+nf);
        xTil = ZZ(:,3*nu+nf+1:3*nu+nf+nx);
        xTilp = ZZ(:,3*nu+nf+nx+1:3*nu+nf+2*nx);
        xTilpp = ZZ(:,3*nu+nf+2*nx+1:3*nu+nf+3*nx);
        Y1d = ZZ(:,3*nu+nf+3*nx+1:3*nu+nf+3*nx+ny);
        Zd = ZZ(:,3*nu+nf+3*nx+ny+1:3*nu+nf+3*nx+ny+nz);
        Uc = ZZ(:,3*nu+nf+3*nx+ny+nz+1:3*nu+nf+3*nx+ny+nz+nt);
        S1 = ZZ(:,3*nu+nf+3*nx+ny+nz+nt+1:3*nu+nf+3*nx+ny+nz+nt+ns);
        
        u1 = ubTilU(:,1);
        u2 = ubTilU(:,2);
        u3 = ubTilU(:,3);
        up1 = ubTilpU(:,1);
        up2 = ubTilpU(:,2);
        up3 = ubTilpU(:,3);
        upp1 = ubTilppU(:,1);
        upp2 = ubTilppU(:,2);
        upp3 = ubTilppU(:,3);
        
        UtilT = xTil(:,1:nu);
        Mid = xTil(:,nu+1);
        FiBd = xTil(:,nu+2);
        UtilpT = xTilp(:,1:nu);
        Mipd = xTilp(:,nu+1);
        FipBd = xTilp(:,nu+2);
        UtilppT = xTilpp(:,1:nu);
        Mippd = xTilpp(:,nu+1);
        FippBd = xTilpp(:,nu+2);
        
        y2 = fNub3(u3);
        
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        %%
        %Entrada:
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        %Variável de saída:
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,u1*10^3,'k',T,u2*10^3,'r',T,u3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,up1*10^3,'k',T,up2*10^3,'r',T,up3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,upp1,'k',T,upp2,'r',T,upp3,'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Vetor deslocamento e derivadas:
        figure
        plot(T,UtilT(:,1),'k',T,UtilT(:,2),'r',T,UtilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        figure
        plot(T,UtilpT(:,1)*30/pi,'k',T,UtilpT(:,2)*30/pi,'r',T,UtilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        figure
        plot(T,UtilppT(:,1),'k',T,UtilppT(:,2),'r',T,UtilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Variável de controle:
        figure
        plot(T,Uc,'k')
        xlabel('Tempo (s)')
        ylabel('Uc (Nm)')
        title('Variável de controle - Tx')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
    case 12
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        nf = length(Fbu0);
        ns = length(s10);
        nu = length(ubTilAc0);
        
        %Newmark:
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        Y1d = ZZ(:,3*nx+1:3*nx+ny);
        Zd = ZZ(:,3*nx+ny+1:3*nx+ny+nz);
        TxUc = ZZ(:,3*nx+ny+nz+1);
        FbuUc = ZZ(:,3*nx+ny+nz+2);
        S1 = ZZ(:,3*nx+ny+nz+2+1:3*nx+ny+nz+2+ns);
        
        UbTilAc = xTil*AuAcx;
        UbpTilAc = xTilp*AuAcx;
        UbppTilAc = xTilpp*AuAcx;
        Mid = xTil*Amix;
        Mipd = xTilp*Amix;
        Mippd = xTilpp*Amix;
        FiBd = xTil*Afix;
        FipBd = xTilp*Afix;
        FippBd = xTilpp*Afix;
        
        uTilT = UbTilAc*Atu;
        uTilpT = UbpTilAc*Atu;
        uTilppT = UbppTilAc*Atu;
        ubTilU = UbTilAc*Auu;
        ubpTilU = UbpTilAc*Auu;
        ubppTilU = UbppTilAc*Auu;
        vL = UbTilAc*Avu;
        wL = UbTilAc*Awu;
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTilU(:,3);
        y2 = fNub3(ub3);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,FbuUc,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        figure
        plot(T,TxUc,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
    case 13
        %ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, Fbu0, s10];
        %xTil0 = [ubTilAc0; mid0; fiBdT0];
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        nf = length(Fbu0);
        ns = length(s10);
        nu = length(ubTilAc0);
        nuU = length(ubTilU0);
        
        %Newmark:
        ubTilUvt = ZZ(:,1:nuU);
        ubTilpUvt = ZZ(:,nuU+1:2*nuU);
        ubTilppUvt = ZZ(:,2*nuU+1:3*nuU);
        xTil = ZZ(:,3*nuU+1:3*nuU+nx);
        xTilp = ZZ(:,3*nuU+nx+1:3*nuU+2*nx);
        xTilpp = ZZ(:,3*nuU+2*nx+1:3*nuU+3*nx);
        Y1d = ZZ(:,3*nuU+3*nx+1:3*nuU+3*nx+ny);
        Zd = ZZ(:,3*nuU+3*nx+ny+1:3*nuU+3*nx+ny+nz);
        TxUc = ZZ(:,3*nuU+3*nx+ny+nz+1);
        FbuUc = ZZ(:,3*nuU+3*nx+ny+nz+2);
        S1 = ZZ(:,3*nuU+3*nx+ny+nz+2+1:3*nuU+3*nx+ny+nz+2+ns);
        
        UbTilAc = xTil*AuAcx;
        UbpTilAc = xTilp*AuAcx;
        UbppTilAc = xTilpp*AuAcx;
        Mid = xTil*Amix;
        Mipd = xTilp*Amix;
        Mippd = xTilpp*Amix;
        FiBd = xTil*Afix;
        FipBd = xTilp*Afix;
        FippBd = xTilpp*Afix;
        
        uTilT = UbTilAc*Atu;
        uTilpT = UbpTilAc*Atu;
        uTilppT = UbppTilAc*Atu;
        ubTilU = UbTilAc*Auu;
        ubpTilU = UbpTilAc*Auu;
        ubppTilU = UbppTilAc*Auu;
        vL = UbTilAc*Avu;
        wL = UbTilAc*Awu;
        uL = (vL.^2 + wL.^2).^(1/2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTilU(:,3);
        y2 = fNub3(ub3);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,FbuUc,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        figure
        plot(T,TxUc,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilUvt(:,1)*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub1 (mm)')
        title('Posições longitudinais (ub1): planta e controlador')
        legend('planta','controlador')
        grid
        figure
        plot(T,ubTilU(:,2)*10^3,'k',T,ubTilUvt(:,2)*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub2 (mm)')
        title('Posições longitudinais (ub2: planta e controlador')
        legend('planta','controlador')
        grid
        figure
        plot(T,ubTilU(:,3)*10^3,'k',T,ubTilUvt(:,3)*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub3 (mm)')
        title('Posições longitudinais (ub3): planta e controlador')
        legend('planta','controlador')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
    case 14
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        nf = length(Fbu0);
        ns = length(s10);
        nu = length(ubTilAc0);
        
        %Newmark:
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        Y1d = ZZ(:,3*nx+1:3*nx+ny);
        Zd = ZZ(:,3*nx+ny+1:3*nx+ny+nz);
        TxUc = ZZ(:,3*nx+ny+nz+1);
        FbuUc = ZZ(:,3*nx+ny+nz+2);
        S1 = ZZ(:,3*nx+ny+nz+2+1:3*nx+ny+nz+2+ns);
        
        UbTilAc = xTil*AuAcx;
        UbpTilAc = xTilp*AuAcx;
        UbppTilAc = xTilpp*AuAcx;
        Mid = xTil*Amix;
        Mipd = xTilp*Amix;
        Mippd = xTilpp*Amix;
        FiBd = xTil*Afix;
        FipBd = xTilp*Afix;
        FippBd = xTilpp*Afix;
        
        uTilT = UbTilAc*Atu;
        uTilpT = UbpTilAc*Atu;
        uTilppT = UbppTilAc*Atu;
        ubTilU = UbTilAc*Auu;
        ubpTilU = UbpTilAc*Auu;
        ubppTilU = UbppTilAc*Auu;
        vL = UbTilAc*Avu;
        wL = UbTilAc*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTilU(:,3);
        y2 = fNub3(ub3);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,FbuUc,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        figure
        plot(T,TxUc,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
    case 15
        %ZZ(1,:) = [ubTilAc0.', ubpTilAc0.', ubppTilAc0.', Tx0, Fbu0];
        nu = length(ubTilAc0);
        nt = length(Tx0);
        nf = length(Fbu0);
        %Newmark:
        ubTil = ZZ(:,1:nu);
        ubTilp = ZZ(:,nu+1:2*nu);
        ubTilpp = ZZ(:,2*nu+1:3*nu);
        Tx = ZZ(:,3*nu+nt);
        Fbu = ZZ(:,3*nu+nt+nf);
        
        uTilT = ubTil*Atu;
        uTilpT = ubTilp*Atu;
        uTilppT = ubTilpp*Atu;
        ubTilU = ubTil*Auu;
        ubTilpU = ubTilp*Auu;
        ubTilppU = ubTilpp*Auu;
        vL = ubTil*Avu;
        wL = ubTil*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        teta1 = uTilT(:,1);
        teta2 = uTilT(:,2);
        teta3 = uTilT(:,3);
        tetap1 = uTilpT(:,1);
        tetap2 = uTilpT(:,2);
        tetap3 = uTilpT(:,3);
        tetapp1 = uTilppT(:,1);
        tetapp2 = uTilppT(:,2);
        tetapp3 = uTilppT(:,3);
        
        ub1 = ubTilU(:,1);
        ub2 = ubTilU(:,2);
        ub3 = ubTilU(:,3);
        ubp1 = ubTilpU(:,1);
        ubp2 = ubTilpU(:,2);
        ubp3 = ubTilpU(:,3);
        ubpp1 = ubTilppU(:,1);
        ubpp2 = ubTilppU(:,2);
        ubpp3 = ubTilppU(:,3);
        
        y1 = tetap3;
        y2 = fNub3(ub3);
        
        %Entrada:
        figure
        plot(T,Tx,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de entrada')
        grid
        %Saídas:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('teta3p (RPM)')
        title('Velocidade teta3p')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Força normal no contato - Nr(ub3)')
        grid
        %Posições angulares:
        figure
        plot(T,teta1,'k',T,teta2,'r',T,teta3,'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Posições angulares (teta)')
        legend('teta1','teta2','teta3')
        grid
        %Velocidades angulares:
        figure
        plot(T,tetap1*(30/pi),'k',T,tetap2*(30/pi),'r',T,tetap3*(30/pi),'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidades angulares (tetap)')
        legend('tetap1','tetap2','tetap3')
        grid
        %Acelerações angulares:
        figure
        plot(T,tetapp1,'k',T,tetapp2,'r',T,tetapp3,'b')
        xlabel('Tempo (s)')
        ylabel('tetapp (rad/s2)')
        title('Acelerações angulares (tetapp)')
        legend('tetapp1','tetapp2','tetapp3')
        grid
        %Posições longitudinais:
        figure
        plot(T,ub1*10^3,'k',T,ub2*10^3,'r',T,ub3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        title('Posições longitudinais (ub)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubp1*10^3,'k',T,ubp2*10^3,'r',T,ubp3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubp (mm/s)')
        title('Velocidades longitudinais (ubp)')
        legend('ubp1','ubp2','ubp3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubpp1*10^3,'k',T,ubpp2*10^3,'r',T,ubpp3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        title('Acelerações longitudinais (ubpp)')
        legend('ubpp1','ubpp2','ubpp3')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
    case 16
        %ZZ(1,:) = [ubTilAc0.', ubpTilAc0.', ubppTilAc0.', Tx0, Fbu0];
        nu = length(ubTilAc0);
        nt = length(Tx0);
        nf = length(Fbu0);
        %Newmark:
        ubTil = ZZ(:,1:nu);
        ubTilp = ZZ(:,nu+1:2*nu);
        ubTilpp = ZZ(:,2*nu+1:3*nu);
        Tx = ZZ(:,3*nu+nt);
        Fbu = ZZ(:,3*nu+nt+nf);
        
        uTilT = ubTil*Atu;
        uTilpT = ubTilp*Atu;
        uTilppT = ubTilpp*Atu;
        ubTilU = ubTil*Auu;
        ubTilpU = ubTilp*Auu;
        ubTilppU = ubTilpp*Auu;
        vL = ubTil*Avu;
        wL = ubTil*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        teta1 = uTilT(:,1);
        teta2 = uTilT(:,2);
        teta3 = uTilT(:,3);
        tetap1 = uTilpT(:,1);
        tetap2 = uTilpT(:,2);
        tetap3 = uTilpT(:,3);
        tetapp1 = uTilppT(:,1);
        tetapp2 = uTilppT(:,2);
        tetapp3 = uTilppT(:,3);
        
        ub1 = ubTilU(:,1);
        ub2 = ubTilU(:,2);
        ub3 = ubTilU(:,3);
        ubp1 = ubTilpU(:,1);
        ubp2 = ubTilpU(:,2);
        ubp3 = ubTilpU(:,3);
        ubpp1 = ubTilppU(:,1);
        ubpp2 = ubTilppU(:,2);
        ubpp3 = ubTilppU(:,3);
        
        y1 = tetap3;
        y2 = fNub3(ub3);
        
        %Entrada:
        figure
        plot(T,Tx,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de entrada')
        grid
        %Saídas:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('teta3p (RPM)')
        title('Velocidade teta3p')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Força normal no contato - Nr(ub3)')
        grid
        %Posições angulares:
        figure
        plot(T,teta1,'k',T,teta2,'r',T,teta3,'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Posições angulares (teta)')
        legend('teta1','teta2','teta3')
        grid
        %Velocidades angulares:
        figure
        plot(T,tetap1*(30/pi),'k',T,tetap2*(30/pi),'r',T,tetap3*(30/pi),'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidades angulares (tetap)')
        legend('tetap1','tetap2','tetap3')
        grid
        %Acelerações angulares:
        figure
        plot(T,tetapp1,'k',T,tetapp2,'r',T,tetapp3,'b')
        xlabel('Tempo (s)')
        ylabel('tetapp (rad/s2)')
        title('Acelerações angulares (tetapp)')
        legend('tetapp1','tetapp2','tetapp3')
        grid
        %Posições longitudinais:
        figure
        plot(T,ub1*10^3,'k',T,ub2*10^3,'r',T,ub3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        title('Posições longitudinais (ub)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubp1*10^3,'k',T,ubp2*10^3,'r',T,ubp3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubp (mm/s)')
        title('Velocidades longitudinais (ubp)')
        legend('ubp1','ubp2','ubp3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubpp1*10^3,'k',T,ubpp2*10^3,'r',T,ubpp3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        title('Acelerações longitudinais (ubpp)')
        legend('ubpp1','ubpp2','ubpp3')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
    case 17
    case 18
        %Newmark:
        ubTilU = ZZ(:,1:3);
        ubTilpU = ZZ(:,4:6);
        ubTilppU = ZZ(:,7:9);
        Fbu = ZZ(:,10);
        u1 = ubTilU(:,1);
        u2 = ubTilU(:,2);
        u3 = ubTilU(:,3);
        up1 = ubTilpU(:,1);
        up2 = ubTilpU(:,2);
        up3 = ubTilpU(:,3);
        upp1 = ubTilppU(:,1);
        upp2 = ubTilppU(:,2);
        upp3 = ubTilppU(:,3);
        y2 = fNub3(u3);
        %Entrada:
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        %Variável de saída:
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,u1*10^3,'k',T,u2*10^3,'r',T,u3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,up1*10^3,'k',T,up2*10^3,'r',T,up3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,upp1,'k',T,upp2,'r',T,upp3,'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        
        %Newmark:
        uTilT = ZZ(:,10+1:10+3);
        uTilpT = ZZ(:,10+4:10+6);
        uTilppT = ZZ(:,10+7:10+9);
        Tx = ZZ(:,10+10);
        y1 = ZZ(:,10+1:10+6)*C1.';
        teta1 = uTilT(:,1);
        teta2 = uTilT(:,2);
        teta3 = uTilT(:,3);
        tetap1 = uTilpT(:,1);
        tetap2 = uTilpT(:,2);
        tetap3 = uTilpT(:,3);
        tetapp1 = uTilppT(:,1);
        tetapp2 = uTilppT(:,2);
        tetapp3 = uTilppT(:,3);
        Tatr = fTatc(tetap3,y2);
        %Entrada:
        figure
        plot(T,Tx,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variável de saída:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = tetap3')
        grid
        %Posições angulares:
        figure
        plot(T,teta1,'k',T,teta2,'r',T,teta3,'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Posições angulares (teta)')
        legend('teta1','teta2','teta3')
        grid
        %Velocidades angulares:
        figure
        plot(T,tetap1*(30/pi),'k',T,tetap2*(30/pi),'r',T,tetap3*(30/pi),'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidades angulares (tetap)')
        legend('tetap1','tetap2','tetap3')
        grid
        %Acelerações angulares:
        figure
        plot(T,tetapp1,'k',T,tetapp2,'r',T,tetapp3,'b')
        xlabel('Tempo (s)')
        ylabel('tetapp (rad/s2)')
        title('Acelerações angulares (tetapp)')
        legend('tetapp1','tetapp2','tetapp3')
        grid
        %Torque de atrito:
        figure
        plot(T,Tx,'k',T,Tatr,'r')
        xlabel('Tempo (s)')
        ylabel('Tx, Tatr (Nm)')
        title('Torques de entrada e de atrito')
        grid
    case 19
        %ZZ(i,:) = [ubTilrnM1.', ubTilprnM1.', ubTilpprnM1.', TxnM1, FbunM1];
        ubr = ZZ(:,1:nIdNgr);
        ubpr = ZZ(:,nIdNgr+1:2*nIdNgr);
        ubppr = ZZ(:,2*nIdNgr+1:3*nIdNgr);
        Tx = ZZ(:,3*nIdNgr+1);
        Fbu = ZZ(:,3*nIdNgr+2);
        
        ub = ubr*SigmaReF.';
        ubp = ubpr*SigmaReF.';
        ubpp = ubppr*SigmaReF.';
        
        ubU = ub*IdU.';
        ubV = ub*IdV.';
        ubTz = ub*IdTetaZ.';
        ubW = ub*IdW.';
        ubTy = ub*IdTetaY.';
        ubTx = ub*IdTetaX.';
        ubpU = ubp*IdU.';
        ubpV = ubp*IdV.';
        ubpTz = ubp*IdTetaZ.';
        ubpW = ubp*IdW.';
        ubpTy = ubp*IdTetaY.';
        ubpTx = ubp*IdTetaX.';
        ubppU = ubpp*IdU.';
        ubppV = ubpp*IdV.';
        ubppTz = ubpp*IdTetaZ.';
        ubppW = ubpp*IdW.';
        ubppTy = ubpp*IdTetaY.';
        ubppTx = ubpp*IdTetaX.';
        
        y1 = ubpTx(:,N+1);
        ubNM1 = ubU(:,N+1);
        y2 = FnbNM1(ubNM1);
        vL = ubV(:,N1+1);
        wL = ubW(:,N1+1);
        uL = sqrt(vL.^2 + wL.^2);
        
        %Entradas:
        figure
        plot(T,Tx,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de entrada')
        grid
        %Saídas:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Velocidade TetaXp do rotor principal')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Força normal de contato')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
    case 20
        %ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, Fbu0, s10];
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        nf = length(Fbu0);
        ns = length(s10);
        
        %Newmark:
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        Y1d = ZZ(:,3*nx+1:3*nx+ny);
        Zd = ZZ(:,3*nx+ny+1:3*nx+ny+nz);
        TxUc = ZZ(:,3*nx+ny+nz+1);
        FbuUc = ZZ(:,3*nx+ny+nz+2);
        S1 = ZZ(:,3*nx+ny+nz+2+1:3*nx+ny+nz+2+ns);
        
        UbTilAc = xTil*AuAcx;
        UbpTilAc = xTilp*AuAcx;
        UbppTilAc = xTilpp*AuAcx;
        Mid = xTil*Amix;
        Mipd = xTilp*Amix;
        Mippd = xTilpp*Amix;
        FiBd = xTil*Afix;
        FipBd = xTilp*Afix;
        FippBd = xTilpp*Afix;
        
        uTilT = UbTilAc*Atu;
        uTilpT = UbpTilAc*Atu;
        uTilppT = UbppTilAc*Atu;
        ubTilU = UbTilAc*Auu;
        ubpTilU = UbpTilAc*Auu;
        ubppTilU = UbppTilAc*Auu;
        vL = UbTilAc*Avu;
        wL = UbTilAc*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTilU(:,3);
        y2 = fNub3(ub3);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,FbuUc,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        figure
        plot(T,TxUc,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
    case 21
        %ZZ(i,:) = [ubTilnM1.', ubTilpnM1.', ubTilppnM1.', lmbnM1.', TxnM1, FbunM1];
        ub = ZZ(:,1:Ng);
        ubp = ZZ(:,Ng+1:2*Ng);
        ubpp = ZZ(:,2*Ng+1:3*Ng);
        lmb = kNmk*ZZ(:,3*Ng+1:3*Ng+rEstr);
        Tx = ZZ(:,3*Ng+rEstr+1);
        Fbu = ZZ(:,3*Ng+rEstr+2);
        
        ubU = ub*IdU.';
        ubV = ub*IdV.';
        ubTz = ub*IdTetaZ.';
        ubW = ub*IdW.';
        ubTy = ub*IdTetaY.';
        ubTx = ub*IdTetaX.';
        ubpU = ubp*IdU.';
        ubpV = ubp*IdV.';
        ubpTz = ubp*IdTetaZ.';
        ubpW = ubp*IdW.';
        ubpTy = ubp*IdTetaY.';
        ubpTx = ubp*IdTetaX.';
        ubppU = ubpp*IdU.';
        ubppV = ubpp*IdV.';
        ubppTz = ubpp*IdTetaZ.';
        ubppW = ubpp*IdW.';
        ubppTy = ubpp*IdTetaY.';
        ubppTx = ubpp*IdTetaX.';
        
        y1 = ubpTx(:,N+1);
        ubNM1 = ubU(:,N+1);
        y2 = FnbNM1(ubNM1);
        vL = ubV(:,N1+1);
        wL = ubW(:,N1+1);
        uL = sqrt(vL.^2 + wL.^2);
        
        %Entradas:
        figure
        plot(T,Tx,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de entrada')
        grid
        %Saídas:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Velocidade TetaXp do rotor principal')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Força normal de contato')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,lmb(:,1),'k',T,lmb(:,3),'r',T,lmb(:,5),'b',T,lmb(:,7),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,lmb(:,2),'k',T,lmb(:,4),'r',T,lmb(:,6),'b',T,lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 22
        %ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', lmb0.', y1dTil0.', Zd0.', Tx0, Fbu0, s10];
        nx = length(xTil0);
        nLmb = length(lmb0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        nf = length(Fbu0);
        ns = length(s10);
        
        %Newmark:
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        lmb = kNmk*ZZ(:,3*nx+1:3*nx+nLmb);
        Y1d = ZZ(:,3*nx+nLmb+1:3*nx+nLmb+ny);
        Zd = ZZ(:,3*nx+nLmb+ny+1:3*nx+nLmb+ny+nz);
        TxUc = ZZ(:,3*nx+nLmb+ny+nz+1);
        FbuUc = ZZ(:,3*nx+nLmb+ny+nz+2);
        S1 = ZZ(:,3*nx+nLmb+ny+nz+2+1:3*nx+nLmb+ny+nz+2+ns);
        
        UbTil = xTil*Aux;
        UbpTil = xTilp*Aux;
        UbppTil = xTilpp*Aux;
        Mid = xTil*Amix;
        Mipd = xTilp*Amix;
        Mippd = xTilpp*Amix;
        FiBd = xTil*Afix;
        FipBd = xTilp*Afix;
        FippBd = xTilpp*Afix;
        
        uTilT = UbTil*Atu;
        uTilpT = UbpTil*Atu;
        uTilppT = UbppTil*Atu;
        ubTilU = UbTil*Auu;
        ubpTilU = UbpTil*Auu;
        ubppTilU = UbppTil*Auu;
        vL = UbTil*Avu;
        wL = UbTil*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTilU(:,3);
        y2 = fNub3(ub3);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,FbuUc,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        figure
        plot(T,TxUc,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,lmb(:,1),'k',T,lmb(:,3),'r',T,lmb(:,5),'b',T,lmb(:,7),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,lmb(:,2),'k',T,lmb(:,4),'r',T,lmb(:,6),'b',T,lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 23
        %ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', lmb0.', y1dTil0.', Zd0.', Tx0, s10];
        nx = length(xTil0);
        nLmb = length(lmb0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        ns = length(s10);
        
        %Newmark:
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        lmb = kNmk*ZZ(:,3*nx+1:3*nx+nLmb);
        Y1d = ZZ(:,3*nx+nLmb+1:3*nx+nLmb+ny);
        Zd = ZZ(:,3*nx+nLmb+ny+1:3*nx+nLmb+ny+nz);
        TxUc = ZZ(:,3*nx+nLmb+ny+nz+1);
        S1 = ZZ(:,3*nx+nLmb+ny+nz+1+1:3*nx+nLmb+ny+nz+1+ns);
        
        UbTil = xTil*Aux;
        UbpTil = xTilp*Aux;
        UbppTil = xTilpp*Aux;
        Mid = xTil*Amix;
        Mipd = xTilp*Amix;
        Mippd = xTilpp*Amix;
        FiBd = xTil*Afix;
        FipBd = xTilp*Afix;
        FippBd = xTilpp*Afix;
        
        uTilT = UbTil*Atu;
        uTilpT = UbpTil*Atu;
        uTilppT = UbppTil*Atu;
        ubTilU = UbTil*Auu;
        ubpTilU = UbpTil*Auu;
        ubppTilU = UbppTil*Auu;
        vL = UbTil*Avu;
        wL = UbTil*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTilU(:,3);
        y2mm = fNub3(ub3);
        y2 = FnbNM1(ub3);
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbUb1 (N)')
        title('Força de reação longitudinal')
        grid
        figure
        plot(T,TxUc,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k',T,y2mm,'r')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        legend('Planta','Metamodelo')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,lmb(:,2),'k',T,lmb(:,4),'r',T,lmb(:,6),'b',T,lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,lmb(:,3),'k',T,lmb(:,5),'r',T,lmb(:,7),'b',T,lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 24
        nu = length(ubTilU0);
        nf = length(Fbu0);
        nx = length(xTil0);
        ny = length(y1dTil0);
        nz = length(Zd0);
        nt = length(Tx0);
        ns = length(s10);
        
        %Newmark:
        ubTilU = ZZ(:,1:nu);
        ubTilpU = ZZ(:,nu+1:2*nu);
        ubTilppU = ZZ(:,2*nu+1:3*nu);
        Fbu = ZZ(:,3*nu+1:3*nu+nf);
        xTil = ZZ(:,3*nu+nf+1:3*nu+nf+nx);
        xTilp = ZZ(:,3*nu+nf+nx+1:3*nu+nf+2*nx);
        xTilpp = ZZ(:,3*nu+nf+2*nx+1:3*nu+nf+3*nx);
        Y1d = ZZ(:,3*nu+nf+3*nx+1:3*nu+nf+3*nx+ny);
        Zd = ZZ(:,3*nu+nf+3*nx+ny+1:3*nu+nf+3*nx+ny+nz);
        Uc = ZZ(:,3*nu+nf+3*nx+ny+nz+1:3*nu+nf+3*nx+ny+nz+nt);
        S1 = ZZ(:,3*nu+nf+3*nx+ny+nz+nt+1:3*nu+nf+3*nx+ny+nz+nt+ns);
        
        u1 = ubTilU(:,1);
        u2 = ubTilU(:,2);
        u3 = ubTilU(:,3);
        up1 = ubTilpU(:,1);
        up2 = ubTilpU(:,2);
        up3 = ubTilpU(:,3);
        upp1 = ubTilppU(:,1);
        upp2 = ubTilppU(:,2);
        upp3 = ubTilppU(:,3);
        
        UtilT = xTil(:,1:nu);
        Mid = xTil(:,nu+1);
        FiBd = xTil(:,nu+2);
        UtilpT = xTilp(:,1:nu);
        Mipd = xTilp(:,nu+1);
        FipBd = xTilp(:,nu+2);
        UtilppT = xTilpp(:,1:nu);
        Mippd = xTilpp(:,nu+1);
        FippBd = xTilpp(:,nu+2);
        
        y1 = [UtilT, UtilpT]*C1.';
        y2 = fNub3(u3);
        
        teta1d = Zd(:,1);
        teta2d = Zd(:,2);
        teta3d = Zd(:,3);
        tetap1d = Zd(:,4);
        tetap2d = Zd(:,5);
        tetap3d = Zd(:,6);
        %%
        %Entrada:
        figure
        plot(T,Fbu,'k')
        xlabel('Tempo (s)')
        ylabel('Fbu (N)')
        title('Força de atuação longitudinal')
        grid
        %Variável de saída:
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2 = Nr(ub3)')
        grid
        %Posições longitudinais:
        figure
        plot(T,u1*10^3,'k',T,u2*10^3,'r',T,u3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,up1*10^3,'k',T,up2*10^3,'r',T,up3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,upp1,'k',T,upp2,'r',T,upp3,'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Vetor deslocamento e derivadas:
        figure
        plot(T,UtilT(:,1),'k',T,UtilT(:,2),'r',T,UtilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        figure
        plot(T,UtilpT(:,1)*30/pi,'k',T,UtilpT(:,2)*30/pi,'r',T,UtilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        figure
        plot(T,UtilppT(:,1),'k',T,UtilppT(:,2),'r',T,UtilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Variável de controle:
        figure
        plot(T,Uc,'k')
        xlabel('Tempo (s)')
        ylabel('Uc (Nm)')
        title('Variável de controle - Tx')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBd,'k',T,-FiBd,'k',T,S1,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBd,'k',T,FipBd,'r',T,FippBd,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,Mid,'k',T,Mipd,'r',T,Mippd,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória de desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,Y1d(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        figure
        %plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b',T,Y1d(:,4),'b',T,Y1d(:,5),'b')
        plot(T,Y1d(:,1),'k',T,Y1d(:,2),'r',T,Y1d(:,3),'b')
        xlabel('Tempo (s)')
        %ylabel('y1d')
        title('Saída desejada e derivadas - Y1d')
        %legend('y1d','yp1d','ypp1d','yppp1d','ypppp1d')
        legend('y1d','yp1d','ypp1d')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
    case 25
        %ZZ(i,:) = [ubTilnM1.', ubTilpnM1.', ubTilppnM1.', lmbnM1.'];
        ub = ZZ(:,1:Ng);
        ubp = ZZ(:,Ng+1:2*Ng);
        ubpp = ZZ(:,2*Ng+1:3*Ng);
        lmb = kNmk*ZZ(:,3*Ng+1:3*Ng+rEstrH);
        
        ubU = ub*IdU.';
        ubV = ub*IdV.';
        ubTz = ub*IdTetaZ.';
        ubW = ub*IdW.';
        ubTy = ub*IdTetaY.';
        ubTx = ub*IdTetaX.';
        ubpU = ubp*IdU.';
        ubpV = ubp*IdV.';
        ubpTz = ubp*IdTetaZ.';
        ubpW = ubp*IdW.';
        ubpTy = ubp*IdTetaY.';
        ubpTx = ubp*IdTetaX.';
        ubppU = ubpp*IdU.';
        ubppV = ubpp*IdV.';
        ubppTz = ubpp*IdTetaZ.';
        ubppW = ubpp*IdW.';
        ubppTy = ubpp*IdTetaY.';
        ubppTx = ubpp*IdTetaX.';
        
        y1 = ubpTx(:,N+1);
        ubNM1 = ubU(:,N+1);
        y2 = FnbNM1(ubNM1);
        vL = ubV(:,N1+1);
        wL = ubW(:,N1+1);
        uL = sqrt(vL.^2 + wL.^2);
        
        %Entradas:
        figure
        plot(T,lmb(:,2),'k')
        xlabel('Tempo (s)')
        ylabel('LmbT (Nm)')
        title('Torque de reação')
        grid
        figure
        plot(T,lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbU (N)')
        title('Força de reação')
        grid
        %Saídas:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Velocidade TetaXp do rotor principal')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Força normal de contato')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,lmb(:,3),'k',T,lmb(:,5),'r',T,lmb(:,7),'b',T,lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,lmb(:,4),'k',T,lmb(:,6),'r',T,lmb(:,8),'b',T,lmb(:,10),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 26
        %ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmb0));
        %ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmbnM1.'];
        nu = length(ubTilU0);
        nx = length(xTil0);
        nyd = length(y1dTil0);
        nZ = length(Zd0);
        nTx = length(Tx0);
        ns = length(s10);
        nlmb = length(lmb0);
        ubTil = ZZ(:,1:nu);
        ubpTil = ZZ(:,nu+1:2*nu);
        ubppTil = ZZ(:,2*nu+1:3*nu);
        xTil = ZZ(:,3*nu+1:3*nu+nx);
        xTilp = ZZ(:,3*nu+nx+1:3*nu+2*nx);
        xTilpp = ZZ(:,3*nu+2*nx+1:3*nu+3*nx);
        y1dVt = ZZ(:,3*nu+3*nx+1:3*nu+3*nx+nyd);
        ZdDInv = ZZ(:,3*nu+3*nx+nyd+1:3*nu+3*nx+nyd+nZ);
        TxIn = ZZ(:,3*nu+3*nx+nyd+nZ+1);
        s1Tx = ZZ(:,3*nu+3*nx+nyd+nZ+2);
        Lmb = kNmk*ZZ(:,3*nu+3*nx+nyd+nZ+2+1:3*nu+3*nx+nyd+nZ+2+nlmb);
        
        ubEF = xTil(:,1:Ng);
        midVc = xTil(:,Ng+1);
        FiBdVc = xTil(:,Ng+2);
        ubpEF = xTilp(:,1:Ng);
        mipdVc = xTilp(:,Ng+1);
        FipBdVc = xTilp(:,Ng+2);
        ubppEF = xTilpp(:,1:Ng);
        mippdVc = xTilpp(:,Ng+1);
        FippBdVc = xTilpp(:,Ng+2);
        
        uTilT = ubEF*Atu;
        uTilpT = ubpEF*Atu;
        uTilppT = ubppEF*Atu;
        ubTilU = ubEF*Auu;
        ubpTilU = ubpEF*Auu;
        ubppTilU = ubppEF*Auu;
        vL = ubEF*Avu;
        wL = ubEF*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTil(:,3);
        y2mm = fNub3(ub3);
        ubNM1 = ubEF(:,6*(N+1)-5);
        y2 = FnbNM1(ubNM1);
        teta1d = ZdDInv(:,1);
        teta2d = ZdDInv(:,2);
        teta3d = ZdDInv(:,3);
        tetap1d = ZdDInv(:,4);
        tetap2d = ZdDInv(:,5);
        tetap3d = ZdDInv(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,Lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbUb1 (N)')
        title('Força de reação longitudinal')
        grid
        figure
        plot(T,TxIn,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k',T,y2mm,'r')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2')
        legend('Planta','Metamodelo')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBdVc,'k',T,-FiBdVc,'k',T,s1Tx,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBdVc,'k',T,FipBdVc,'r',T,FippBdVc,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,midVc,'k',T,mipdVc,'r',T,mippdVc,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,Lmb(:,2),'k',T,Lmb(:,4),'r',T,Lmb(:,6),'b',T,Lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,Lmb(:,3),'k',T,Lmb(:,5),'r',T,Lmb(:,7),'b',T,Lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 27
        %ZZ(i,:) = [uTilTxnM1.', uTilpTxnM1.', uTilppTxnM1.', lmbnM1.',tX1,tpX1,tppX1];
        uTx = ZZ(:,1:Ng-1);
        upTx = ZZ(:,Ng-1+1:2*(Ng-1));
        uppTx = ZZ(:,2*(Ng-1)+1:3*(Ng-1));
        lmb = kNmk*ZZ(:,3*(Ng-1)+1:3*(Ng-1)+rEstrH);
        tX1 = ZZ(:,3*(Ng-1)+rEstrH+1);
        tpX1 = ZZ(:,3*(Ng-1)+rEstrH+2);
        tppX1 = ZZ(:,3*(Ng-1)+rEstrH+3);
        ub = uTx*Auru.' + tX1*AtXu.';
        ubp = upTx*Auru.' + tpX1*AtXu.';
        ubpp = uppTx*Auru.' + tppX1*AtXu.';
        
        ubU = ub*IdU.';
        ubV = ub*IdV.';
        ubTz = ub*IdTetaZ.';
        ubW = ub*IdW.';
        ubTy = ub*IdTetaY.';
        ubTx = ub*IdTetaX.';
        ubpU = ubp*IdU.';
        ubpV = ubp*IdV.';
        ubpTz = ubp*IdTetaZ.';
        ubpW = ubp*IdW.';
        ubpTy = ubp*IdTetaY.';
        ubpTx = ubp*IdTetaX.';
        ubppU = ubpp*IdU.';
        ubppV = ubpp*IdV.';
        ubppTz = ubpp*IdTetaZ.';
        ubppW = ubpp*IdW.';
        ubppTy = ubpp*IdTetaY.';
        ubppTx = ubpp*IdTetaX.';
        
        y1 = ubpTx(:,N+1);
        ubNM1 = ubU(:,N+1);
        y2 = FnbNM1(ubNM1);
        vL = ubV(:,N1+1);
        wL = ubW(:,N1+1);
        uL = sqrt(vL.^2 + wL.^2);
        
        %Entradas:
        figure
        plot(T,lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbU (N)')
        %title('Força de reação')
        grid
        %saveas(gcf,'figLmbU.fig')
        %Saídas:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        %title('Velocidade TetaXp do rotor principal')
        grid
        %saveas(gcf,'figy1.fig')
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        %title('Força normal de contato')
        grid
        %saveas(gcf,'figy2.fig')
        %Posições longitudinais:
        figure
        plot(T,ubU(:,1)*10^3,'k',T,ubU(:,N1+1)*10^3,'r',T,ubU(:,N+1)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        %saveas(gcf,'figvwt.fig')
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %saveas(gcf,'figvw.fig')
        %Reações referentes às restrições:
        figure
        plot(T,lmb(:,2),'k',T,lmb(:,4),'r',T,lmb(:,6),'b',T,lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        %title('Reações de força nas extremidades')
        grid
        %saveas(gcf,'figlmbN.fig')
        figure
        plot(T,lmb(:,3),'k',T,lmb(:,5),'r',T,lmb(:,7),'b',T,lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
        %saveas(gcf,'figlmbNm.fig')
    case 28
        %ZZ(i,:) = [uTilTxnM1.', uTilpTxnM1.', uTilppTxnM1.', lmbnM1.',tX1,tpX1,tppX1];
        uTx = ZZ(:,1:Ng-1);
        upTx = ZZ(:,Ng-1+1:2*(Ng-1));
        uppTx = ZZ(:,2*(Ng-1)+1:3*(Ng-1));
        lmb = kNmk*ZZ(:,3*(Ng-1)+1:3*(Ng-1)+rEstrH);
        tX1 = ZZ(:,3*(Ng-1)+rEstrH+1);
        tpX1 = ZZ(:,3*(Ng-1)+rEstrH+2);
        tppX1 = ZZ(:,3*(Ng-1)+rEstrH+3);
        ub = uTx*Auru.' + tX1*AtXu.';
        ubp = upTx*Auru.' + tpX1*AtXu.';
        ubpp = uppTx*Auru.' + tppX1*AtXu.';
        
        ubU = ub*IdU.';
        ubV = ub*IdV.';
        ubTz = ub*IdTetaZ.';
        ubW = ub*IdW.';
        ubTy = ub*IdTetaY.';
        ubTx = ub*IdTetaX.';
        ubpU = ubp*IdU.';
        ubpV = ubp*IdV.';
        ubpTz = ubp*IdTetaZ.';
        ubpW = ubp*IdW.';
        ubpTy = ubp*IdTetaY.';
        ubpTx = ubp*IdTetaX.';
        ubppU = ubpp*IdU.';
        ubppV = ubpp*IdV.';
        ubppTz = ubpp*IdTetaZ.';
        ubppW = ubpp*IdW.';
        ubppTy = ubpp*IdTetaY.';
        ubppTx = ubpp*IdTetaX.';
        
        y1 = ubpTx(:,N+1);
        ubNM1 = ubU(:,N+1);
        y2 = FnbNM1(ubNM1);
        vL = ubV(:,N1+1);
        wL = ubW(:,N1+1);
        uL = sqrt(vL.^2 + wL.^2);
        
        %Entradas:
        figure
        plot(T,lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbU (N)')
        title('Força de reação')
        grid
        %Saídas:
        figure
        plot(T,y1*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Velocidade TetaXp do rotor principal')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Força normal de contato')
        grid
        %Velocidades angulares:
        figure
        plot(T,ubpTx(:,1)*30/pi,'k',T,ubpTx(:,N+1)*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        legend('servo-motor','rotor principal')
        grid
        %Posições longitudinais:
        figure
        plot(T,u1*10^3,'k',T,u2*10^3,'r',T,u3*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,lmb(:,2),'k',T,lmb(:,4),'r',T,lmb(:,6),'b',T,lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,lmb(:,3),'k',T,lmb(:,5),'r',T,lmb(:,7),'b',T,lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 29
        %ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', lmblnM1,...
        %    ubAcnM1.', ubpAcnM1.', ubppAcnM1.', lmbAcnM1.',...
        %    xTilnM1.', xTilpnM1.', xTilppnM1.', lmbnM1.'];
        ubl = ZZ(:,1:3);
        ubpl = ZZ(:,3+1:2*3);
        ubppl = ZZ(:,2*3+1:3*3);
        lmbl = kNmk*ZZ(:,3*3+1);
        ubAc = ZZ(:,3*3+1+1:3*3+1+nAc);
        ubpAc = ZZ(:,3*3+1+nAc+1:3*3+1+2*nAc);
        ubppAc = ZZ(:,3*3+1+2*nAc+1:3*3+1+3*nAc);
        lmbAc = kNmk*ZZ(:,3*3+1+3*nAc+1:3*3+1+3*nAc+2);
        ubEf = ZZ(:,3*3+1+3*nAc+2+1:3*3+1+3*nAc+2+Ng);
        ubpEf = ZZ(:,3*3+1+3*nAc+2+Ng+1:3*3+1+3*nAc+2+2*Ng);
        ubppEf = ZZ(:,3*3+1+3*nAc+2+2*Ng+1:3*3+1+3*nAc+2+3*Ng);
        lmbEF = kNmk*ZZ(:,3*3+1+3*nAc+2+3*Ng+1:3*3+1+3*nAc+2+3*Ng+rEstr+2);
        
        ub1long = ubl(:,1);
        ub2long = ubl(:,2);
        ub3long = ubl(:,3);
        y2long = fNub3(ub3long);
        
        y1Ac = ubpAc(:,nAc);
        ub1Ac = ubAc(:,1);
        ub2Ac = ubAc(:,3);
        ub3Ac = ubAc(:,7);
        y2Ac = fNbu3Ac(ub3Ac);
        vL = ubAc(:,4);
        wL = ubAc(:,5);
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = ubpEf(:,6*(N+1));
        ub1 = ubEf(:,6*(1)-5);
        ubN1M1 = ubEf(:,6*(N1+1)-5);
        ubNM1 = ubEf(:,6*(N+1)-5);
        y2 = FnbNM1(ubNM1);
        vN1M1 = ubEf(:,6*(N1+1)-4);
        wN1M1 = ubEf(:,6*(N1+1)-2);
        uN1M1 = sqrt(vN1M1.^2 + wN1M1.^2);
        
        %Entradas:
        figure
        plot(T,lmbl,'k',T,lmbAc(:,1),'b',T,lmbEF(:,1),'r')
        xlabel('Tempo (s)')
        ylabel('LmbU (N)')
        legend('PC long','PC c/ acopl','EF')
        %title('Força de reação')
        grid
        figure
        plot(T,lmbAc(:,2),'b',T,lmbEF(:,2),'r')
        xlabel('Tempo (s)')
        ylabel('LmbtX (Nm)')
        legend('PC c/ acopl','EF')
        %title('Força de reação')
        grid
        %Saídas:
        figure
        plot(T,y1Ac*30/pi,'b',T,y1*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        %title('Velocidade TetaXp do rotor principal')
        legend('PC c/ acopl','EF')
        grid
        figure
        plot(T,y2long,'k',T,y2Ac,'b',T,y2,'r')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        %title('Força normal de contato')
        legend('PC long','PC c/ acopl','EF')
        grid
        %Posições longitudinais:
        figure
        plot(T,ub1long*10^3,'k',T,ub1Ac*10^3,'b',T,ub1*10^3,'r',...
             T,ub2long*10^3,'k',T,ub2Ac*10^3,'b',T,ub2*10^3,'r',...
             T,ub3long*10^3,'k',T,ub3Ac*10^3,'b',T,ub3*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        %title('Posições longitudinais (u)')
        legend('PC long','PC c/ acopl','EF')
        grid
        %Posições laterais:
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        title('PC c/ acopl')
        legend('vL','wL','uL')
        grid
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        title('PC c/ acopl')
        legend('órbita','folga')
        grid
        figure
        plot(T,vN1M1*10^3,'b',T,wN1M1*10^3,'r',T,uN1M1*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        title('EF')
        legend('vN1M1','wN1M1','uN1M1')
        grid
        figure
        plot(vN1M1*10^3,wN1M1*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yN1M1 (mm)')
        ylabel('zN1M1 (mm)')
        title('EF')
        axis equal
        legend('órbita','folga')
        grid
        
        %Reações referentes às restrições (EF):
        figure
        plot(T,lmbEF(:,3),'k',T,lmbEF(:,5),'r',T,lmbEF(:,7),'b',T,lmbEF(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,lmbEF(:,4),'k',T,lmbEF(:,6),'r',T,lmbEF(:,8),'b',T,lmbEF(:,10),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 30
        %ZZ(i,:) = [ubAcnM1.', ubpAcnM1.', ubppAcnM1.', lmbAcnM1.',...
        %    ubEFnM1.', ubpEFnM1.', ubppEFnM1.', lmbnM1.',...
        %    ubrknM1.', ubprknM1.', ubpprknM1.', ub1nM1, ubp1nM1, ubpp1nM1, TxnM1];
        ubAc = ZZ(:,1:nAc);
        ubpAc = ZZ(:,nAc+1:2*nAc);
        ubppAc = ZZ(:,2*nAc+1:3*nAc);
        lmbAc = kNmk*ZZ(:,3*nAc+1:3*nAc+2);
        ubEf = ZZ(:,3*nAc+2+1:3*nAc+2+Ng);
        ubpEf = ZZ(:,3*nAc+2+Ng+1:3*nAc+2+2*Ng);
        ubppEf = ZZ(:,3*nAc+2+2*Ng+1:3*nAc+2+3*Ng);
        lmbEF = kNmk*ZZ(:,3*nAc+2+3*Ng+1:3*nAc+2+3*Ng+rEstr+2);
        ubrk = ZZ(:,3*nAc+2+3*Ng+rEstr+2+1:3*nAc+2+3*Ng+rEstr+2+cPk);
        ubprk = ZZ(:,3*nAc+2+3*Ng+rEstr+2+cPk+1:3*nAc+2+3*Ng+rEstr+2+2*cPk);
        ubpprk = ZZ(:,3*nAc+2+3*Ng+rEstr+2+2*cPk+1:3*nAc+2+3*Ng+rEstr+2+3*cPk);
        
        y1Ac = ubpAc(:,nAc);
        ub1Ac = ubAc(:,1);
        ub2Ac = ubAc(:,3);
        ub3Ac = ubAc(:,7);
        y2Ac = fNbu3Ac(ub3Ac);
        vL = ubAc(:,4);
        wL = ubAc(:,5);
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = ubpEf(:,6*(N+1));
        ub1 = ubEf(:,6*(1)-5);
        ubN1M1 = ubEf(:,6*(N1+1)-5);
        ubNM1 = ubEf(:,6*(N+1)-5);
        y2 = FnbNM1(ubNM1);
        vN1M1 = ubEf(:,6*(N1+1)-4);
        wN1M1 = ubEf(:,6*(N1+1)-2);
        uN1M1 = sqrt(vN1M1.^2 + wN1M1.^2);
        
        Zrk = [ubrk, ubprk];
        
        ub1rk = Zrk*Cu1.';
        tpX1rk = Zrk*CtpX1.';
        ubN1M1rk = Zrk*CuN1M1.';
        vN1M1rk = Zrk*CvN1M1.';
        wN1M1rk = Zrk*CwN1M1.';
        tpXN1M1rk = Zrk*CtpXN1M1.';
        ubNM1rk = Zrk*CuNM1.';
        tpXNM1rk = Zrk*CtpXNM1.';
        
        %Entradas:
        figure
        plot(T,lmbAc(:,1),'b',T,lmbEF(:,1),'r')
        xlabel('Tempo (s)')
        ylabel('LmbU (N)')
        legend('PC c/ acopl','EF')
        %title('Força de reação')
        grid
        figure
        plot(T,lmbAc(:,2),'b',T,lmbEF(:,2),'r')
        xlabel('Tempo (s)')
        ylabel('LmbtX (Nm)')
        legend('PC c/ acopl','EF')
        %title('Força de reação')
        grid
        %Saídas:
        figure
        plot(T,y1Ac*30/pi,'b',T,y1*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        %title('Velocidade TetaXp do rotor principal')
        legend('PC c/ acopl','EF')
        grid
        figure
        plot(T,y2Ac,'k',T,y2,'r')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        %title('Força normal de contato')
        legend('PC c/ acopl','EF')
        grid
        %Posições longitudinais:
        figure
        plot(T,ub1Ac*10^3,'k',T,ub1*10^3,'r',...
             T,ub2Ac*10^3,'k',T,ubN1M1*10^3,'r',...
             T,ub3Ac*10^3,'k',T,ubNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        %title('Posições longitudinais (u)')
        legend('PC c/ acopl','EF')
        grid
        figure
        plot(T,ub3Ac*10^3,'k',T,ubNM1*10^3,'r',T,ubNM1rk*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('ubNM1 (mm)')
        legend('PC c/ acopl','EF','EF red')
        grid
        %Posições laterais:
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        title('PC c/ acopl')
        legend('vL','wL','uL')
        grid
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        title('PC c/ acopl')
        legend('órbita','folga')
        grid
        figure
        plot(T,vN1M1*10^3,'b',T,wN1M1*10^3,'r',T,uN1M1*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        title('EF')
        legend('vN1M1','wN1M1','uN1M1')
        grid
        figure
        plot(vN1M1*10^3,wN1M1*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yN1M1 (mm)')
        ylabel('zN1M1 (mm)')
        title('EF')
        axis equal
        legend('órbita','folga')
        grid
        
        %Reações referentes às restrições (EF):
        figure
        plot(T,lmbEF(:,3),'k',T,lmbEF(:,5),'r',T,lmbEF(:,7),'b',T,lmbEF(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,lmbEF(:,4),'k',T,lmbEF(:,6),'r',T,lmbEF(:,8),'b',T,lmbEF(:,10),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 31
        %ZZ(i,:) = [ubrknM1.', ubprknM1.', ubpprknM1.', lmbnM1, TxnM1];
        ubk = ZZ(:,1:cPk);
        ubpk = ZZ(:,cPk+1:2*cPk);
        ubppk = ZZ(:,2*cPk+1:3*cPk);
        lmbrk = kNmk*ZZ(:,3*cPk+1);
        Tx = ZZ(:,3*cPk+2);
        Zrk = [ubk, ubpk];
        
        ub1 = Zrk*Cu1.';
        tpX1 = Zrk*CtpX1.';
        ubN1M1 = Zrk*CuN1M1.';
        vN1M1 = Zrk*CvN1M1.';
        wN1M1 = Zrk*CwN1M1.';
        tpXN1M1 = Zrk*CtpXN1M1.';
        ubNM1 = Zrk*CuNM1.';
        tpXNM1 = Zrk*CtpXNM1.';
        
        uL = sqrt(vN1M1.^2 + wN1M1.^2);
        teta_fg = 0:0.001:(2*pi);
        
        %Entradas:
        figure
        plot(T,lmbrk,'k')
        xlabel('Tempo (s)')
        ylabel('lmbU1 (N)')
        grid
        figure
        plot(T,Tx,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        grid
        %Posições longitudinais:
        figure
        plot(T,ub1*10^3,'k',T,ubN1M1*10^3,'r',T,ubNM1*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ubN1M1','ubN1M1')
        grid
        %Posições laterais:
        figure
        plot(T,vN1M1*10^3,'b',T,wN1M1*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        title('EF')
        legend('vN1M1','wN1M1','uN1M1')
        grid
        figure
        plot(vN1M1*10^3,wN1M1*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yN1M1 (mm)')
        ylabel('zN1M1 (mm)')
        title('EF')
        axis equal
        legend('órbita','folga')
        grid
        %Velocidades de rotação:
        figure
        plot(T,tpX1*30/pi,'k',T,tpXN1M1*30/pi,'r',T,tpXNM1*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tpX (RPM)')
        title('Velocidades de rotação (tpX)')
        legend('tpX1','tpXN1M1','tpXNM1')
        grid
    case 32
        %ZZ(i,:) = [ubrknM1.', ubprknM1.', ubpprknM1.', ub1nM1, ubp1nM1, ubpp1nM1, TxnM1];
        ubk = ZZ(:,1:cPk);
        ubpk = ZZ(:,cPk+1:2*cPk);
        ubppk = ZZ(:,2*cPk+1:3*cPk);
        ub1 = ZZ(:,3*cPk+1);
        ubp1 = ZZ(:,3*cPk+2);
        ubpp1 = ZZ(:,3*cPk+3);
        Tx = ZZ(:,3*cPk+4);
        Zrk = [ubk, ubpk];
        
        ub1test = Zrk*Cu1.';
        tpX1 = Zrk*CtpX1.';
        ubN1M1 = Zrk*CuN1M1.';
        vN1M1 = Zrk*CvN1M1.';
        wN1M1 = Zrk*CwN1M1.';
        tpXN1M1 = Zrk*CtpXN1M1.';
        ubNM1 = Zrk*CuNM1.';
        tpXNM1 = Zrk*CtpXNM1.';
        
        uL = sqrt(vN1M1.^2 + wN1M1.^2);
        teta_fg = 0:0.001:(2*pi);
        
        %Entradas:
        figure
        plot(T,Tx,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        grid
        figure
        plot(T,ub1,'k')
        xlabel('Tempo (s)')
        ylabel('ub1In (mm)')
        grid
        figure
        plot(T,ubp1,'k')
        xlabel('Tempo (s)')
        ylabel('ubp1In (mm/s)')
        grid
        figure
        plot(T,ubpp1,'k')
        xlabel('Tempo (s)')
        ylabel('ubpp1In (mm/s2)')
        grid
        %Posições longitudinais:
        figure
        plot(T,ub1test*10^3,'k',T,ubN1M1*10^3,'r',T,ubNM1*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1Out','ubN1M1','ubNM1')
        grid
        %Posições laterais:
        figure
        plot(T,vN1M1*10^3,'b',T,wN1M1*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        title('EF')
        legend('vN1M1','wN1M1','uN1M1')
        grid
        figure
        plot(vN1M1*10^3,wN1M1*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yN1M1 (mm)')
        ylabel('zN1M1 (mm)')
        title('EF')
        axis equal
        legend('órbita','folga')
        grid
        %Velocidades de rotação:
        figure
        plot(T,tpX1*30/pi,'k',T,tpXN1M1*30/pi,'r',T,tpXNM1*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tpX (RPM)')
        title('Velocidades de rotação (tpX)')
        legend('tpX1','tpXN1M1','tpXNM1')
        grid
    case 33
        tX1 = ZZ(:,1);
        tX2 = ZZ(:,2);
        tX3 = ZZ(:,3);
        tpX1 = ZZ(:,4);
        tpX2 = ZZ(:,5);
        tpX3 = ZZ(:,6);
        figure
        plot(T,tX1,'k',T,tX2,'r',T,tX3,'b')
        xlabel('Tempo (s)')
        ylabel('tX (rad)')
        legend('tX1','tX2','tX3')
        grid
        figure
        plot(T,tpX1*(30/pi),'k',T,tpX2*(30/pi),'r',T,tpX3*(30/pi),'b')
        xlabel('Tempo (s)')
        ylabel('tpX (RPM)')
        legend('tpX1','tpX2','tpX3')
        grid
    case 34
        tX1 = fTeta1(T);
        tX2 = ZZ(:,1);
        tX3 = ZZ(:,2);
        tpX1 = fTetap1(T);
        tpX2 = ZZ(:,3);
        tpX3 = ZZ(:,4);
        figure
        plot(T,tX1,'k',T,tX2,'r',T,tX3,'b')
        xlabel('Tempo (s)')
        ylabel('tX (rad)')
        legend('tX1','tX2','tX3')
        grid
        figure
        plot(T,tpX1*(30/pi),'k',T,tpX2*(30/pi),'r',T,tpX3*(30/pi),'b')
        xlabel('Tempo (s)')
        ylabel('tpX (RPM)')
        legend('tpX1','tpX2','tpX3')
        grid
    case 35
        %ZZ = zeros(3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmb0));
        %ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmbnM1.'];
        nx = length(xTil0);
        nyd = length(y1dTil0);
        nZ = length(Zd0);
        nTx = length(Tx0);
        ns = length(s10);
        nlmb = length(lmb0);
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        y1dVt = ZZ(:,3*nx+1:3*nx+nyd);
        ZdDInv = ZZ(:,3*nx+nyd+1:3*nx+nyd+nZ);
        TxIn = ZZ(:,3*nx+nyd+nZ+1);
        s1Tx = ZZ(:,3*nx+nyd+nZ+2);
        Lmb = kNmk*ZZ(:,3*nx+nyd+nZ+2+1:3*nx+nyd+nZ+2+nlmb);
        
        ubEF = xTil(:,1:Ng);
        midVc = xTil(:,Ng+1);
        FiBdVc = xTil(:,Ng+2);
        ubpEF = xTilp(:,1:Ng);
        mipdVc = xTilp(:,Ng+1);
        FipBdVc = xTilp(:,Ng+2);
        ubppEF = xTilpp(:,1:Ng);
        mippdVc = xTilpp(:,Ng+1);
        FippBdVc = xTilpp(:,Ng+2);
        
        uTilT = ubEF*Atu;
        uTilpT = ubpEF*Atu;
        uTilppT = ubppEF*Atu;
        ubTilU = ubEF*Auu;
        ubpTilU = ubpEF*Auu;
        ubppTilU = ubppEF*Auu;
        vL = ubEF*Avu;
        wL = ubEF*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTil(:,3);
        y2mm = fNub3(ub3);
        ubNM1 = ubEF(:,6*(N+1)-5);
        y2 = FnbNM1(ubNM1);
        teta1d = ZdDInv(:,1);
        teta2d = ZdDInv(:,2);
        teta3d = ZdDInv(:,3);
        tetap1d = ZdDInv(:,4);
        tetap2d = ZdDInv(:,5);
        tetap3d = ZdDInv(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,Lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbUb1 (N)')
        title('Força de reação longitudinal')
        grid
        figure
        plot(T,TxIn,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k',T,y2mm,'r')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2')
        legend('Planta','Metamodelo')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBdVc,'k',T,-FiBdVc,'k',T,s1Tx,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBdVc,'k',T,FipBdVc,'r',T,FippBdVc,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,midVc,'k',T,mipdVc,'r',T,mippdVc,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,Lmb(:,2),'k',T,Lmb(:,4),'r',T,Lmb(:,6),'b',T,Lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,Lmb(:,3),'k',T,Lmb(:,5),'r',T,Lmb(:,7),'b',T,Lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 36
        %ZZ = zeros(3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmb0));
        %ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmbnM1.'];
        nx = length(xTil0);
        nyd = length(y1dTil0);
        nZ = length(Zd0);
        nTx = length(Tx0);
        ns = length(s10);
        nlmb = length(lmb0);
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        y1dVt = ZZ(:,3*nx+1:3*nx+nyd);
        ZdDInv = ZZ(:,3*nx+nyd+1:3*nx+nyd+nZ);
        TxIn = ZZ(:,3*nx+nyd+nZ+1);
        s1Tx = ZZ(:,3*nx+nyd+nZ+2);
        Lmb = kNmk*ZZ(:,3*nx+nyd+nZ+2+1:3*nx+nyd+nZ+2+nlmb);
        
        ubEF = xTil(:,1:Ng);
        midVc = xTil(:,Ng+1);
        FiBdVc = xTil(:,Ng+2);
        ubpEF = xTilp(:,1:Ng);
        mipdVc = xTilp(:,Ng+1);
        FipBdVc = xTilp(:,Ng+2);
        ubppEF = xTilpp(:,1:Ng);
        mippdVc = xTilpp(:,Ng+1);
        FippBdVc = xTilpp(:,Ng+2);
        
        uTilT = ubEF*Atu;
        uTilpT = ubpEF*Atu;
        uTilppT = ubppEF*Atu;
        ubTilU = ubEF*Auu;
        ubpTilU = ubpEF*Auu;
        ubppTilU = ubppEF*Auu;
        vL = ubEF*Avu;
        wL = ubEF*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ubNM1 = ubEF(:,6*(N+1)-5);
        y2 = FnbNM1(ubNM1);
        teta1d = ZdDInv(:,1);
        teta2d = ZdDInv(:,2);
        teta3d = ZdDInv(:,3);
        tetap1d = ZdDInv(:,4);
        tetap2d = ZdDInv(:,5);
        tetap3d = ZdDInv(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,Lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbUb1 (N)')
        title('Força de reação longitudinal')
        grid
        figure
        plot(T,TxIn,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        grid
        figure
        plot(T,y2,'k')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2')
        legend('Planta','Metamodelo')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Deslocamentos longitudinais:
        figure
        plot(T,ubTilU(:,1),'k',T,ubTilU(:,2),'r',T,ubTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBdVc,'k',T,-FiBdVc,'k',T,s1Tx,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBdVc,'k',T,FipBdVc,'r',T,FippBdVc,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,midVc,'k',T,mipdVc,'r',T,mippdVc,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,Lmb(:,2),'k',T,Lmb(:,4),'r',T,Lmb(:,6),'b',T,Lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,Lmb(:,3),'k',T,Lmb(:,5),'r',T,Lmb(:,7),'b',T,Lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 37
        %ZZ(i,:) = [ubTilUnM1.', ubTilpUnM1.', ubTilppUnM1.', xTilnM1.', xTilpnM1.', xTilppnM1.', y1dTilnM1.', ZdnM1.', TxnM1, s1nM1, lmblnM1, lmbnM1];
        nu = length(ubTilU0);
        nx = length(xTil0);
        nyd = length(y1dTil0);
        nZ = length(Zd0);
        nTx = length(Tx0);
        ns = length(s10);
        nlmbl = length(lmbl0);
        nlmb = length(lmb0);
        ubTil = ZZ(:,1:nu);
        ubpTil = ZZ(:,nu+1:2*nu);
        ubppTil = ZZ(:,2*nu+1:3*nu);
        xTil = ZZ(:,3*nu+1:3*nu+nx);
        xTilp = ZZ(:,3*nu+nx+1:3*nu+2*nx);
        xTilpp = ZZ(:,3*nu+2*nx+1:3*nu+3*nx);
        y1dVt = ZZ(:,3*nu+3*nx+1:3*nu+3*nx+nyd);
        ZdDInv = ZZ(:,3*nu+3*nx+nyd+1:3*nu+3*nx+nyd+nZ);
        TxIn = ZZ(:,3*nu+3*nx+nyd+nZ+1);
        s1Tx = ZZ(:,3*nu+3*nx+nyd+nZ+2);
        Lmbl = kNmk*ZZ(:,3*nu+3*nx+nyd+nZ+2+1:3*nu+3*nx+nyd+nZ+2+nlmbl);
        Lmb = kNmk*ZZ(:,3*nu+3*nx+nyd+nZ+2+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb);
        
        ubAc = xTil(:,1:8);
        midVc = xTil(:,8+1);
        FiBdVc = xTil(:,8+2);
        ubpAc = xTilp(:,1:8);
        mipdVc = xTilp(:,8+1);
        FipBdVc = xTilp(:,8+2);
        ubppAc = xTilpp(:,1:8);
        mippdVc = xTilpp(:,8+1);
        FippBdVc = xTilpp(:,8+2);
        
        uTilT = ubAc*Atu;
        uTilpT = ubpAc*Atu;
        uTilppT = ubppAc*Atu;
        ubTilU = ubAc*Auu;
        ubpTilU = ubpAc*Auu;
        ubppTilU = ubppAc*Auu;
        vL = ubAc*Avu;
        wL = ubAc*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTil(:,3);
        y2mm = fNub3(ub3);
        ub3Ac = ubTilU(:,3);
        y2Ac = fNbu3Ac(ub3Ac);
        
        teta1d = ZdDInv(:,1);
        teta2d = ZdDInv(:,2);
        teta3d = ZdDInv(:,3);
        tetap1d = ZdDInv(:,4);
        tetap2d = ZdDInv(:,5);
        tetap3d = ZdDInv(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,Lmbl(:,1),'k',T,Lmb(:,1),'r')
        xlabel('Tempo (s)')
        ylabel('LmbUb1 (N)')
        title('Força de reação longitudinal')
        legend('PC','PCAcpl')
        grid
        figure
        plot(T,TxIn,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k--',T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        legend('Traj real','Traj desejada')
        grid
        figure
        plot(T,y2Ac,'k',T,y2mm,'k--')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2')
        legend('Planta','Metamodelo')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBdVc,'k',T,-FiBdVc,'k',T,s1Tx,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBdVc,'k',T,FipBdVc,'r',T,FippBdVc,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,midVc,'k',T,mipdVc,'r',T,mippdVc,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2Ac,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
    case 38
        %ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmbl0) + length(lmb0));
        %ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmbl0, lmb0.'];
        nu = length(ubTilU0);
        nx = length(xTil0);
        nyd = length(y1dTil0);
        nZ = length(Zd0);
        nTx = length(Tx0);
        ns = length(s10);
        nlmbl = length(lmbl0);
        nlmb = length(lmb0);
        ubTil = ZZ(:,1:nu);
        ubpTil = ZZ(:,nu+1:2*nu);
        ubppTil = ZZ(:,2*nu+1:3*nu);
        xTil = ZZ(:,3*nu+1:3*nu+nx);
        xTilp = ZZ(:,3*nu+nx+1:3*nu+2*nx);
        xTilpp = ZZ(:,3*nu+2*nx+1:3*nu+3*nx);
        y1dVt = ZZ(:,3*nu+3*nx+1:3*nu+3*nx+nyd);
        ZdDInv = ZZ(:,3*nu+3*nx+nyd+1:3*nu+3*nx+nyd+nZ);
        TxIn = ZZ(:,3*nu+3*nx+nyd+nZ+1);
        s1Tx = ZZ(:,3*nu+3*nx+nyd+nZ+2);
        Lmb = kNmk*ZZ(:,3*nu+3*nx+nyd+nZ+2+1:3*nu+3*nx+nyd+nZ+2+nlmb);
        
        ubEF = xTil(:,1:Ng);
        midVc = xTil(:,Ng+1);
        FiBdVc = xTil(:,Ng+2);
        ubpEF = xTilp(:,1:Ng);
        mipdVc = xTilp(:,Ng+1);
        FipBdVc = xTilp(:,Ng+2);
        ubppEF = xTilpp(:,1:Ng);
        mippdVc = xTilpp(:,Ng+1);
        FippBdVc = xTilpp(:,Ng+2);
        
        uTilT = ubEF*Atu;
        uTilpT = ubpEF*Atu;
        uTilppT = ubppEF*Atu;
        ubTilU = ubEF*Auu;
        ubpTilU = ubpEF*Auu;
        ubppTilU = ubppEF*Auu;
        vL = ubEF*Avu;
        wL = ubEF*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTil(:,3);
        y2mm = fNub3(ub3);
        ubNM1 = ubEF(:,6*(N+1)-5);
        y2 = FnbNM1(ubNM1);
        teta1d = ZdDInv(:,1);
        teta2d = ZdDInv(:,2);
        teta3d = ZdDInv(:,3);
        tetap1d = ZdDInv(:,4);
        tetap2d = ZdDInv(:,5);
        tetap3d = ZdDInv(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,Lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbUb1 (N)')
        title('Força de reação longitudinal')
        grid
        figure
        plot(T,TxIn,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k--',T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        legend('Traj real','Traj desejada')
        grid
        figure
        plot(T,y2,'k',T,y2mm,'k--')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2')
        legend('Planta','Metamodelo')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBdVc,'k',T,-FiBdVc,'k',T,s1Tx,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBdVc,'k',T,FipBdVc,'r',T,FippBdVc,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,midVc,'k',T,mipdVc,'r',T,mippdVc,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        title('Mapa de operação')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,Lmb(:,2),'k',T,Lmb(:,4),'r',T,Lmb(:,6),'b',T,Lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        title('Reações de força nas extremidades')
        grid
        figure
        plot(T,Lmb(:,3),'k',T,Lmb(:,5),'r',T,Lmb(:,7),'b',T,Lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        title('Reações de momento nas extremidades')
        grid
    case 39
        %ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmbl0) + length(lmb0) + 3*length(FiTil0));
        %ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmbl0, lmb0.', FiTil0.', FipTil0.', FippTil0.'];
        nu = length(ubTilU0);
        nx = length(xTil0);
        nyd = length(y1dTil0);
        nZ = length(Zd0);
        nTx = length(Tx0);
        ns = length(s10);
        nlmbl = length(lmbl0);
        nlmb = length(lmb0);
        nFiTil = length(FiTil0);
        ubTil = ZZ(:,1:nu);
        ubpTil = ZZ(:,nu+1:2*nu);
        ubppTil = ZZ(:,2*nu+1:3*nu);
        xTil = ZZ(:,3*nu+1:3*nu+nx);
        xTilp = ZZ(:,3*nu+nx+1:3*nu+2*nx);
        xTilpp = ZZ(:,3*nu+2*nx+1:3*nu+3*nx);
        y1dVt = ZZ(:,3*nu+3*nx+1:3*nu+3*nx+nyd);
        ZdDInv = ZZ(:,3*nu+3*nx+nyd+1:3*nu+3*nx+nyd+nZ);
        TxIn = ZZ(:,3*nu+3*nx+nyd+nZ+1);
        s1Tx = ZZ(:,3*nu+3*nx+nyd+nZ+2);
        Lmbl = kNmk*ZZ(:,3*nu+3*nx+nyd+nZ+2+1:3*nu+3*nx+nyd+nZ+2+nlmbl);
        Lmb = kNmk*ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb);
        FiEuler = ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+nFiTil);
        FipEuler = ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+nFiTil+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+2*nFiTil);
        FippEuler = ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+2*nFiTil+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+3*nFiTil);
        %y1nm1 = ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+3*nFiTil+1); %y1 no passo anterior
        
        ubEF = xTil(:,1:Ng);
        midVc = xTil(:,Ng+1);
        FiBdVc = xTil(:,Ng+2);
        ubpEF = xTilp(:,1:Ng);
        mipdVc = xTilp(:,Ng+1);
        FipBdVc = xTilp(:,Ng+2);
        ubppEF = xTilpp(:,1:Ng);
        mippdVc = xTilpp(:,Ng+1);
        FippBdVc = xTilpp(:,Ng+2);
        
        uTilT = ubEF*Atu;
        uTilpT = ubpEF*Atu;
        uTilppT = ubppEF*Atu;
        ubTilU = ubEF*Auu;
        ubpTilU = ubpEF*Auu;
        ubppTilU = ubppEF*Auu;
        vL = ubEF*Avu;
        wL = ubEF*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTil(:,3);
        y2mm = fNub3(ub3);
        ubNM1 = ubEF(:,6*(N+1)-5);
        y2 = FnbNM1(ubNM1);
        teta1d = ZdDInv(:,1);
        teta2d = ZdDInv(:,2);
        teta3d = ZdDInv(:,3);
        tetap1d = ZdDInv(:,4);
        tetap2d = ZdDInv(:,5);
        tetap3d = ZdDInv(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,Lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbUb1 (N)')
        title('Força de reação longitudinal')
        grid
        figure
        plot(T,TxIn,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        %plot(T,y1*(30/pi),'k--',T,y1dVt(:,1)*30/pi,'k',T,y1nm1*(30/pi),'r')
        plot(T,y1*(30/pi),'k--',T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        legend('Traj real','Traj desejada')
        grid
        figure
        plot(T,y2,'k',T,y2mm,'k--')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2')
        legend('Planta','Metamodelo')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBdVc,'k',T,-FiBdVc,'k',T,s1Tx,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBdVc,'k',T,FipBdVc,'r',T,FippBdVc,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,midVc,'k',T,mipdVc,'r',T,mippdVc,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        %title('Mapa de operação')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,Lmb(:,2),'k',T,Lmb(:,4),'r',T,Lmb(:,6),'b',T,Lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        %title('Reações de força nas extremidades')
        grid
        figure
        plot(T,Lmb(:,3),'k',T,Lmb(:,5),'r',T,Lmb(:,7),'b',T,Lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        %title('Reações de momento nas extremidades')
        grid
        %Ângulos de Euler para o rotor intermediario:
        figure
        plot(T,FipEuler(:,1)*30/pi,'k',T,FipEuler(:,2)*30/pi,'r',T,FipEuler(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('RPM')
        legend('precessão','nutação','spin')
        grid
    case 40
        %ZZ = zeros(nT,3*length(ubTilU0) + 3*length(xTil0) + length(y1dTil0) + length(Zd0) + length(Tx0) + length(s10) + length(lmbl0) + length(lmb0) + 3*length(FiTil0) + length(y10));
        %ZZ(1,:) = [ubTilU0.', ubTilpU0.', ubTilppU0.', xTil0.', xTilp0.', xTilpp0.', y1dTil0.', Zd0.', Tx0, s10, lmbl0, lmb0.', FiTil0.', FipTil0.', FippTil0.', y10];
        nu = length(ubTilU0);
        nx = length(xTil0);
        nyd = length(y1dTil0);
        nZ = length(Zd0);
        nTx = length(Tx0);
        ns = length(s10);
        nlmbl = length(lmbl0);
        nlmb = length(lmb0);
        nFiTil = length(FiTil0);
        ny10 = length(y10);
        ubTil = ZZ(:,1:nu);
        ubpTil = ZZ(:,nu+1:2*nu);
        ubppTil = ZZ(:,2*nu+1:3*nu);
        xTil = ZZ(:,3*nu+1:3*nu+nx);
        xTilp = ZZ(:,3*nu+nx+1:3*nu+2*nx);
        xTilpp = ZZ(:,3*nu+2*nx+1:3*nu+3*nx);
        y1dVt = ZZ(:,3*nu+3*nx+1:3*nu+3*nx+nyd);
        ZdDInv = ZZ(:,3*nu+3*nx+nyd+1:3*nu+3*nx+nyd+nZ);
        TxIn = ZZ(:,3*nu+3*nx+nyd+nZ+1);
        s1Tx = ZZ(:,3*nu+3*nx+nyd+nZ+2);
        Lmbl = kNmk*ZZ(:,3*nu+3*nx+nyd+nZ+2+1:3*nu+3*nx+nyd+nZ+2+nlmbl);
        Lmb = kNmk*ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb);
        FiEuler = ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+nFiTil);
        FipEuler = ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+nFiTil+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+2*nFiTil);
        FippEuler = ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+2*nFiTil+1:3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+3*nFiTil);
        y1nm1 = ZZ(:,3*nu+3*nx+nyd+nZ+2+nlmbl+nlmb+3*nFiTil+1); %y1 no passo anterior
        
        ubEF = xTil(:,1:Ng);
        midVc = xTil(:,Ng+1);
        FiBdVc = xTil(:,Ng+2);
        ubpEF = xTilp(:,1:Ng);
        mipdVc = xTilp(:,Ng+1);
        FipBdVc = xTilp(:,Ng+2);
        ubppEF = xTilpp(:,1:Ng);
        mippdVc = xTilpp(:,Ng+1);
        FippBdVc = xTilpp(:,Ng+2);
        
        uTilT = ubEF*Atu;
        uTilpT = ubpEF*Atu;
        uTilppT = ubppEF*Atu;
        ubTilU = ubEF*Auu;
        ubpTilU = ubpEF*Auu;
        ubppTilU = ubppEF*Auu;
        vL = ubEF*Avu;
        wL = ubEF*Awu;
        uL = sqrt(vL.^2 + wL.^2);
        
        y1 = [uTilT, uTilpT]*C1.';
        ub3 = ubTil(:,3);
        y2mm = fNub3(ub3);
        ubNM1 = ubEF(:,6*(N+1)-5);
        y2 = FnbNM1(ubNM1);
        teta1d = ZdDInv(:,1);
        teta2d = ZdDInv(:,2);
        teta3d = ZdDInv(:,3);
        tetap1d = ZdDInv(:,4);
        tetap2d = ZdDInv(:,5);
        tetap3d = ZdDInv(:,6);
        
        %%
        %Entradas:
        figure
        plot(T,Lmb(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('LmbUb1 (N)')
        title('Força de reação longitudinal')
        grid
        figure
        plot(T,TxIn,'k')
        xlabel('Tempo (s)')
        ylabel('Tx (Nm)')
        title('Torque de entrada')
        grid
        %Variáveis de saída:
        figure
        plot(T,y1*(30/pi),'k--',T,y1dVt(:,1)*30/pi,'k',T,y1nm1*(30/pi),'r')
        xlabel('Tempo (s)')
        ylabel('y1 (RPM)')
        title('Variável de saída: y1 = teta3p')
        legend('Traj real','Traj desejada','Traj real aprox')
        grid
        figure
        plot(T,y2,'k',T,y2mm,'k--')
        xlabel('Tempo (s)')
        ylabel('y2 (N)')
        title('Variável de saída: y2')
        legend('Planta','Metamodelo')
        grid
        %Posições longitudinais:
        figure
        plot(T,ubTilU(:,1)*10^3,'k',T,ubTilU(:,2)*10^3,'r',T,ubTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('u (mm)')
        title('Posições longitudinais (u)')
        legend('ub1','ub2','ub3')
        grid
        %Velocidades longitudinais:
        figure
        plot(T,ubpTilU(:,1)*10^3,'k',T,ubpTilU(:,2)*10^3,'r',T,ubpTilU(:,3)*10^3,'b')
        xlabel('Tempo (s)')
        ylabel('up (mm/s)')
        title('Velocidades longitudinais (up)')
        legend('up1','up2','up3')
        grid
        %Acelerações longitudinais:
        figure
        plot(T,ubppTilU(:,1),'k',T,ubppTilU(:,2),'r',T,ubppTilU(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('upp (mm/s2)')
        title('Acelerações longitudinais (upp)')
        legend('upp1','upp2','upp3')
        grid
        %%
        %Deslocamentos angulares:
        figure
        plot(T,uTilT(:,1),'k',T,uTilT(:,2),'r',T,uTilT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Deslocamento angular - Torção')
        grid
        %Velocidades angulares:
        figure
        plot(T,uTilpT(:,1)*30/pi,'k',T,uTilpT(:,2)*30/pi,'r',T,uTilpT(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        title('Velocidade angular - Torção')
        grid
        %Acelerações angulares:
        figure
        plot(T,uTilppT(:,1),'k',T,uTilppT(:,2),'r',T,uTilppT(:,3),'b')
        xlabel('Tempo (s)')
        ylabel('teta (rad/s2)')
        title('Aceleração angular - Torção')
        grid
        %Camada limite e variável de escorregamento:
        figure
        plot(T,FiBdVc,'k',T,-FiBdVc,'k',T,s1Tx,'r')
        xlabel('Tempo (s)')
        ylabel('FiBd, s1')
        title('Camada limite e variável de escorregamento - Torção')
        grid
        figure
        plot(T,FiBdVc,'k',T,FipBdVc,'r',T,FippBdVc,'b')
        xlabel('Tempo (s)')
        %ylabel('FiBd, FipBd, FippBd')
        title('Camada limite e derivadas - Torção')
        legend('FiBd','FipBd','FippBd')
        grid
        %Dinâmica interna:
        figure
        plot(T,midVc,'k',T,mipdVc,'r',T,mippdVc,'b')
        xlabel('Tempo (s)')
        ylabel('Mid, Mipd, Mippd')
        title('Dinâmica interma e derivadas - Torção')
        legend('Mid','Mipd','Mippd')
        grid
        %Trajetória desejada:
        figure
        plot(T,teta1d,'k',T,teta2d,'r',T,teta3d,'b')
        xlabel('Tempo (s)')
        ylabel('teta1d (rad)')
        title('Posição desejada - Torção - Zd')
        legend('teta1d','teta2d','teta3d')
        grid
        figure
        plot(T,tetap1d*30/pi,'k',T,tetap2d*30/pi,'r',T,tetap3d*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        title('Velocidade desejada - Torção - Zd')
        legend('tetap1d','tetap2d','tetap3d')
        grid
        figure
        plot(T,y1dVt(:,1)*30/pi,'k')
        xlabel('Tempo (s)')
        ylabel('y1d (RPM)')
        title('Saída desejada - Y1d')
        grid
        %Posições laterais:
        figure
        plot(T,vL*10^3,'b',T,wL*10^3,'r',T,uL*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('vL','wL','uL')
        grid
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(vL*10^3,wL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %%
        figure
        plot(y1*(30/pi),y2,'k')
        xlabel('tp3 (RPM)')
        ylabel('Nr (N)')
        %title('Mapa de operação')
        grid
        %Reações referentes às restrições:
        figure
        plot(T,Lmb(:,2),'k',T,Lmb(:,4),'r',T,Lmb(:,6),'b',T,Lmb(:,8),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (N)')
        legend('lmbV1','lmbW1','lmbVNM1','lmbWNM1')
        %title('Reações de força nas extremidades')
        grid
        figure
        plot(T,Lmb(:,3),'k',T,Lmb(:,5),'r',T,Lmb(:,7),'b',T,Lmb(:,9),'g')
        xlabel('Tempo (s)')
        ylabel('Reações (Nm)')
        legend('lmbTz1','lmbTy1','lmbTzNM1','lmbTyNM1')
        %title('Reações de momento nas extremidades')
        grid
        %Ângulos de Euler para o rotor intermediario:
        figure
        plot(T,FipEuler(:,1)*30/pi,'k',T,FipEuler(:,2)*30/pi,'r',T,FipEuler(:,3)*30/pi,'b')
        xlabel('Tempo (s)')
        ylabel('RPM')
        legend('precessão','nutação','spin')
        grid
    otherwise
        disp('other value')
end
%%
ad = 1;