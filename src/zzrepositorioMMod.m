function [ argumentos, repositorio ] = zzrepositorioMMod( argumentos )
%ÍNDICE Summary of this function goes here
%   T_gCj_Aux: linha 128
%   T_gamaBj_Aux: linha 134
%   T_m1j_Aux: linha 143
%   T_m2juppj_Aux: linha 151
%   U_kj_Aux: linha 171
%   U_fuj_Aux: linha 306
%   wj_Aux: linha 309
%   T_dgCj_dujT_Aux: linha 316
%   T_dgCj_dupjT_Aux: linha 325
%   T_dgamaBj_dujT_Aux: linha 331
%   T_dgamaBj_dupjT_Aux: linha 343
%   T_dm2juppj_dujT_Aux: linha 352
%   U_dfuj_dujT_Aux: linha 649

R = argumentos.raio;
E = argumentos.E;
G = argumentos.G;
ks = argumentos.ks;
rho = argumentos.rho;
%L = argumentos.L;
%L1 = argumentos.L1;
A = argumentos.Area;
Jp = argumentos.Jp;
J = argumentos.J;
alfa = argumentos.alfa_m2;
g = argumentos.g;
N = argumentos.N;
N0 = argumentos.N0;
N1 = argumentos.N1;
Ip1 = argumentos.Ip1;
I1 = Ip1/2;
M1 = argumentos.M1rt;
Ip2 = argumentos.Ip2b;
I2 = Ip2/2;
M2b = argumentos.M2b;
Ip3bMt = argumentos.Ip3bMt;
M3pl = argumentos.M3_pl;
I3pl = argumentos.I3_pl;
Ip3pl = argumentos.Ip3_pl;
dm = argumentos.dm;
l2 = argumentos.l2;
nGL = argumentos.nGL;
rEstr = argumentos.rEstr;
bUl = argumentos.bU1;
bTx1 = argumentos.bTx1;
bV = argumentos.bV;
bW = argumentos.bW;
bUNM1 = argumentos.bUNM1;
bTxNM1 = argumentos.bTxNM1;
%ok, 07/07/2022

%Integrando da energia potencial U:
Psi_x_linear = @(du_dx, dv_dx, dtetaZ_dx, dw_dx, dtetaY_dx, dtetaX_dx, tetaZ, tetaY)(...
    pi*R^2*E*...
    (du_dx^2) +...
    pi*R^4*E*...
    ((1/4)*(dtetaY_dx^2 + dtetaZ_dx^2)) +...
    pi*R^2*ks*G*...
    ((tetaZ - dv_dx)^2 + (tetaY + dw_dx)^2) +...
    pi*R^4*ks*G*...
    ((1/2)*...
    (dtetaX_dx^2)));
Psi_x_nlinear = @(du_dx, dv_dx, dtetaZ_dx, dw_dx, dtetaY_dx, dtetaX_dx, tetaZ, tetaY)(pi*E*R^2*...
    (du_dx^3 + (du_dx^4)/4 + (dv_dx^4)/4 + (dw_dx^4)/4 +...
     du_dx*dv_dx^2 + du_dx*dw_dx^2 +...
     (du_dx^2*dv_dx^2)/2 + (du_dx^2*dw_dx^2)/2 + (dv_dx^2*dw_dx^2)/2) +...
    pi*E*R^4*...
    ((tetaY^2*dtetaX_dx^2)/4 + (tetaZ^2*dtetaX_dx^2)/4 +...
     (dtetaX_dx^2*du_dx)/2 + (3*dtetaY_dx^2*du_dx)/4 + (3*dtetaZ_dx^2*du_dx)/4 +...
     (dtetaX_dx^2*du_dx^2)/4 + (3*dtetaY_dx^2*du_dx^2)/8 + (3*dtetaZ_dx^2*du_dx^2)/8 +...
     (dtetaX_dx^2*dv_dx^2)/2 + (dtetaY_dx^2*dv_dx^2)/8 + (dtetaZ_dx^2*dv_dx^2)/8 +...
     (dtetaX_dx^2*dw_dx^2)/2 + (dtetaY_dx^2*dw_dx^2)/8 + (dtetaZ_dx^2*dw_dx^2)/8 -...
     (dtetaX_dx*dtetaY_dx*dv_dx)/2 - (dtetaX_dx*dtetaZ_dx*dw_dx)/2 -...
     (tetaY*dtetaX_dx*dtetaZ_dx)/2 + (tetaZ*dtetaX_dx*dtetaY_dx)/2 +...
     (3*tetaY^2*dtetaX_dx^2*du_dx^2)/8 + (3*tetaZ^2*dtetaX_dx^2*du_dx^2)/8 +...
     (tetaY^2*dtetaX_dx^2*dv_dx^2)/8 + (tetaZ^2*dtetaX_dx^2*dv_dx^2)/8 +...
     (tetaY^2*dtetaX_dx^2*dw_dx^2)/8 + (tetaZ^2*dtetaX_dx^2*dw_dx^2)/8 -...
     (tetaZ*dtetaX_dx^2*dv_dx)/2 + (tetaY*dtetaX_dx^2*dw_dx)/2 +...
     (3*tetaY^2*dtetaX_dx^2*du_dx)/4 + (3*tetaZ^2*dtetaX_dx^2*du_dx)/4 -...
     (3*tetaY*dtetaX_dx*dtetaZ_dx*du_dx^2)/4 + (3*tetaZ*dtetaX_dx*dtetaY_dx*du_dx^2)/4 -...
     (tetaY*dtetaX_dx*dtetaZ_dx*dv_dx^2)/4 + (tetaZ*dtetaX_dx*dtetaY_dx*dv_dx^2)/4 -...
     (tetaZ*dtetaX_dx^2*du_dx*dv_dx)/2 - (tetaY*dtetaX_dx*dtetaZ_dx*dw_dx^2)/4 +...
     (tetaZ*dtetaX_dx*dtetaY_dx*dw_dx^2)/4 + (tetaY*dtetaX_dx^2*du_dx*dw_dx)/2 -...
     (dtetaX_dx*dtetaY_dx*du_dx*dv_dx)/2 - (dtetaX_dx*dtetaZ_dx*du_dx*dw_dx)/2 -...
     (3*tetaY*dtetaX_dx*dtetaZ_dx*du_dx)/2 + (3*tetaZ*dtetaX_dx*dtetaY_dx*du_dx)/2) +...
    pi*E*R^6*...
    ((dtetaX_dx^4)/12 + (dtetaY_dx^4)/32 + (dtetaZ_dx^4)/32 +...
     (tetaY^2*dtetaX_dx^4)/12 + (tetaZ^2*dtetaX_dx^4)/12 +...
     (tetaY^4*dtetaX_dx^4)/32 + (tetaZ^4*dtetaX_dx^4)/32 +...
     (dtetaX_dx^2*dtetaY_dx^2)/12 + (dtetaX_dx^2*dtetaZ_dx^2)/12 +...
     (dtetaY_dx^2*dtetaZ_dx^2)/16 + (tetaY^2*dtetaX_dx^2*dtetaY_dx^2)/16 +...
     (3*tetaY^2*dtetaX_dx^2*dtetaZ_dx^2)/16 + (3*tetaZ^2*dtetaX_dx^2*dtetaY_dx^2)/16 +...
     (tetaZ^2*dtetaX_dx^2*dtetaZ_dx^2)/16 + (tetaY^2*tetaZ^2*dtetaX_dx^4)/16 -...
     (tetaY*dtetaX_dx^3*dtetaZ_dx)/6 - (tetaY*dtetaX_dx*dtetaZ_dx^3)/8 +...
     (tetaZ*dtetaX_dx^3*dtetaY_dx)/6 + (tetaZ*dtetaX_dx*dtetaY_dx^3)/8 -...
     (tetaY^3*dtetaX_dx^3*dtetaZ_dx)/8 + (tetaZ^3*dtetaX_dx^3*dtetaY_dx)/8 -...
     (tetaY*dtetaX_dx*dtetaY_dx^2*dtetaZ_dx)/8 + (tetaZ*dtetaX_dx*dtetaY_dx*dtetaZ_dx^2)/8 +...
     (tetaY^2*tetaZ*dtetaX_dx^3*dtetaY_dx)/8 - (tetaY*tetaZ^2*dtetaX_dx^3*dtetaZ_dx)/8 -...
     (tetaY*tetaZ*dtetaX_dx^2*dtetaY_dx*dtetaZ_dx)/4) +...
    pi*G*R^2*ks*...
    (2*tetaY^2*du_dx + 2*tetaZ^2*du_dx +...
     tetaY^2*du_dx^2 + tetaZ^2*du_dx^2 -...
     2*tetaZ*du_dx*dv_dx + 2*tetaY*du_dx*dw_dx) +...
    pi*G*R^4*ks*...
    ((tetaY^2*dtetaX_dx^2)/2 + (tetaY^2*dtetaY_dx^2)/4 +...
     (tetaZ^2*dtetaX_dx^2)/2 + (tetaY^2*dtetaZ_dx^2)/4 + (tetaZ^2*dtetaY_dx^2)/4 +...
     (tetaY^4*dtetaX_dx^2)/4 + (tetaZ^2*dtetaZ_dx^2)/4 + (tetaZ^4*dtetaX_dx^2)/4 -...
     (tetaY*dtetaX_dx*dtetaZ_dx)/2 + (tetaZ*dtetaX_dx*dtetaY_dx)/2 +...
     (tetaY^2*tetaZ^2*dtetaX_dx^2)/2 - (tetaY^3*dtetaX_dx*dtetaZ_dx)/2 +...
     (tetaZ^3*dtetaX_dx*dtetaY_dx)/2 + (tetaY^2*tetaZ*dtetaX_dx*dtetaY_dx)/2 -...
     (tetaY*tetaZ^2*dtetaX_dx*dtetaZ_dx)/2));
Psi_x = @(du_dx, dv_dx, dtetaZ_dx, dw_dx, dtetaY_dx, dtetaX_dx, tetaZ, tetaY)(...
    Psi_x_linear(du_dx, dv_dx, dtetaZ_dx, dw_dx, dtetaY_dx, dtetaX_dx, tetaZ, tetaY) +...
    Psi_x_nlinear(du_dx, dv_dx, dtetaZ_dx, dw_dx, dtetaY_dx, dtetaX_dx, tetaZ, tetaY));
repositorio.U_Psi = Psi_x;

%%
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
    (l2*((NtetaX(xi, le).')*(Nw(xi, le)*cos(alfa)) +...
           (Nw(xi, le).'*cos(alfa))*(NtetaX(xi, le)))) +...
    dm*...
    (l2*(((Nu(xi, le).')*(NtetaY(xi, le)*sin(alfa) -...
                    NtetaZ(xi, le)*cos(alfa))) +...
           ((NtetaY(xi, le).'*sin(alfa) - NtetaZ(xi, le).'*cos(alfa))*(Nu(xi, le))))) +...
    dm*...
    (l2^2*((NtetaY(xi, le).'*sin(alfa) - NtetaZ(xi, le).'*cos(alfa))*(NtetaY(xi, le)*sin(alfa) - NtetaZ(xi, le)*cos(alfa))) -...
     l2*(Nv(xi, le).'*(NtetaX(xi, le)*sin(alfa)) +...
        (NtetaX(xi, le).'*sin(alfa))*Nv(xi, le))));
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
    (2*l2*(((NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*(((tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa))*Nu(xi, le)) +...
                                         (cos(tetaX(xi, le, uj) + alfa)*Nw(xi, le)))) +...
           ((Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfa) + Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)))*(NtetaX(xi, le) + (NtetaY(xi, le)*tetaZ(xi, le, uj)))))) +...
    dm*...
    (2*l2*(((Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*((NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa)) -...
                                                             (NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa)))) +...
           ((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))*(Nu(xi, le) - (Nw(xi, le)*tetaY(xi, le, uj)))))) +...
    dm*...
    (l2^2*2*((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa))) -...
     2*(Nv(xi, le).'*((Nw(xi, le)*tetaZ(xi, le, uj)*tetaY(xi, le, uj)) +...
              (tetaZ(xi, le, uj)*NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa))*l2 +...
              (NtetaX(xi, le)*sin(tetaX(xi, le, uj) + alfa))*l2) +...
        ((Nw(xi, le).'*tetaZ(xi, le, uj)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*l2*NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa) + l2*NtetaX(xi, le).'*sin(tetaX(xi, le, uj) + alfa))*Nv(xi, le)))) -...
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
    (2*l2*(((NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((upp(xi, le, uppj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa))) +...
                                                     (wpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa)))) +...
           ((tetaXpp(xi, le, uppj) + (tetaYpp(xi, le, uppj)*tetaZ(xi, le, uj)))*(Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfa) +...
                                                                                 Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) +...
                                                                                               tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)))))) +...
    dm*...
    (2*l2*(((Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*((tetaYpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfa)) -...
                                                             (tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa)))) +...
           ((upp(xi, le, uppj) - (wpp(xi, le, uppj)*tetaY(xi, le, uj)))*(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) -...
                                                                         NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))))) +...
    dm*...
    (l2^2*2*(((tetaYpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfa)) - (tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa)))*(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))) -...
     2*(Nv(xi, le).'*((wpp(xi, le, uppj)*tetaZ(xi, le, uj)*tetaY(xi, le, uj)) +...
                      (tetaZ(xi, le, uj)*tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa))*l2 +...
                      (tetaXpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfa))*l2) +...
        (vpp(xi, le, uppj)*(Nw(xi, le).'*tetaZ(xi, le, uj)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*l2*NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa) + l2*NtetaX(xi, le).'*sin(tetaX(xi, le, uj) + alfa))))) -...
    T_dmj_dm_Aux(xi, le)*uppj);
%repositorio.T_m2juppj_m2_Aux = T_m2juppj_m2_Aux;
%ok! 22/12/2021
T_dm2juppjdujT_m2_Aux = @(xi, le, uj, uppj)(dm*...
    (l2^2*(2*(NtetaX(xi, le).'*tetaYpp(xi, le, uppj) + tetaXpp(xi, le, uppj)*NtetaY(xi, le).')*NtetaZ(xi, le) +...
           2*tetaYpp(xi, le, uppj)*NtetaY(xi, le).'*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)) +...
     ((2*vpp(xi, le, uppj)*Nv(xi, le).' + 2*upp(xi, le, uppj)*Nu(xi, le).')*2*tetaZ(xi, le, uj)*NtetaZ(xi, le)) +...
     ((2*wpp(xi, le, uppj)*Nw(xi, le).' + 2*upp(xi, le, uppj)*Nu(xi, le).')*2*tetaY(xi, le, uj)*NtetaY(xi, le))) +...
    dm*...
    (2*l2*((((NtetaY(xi, le).'*NtetaZ(xi, le))*((upp(xi, le, uppj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa))) +...
                                                (wpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa))) +...
             (NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((upp(xi, le, uppj)*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) + (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)))) +...
                                                                      (-wpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))))) +...
           ((Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfa) + Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)))*(tetaYpp(xi, le, uppj)*NtetaZ(xi, le)) +...
            (tetaXpp(xi, le, uppj) + (tetaYpp(xi, le, uppj)*tetaZ(xi, le, uj)))*(-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + Nu(xi, le).'*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) + (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))))))) +...
    dm*...
    (2*l2*(((-Nw(xi, le).'*NtetaY(xi, le))*((tetaYpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfa)) -...
                                            (tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa))) +...
            (Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*(tetaYpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) +...
                                                             tetaZpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))) +...
           ((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))*(-wpp(xi, le, uppj)*NtetaY(xi, le)) + (upp(xi, le, uppj) - (wpp(xi, le, uppj)*tetaY(xi, le, uj)))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))))) +...
    dm*...
    (l2^2*2*((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))*(tetaYpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + tetaZpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
             (tetaYpp(xi, le, uppj)*sin(tetaX(xi, le, uj) + alfa) - tetaZpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))) -...
     2*(Nv(xi, le).'*(wpp(xi, le, uppj)*(NtetaZ(xi, le)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le)) +...
                      tetaZpp(xi, le, uppj)*(NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*l2 +...
              tetaXpp(xi, le, uppj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*l2) +...
        (vpp(xi, le, uppj)*(Nw(xi, le).'*(NtetaZ(xi, le)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le)) +...
              NtetaZ(xi, le).'*(NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*l2 +...
              NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*l2)))));
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
    (2*l2*(((NtetaY(xi, le).'*tetaZp(xi, le, upj))*(up(xi, le, upj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)) + wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)) +...
            (NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((up(xi, le, upj)*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
                                              (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)))) +...
                                         (-wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)))) +...
           (((tetaYp(xi, le, upj)*tetaZp(xi, le, upj)))*(Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfa) + Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa))) +...
            (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) +...
                                     Nu(xi, le).'*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
                                           (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))))))) +...
    dm*...
    (2*l2*(((-Nw(xi, le).'*tetaYp(xi, le, upj))*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)) +...
            (Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*((tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) -...
                                 (-tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)))) +...
           ((-(wp(xi, le, upj)*tetaYp(xi, le, upj)))*(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa)) +...
            (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))))) +...
    dm*...
    (l2^2*2*(((tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) - (-tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)))*(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa)) +...
             (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))) -...
     2*(Nv(xi, le).'*((wp(xi, le, upj)*tetaZp(xi, le, upj)*tetaY(xi, le, uj) + wp(xi, le, upj)*tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) +...
              (tetaZp(xi, le, upj)^2*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))*l2 +...
              (tetaXp(xi, le, upj)^2*cos(tetaX(xi, le, uj) + alfa))*l2) +...
        (vp(xi, le, upj)*(Nw(xi, le).'*(tetaZp(xi, le, upj)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) +...
             l2*NtetaZ(xi, le).'*(tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
             l2*NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))))) -...
    dm*(l2^2*(2*tetaXp(xi, le, upj)*tetaYp(xi, le, upj) + tetaYp(xi, le, upj)^2*2*tetaZ(xi, le, uj))*NtetaZ(xi, le).' +...
    2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*(vp(xi, le, upj)^2 + up(xi, le, upj)^2) +...
    2*tetaY(xi, le, uj)*NtetaY(xi, le).'*(wp(xi, le, upj)^2 + up(xi, le, upj)^2) +...
    2*l2*((tetaYp(xi, le, upj)*NtetaZ(xi, le).')*(up(xi, le, upj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)) + wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)) +...
          (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(up(xi, le, upj)*((NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).') +...
                                       (NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')) - wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')) +...
    2*l2*((-wp(xi, le, upj)*NtetaY(xi, le).')*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)) +...
          (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa))*NtetaX(xi, le).') +...
    l2^2*2*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa))*NtetaX(xi, le).' -...
    2*vp(xi, le, upj)*(wp(xi, le, upj)*(NtetaZ(xi, le).'*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le).') +...
    l2*tetaZp(xi, le, upj)*(NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).') +...
    l2*tetaXp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')));
%repositorio.T_gCj_m2_Aux = T_gCj_m2_Aux;
T_dgCjdujT_m2_Aux = @(xi, le, uj, upj)(dm*...
    (l2^2*(2*(2*tetaYp(xi, le, upj)*tetaZp(xi, le, upj))*(NtetaY(xi, le).'*NtetaZ(xi, le))) +...
     (2*tetaZp(xi, le, upj)*(2*vp(xi, le, upj)*Nv(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')*NtetaZ(xi, le)) +...
     (2*tetaYp(xi, le, upj)*(2*wp(xi, le, upj)*Nw(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')*NtetaY(xi, le))) +...
    dm*...
    (2*l2*(((NtetaY(xi, le).'*tetaZp(xi, le, upj))*(up(xi, le, upj)*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
                                   (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))) -...
                               wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
            ((NtetaY(xi, le).'*NtetaZ(xi, le))*((up(xi, le, upj)*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
                                              (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)))) +...
                                (-wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))) +...
             (NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((up(xi, le, upj)*((tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + (NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*tetaXp(xi, le, upj)) +...
                                               (-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) - (NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*tetaXp(xi, le, upj)))) +...
                                          (-wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj))))) +...
           ((tetaYp(xi, le, upj)*tetaZp(xi, le, upj))*(-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) +...
                             Nu(xi, le).'*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
                                   (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)))) +...
            ((-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) +...
              Nu(xi, le).'*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
                    (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))))*(tetaYp(xi, le, upj)*NtetaZ(xi, le)) +...
             (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(-Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj) +...
                                      Nu(xi, le).'*((tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + (NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*tetaXp(xi, le, upj)) +...
                                            (-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) - (NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*tetaXp(xi, le, upj)))))))) +...
    dm*...
    (2*l2*(((-Nw(xi, le).'*tetaYp(xi, le, upj))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
            ((-Nw(xi, le).'*NtetaY(xi, le))*((tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
                             (tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))) +...
             (Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*((-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj)) -...
                                  (-tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj))))) +...
           ((-wp(xi, le, upj)*tetaYp(xi, le, upj))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
            ((NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))*(-wp(xi, le, upj)*NtetaY(xi, le)) +...
             (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(-NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj)))))) +...
    dm*...
    (l2^2*2*(((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))*(-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj)) +...
              (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))) +...
             ((NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
              (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa))*(-NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj)))) -...
     2*(Nv(xi, le).'*((wp(xi, le, upj)*tetaZp(xi, le, upj)*NtetaY(xi, le) + wp(xi, le, upj)*NtetaZ(xi, le)*tetaYp(xi, le, upj)) +...
              (-tetaZp(xi, le, upj)^2*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) - (NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*tetaZp(xi, le, upj)*tetaXp(xi, le, upj))*l2 +...
              (-tetaXp(xi, le, upj)^2*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*l2) +...
        vp(xi, le, upj)*(Nw(xi, le).'*(tetaZp(xi, le, upj)*NtetaY(xi, le) + NtetaZ(xi, le)*tetaYp(xi, le, upj)) +...
            l2*NtetaZ(xi, le).'*(-tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) - (NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*tetaXp(xi, le, upj)) -...
            l2*NtetaX(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)*tetaXp(xi, le, upj)))) -...
    dm*...
    (l2^2*NtetaZ(xi, le).'*(tetaYp(xi, le, upj)^2*2*NtetaZ(xi, le)) +...
     2*NtetaZ(xi, le).'*(vp(xi, le, upj)^2 + up(xi, le, upj)^2)*NtetaZ(xi, le) +...
     2*NtetaY(xi, le).'*(wp(xi, le, upj)^2 + up(xi, le, upj)^2)*NtetaY(xi, le) +...
     2*l2*((tetaYp(xi, le, upj)*NtetaZ(xi, le).')*(up(xi, le, upj)*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
                                  (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))) -...
                              wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
           ((up(xi, le, upj)*((NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) +...
                  NtetaX(xi, le).'*(NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))) +...
                 (-NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) -...
                  NtetaX(xi, le).'*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)))) -...
             wp(xi, le, upj)*NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*(tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj)) +...
            (up(xi, le, upj)*((NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).') +...
                 (NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')) -...
             wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')*(tetaYp(xi, le, upj)*NtetaZ(xi, le)))) +...
     2*l2*((-wp(xi, le, upj)*NtetaY(xi, le).')*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
           NtetaX(xi, le).'*((-wp(xi, le, upj)*NtetaY(xi, le))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)) +...
                     (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)))) +...
     l2^2*2*NtetaX(xi, le).'*((tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)) +...
                      (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa))*(-tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))) -...
     2*vp(xi, le, upj)*(wp(xi, le, upj)*(NtetaZ(xi, le).'*NtetaY(xi, le) + NtetaY(xi, le).'*NtetaZ(xi, le)) +...
           l2*tetaZp(xi, le, upj)*(-NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) - NtetaX(xi, le).'*(NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))) -...
           l2*tetaXp(xi, le, upj)*NtetaX(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))));
%repositorio.T_dgCj_dujT_m2_Aux = T_dgCj_dujT_m2_Aux;
T_dgCjdupjT_m2_Aux = @(xi, le, uj, upj)(dm*...
    (l2^2*(2*(((NtetaX(xi, le).'*tetaYp(xi, le, upj) + tetaXp(xi, le, upj)*NtetaY(xi, le).')*NtetaZ(xi, le) +...
               tetaZp(xi, le, upj)*(NtetaX(xi, le).'*NtetaY(xi, le) + NtetaY(xi, le).'*NtetaX(xi, le)))) +...
           2*(2*tetaZ(xi, le, uj))*NtetaY(xi, le).'*(NtetaY(xi, le)*tetaZp(xi, le, upj) + tetaYp(xi, le, upj)*NtetaZ(xi, le))) +...
     (2*tetaZ(xi, le, uj)*((2*vp(xi, le, upj)*Nv(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')*NtetaZ(xi, le) + 2*tetaZp(xi, le, upj)*(Nv(xi, le).'*Nv(xi, le) + Nu(xi, le).'*Nu(xi, le)))) +...
     (2*tetaY(xi, le, uj)*((2*wp(xi, le, upj)*Nw(xi, le).' + 2*up(xi, le, upj)*Nu(xi, le).')*NtetaY(xi, le) + 2*tetaYp(xi, le, upj)*(Nw(xi, le).'*Nw(xi, le) + Nu(xi, le).'*Nu(xi, le))))) +...
    dm*...
    (2*l2*((((NtetaY(xi, le).'*NtetaZ(xi, le))*(up(xi, le, upj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)) + wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)) +...
             (NtetaY(xi, le).'*tetaZp(xi, le, upj))*(Nu(xi, le)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)) + Nw(xi, le)*cos(tetaX(xi, le, uj) + alfa))) +...
            (NtetaX(xi, le).' + NtetaY(xi, le).'*tetaZ(xi, le, uj))*((Nu(xi, le)*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
                                              (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))) +...
                                          up(xi, le, upj)*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
                                              (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)))) +...
                                         (-sin(tetaX(xi, le, uj) + alfa)*(Nw(xi, le)*tetaXp(xi, le, upj) + wp(xi, le, upj)*NtetaX(xi, le))))) +...
           ((Nw(xi, le).'*cos(tetaX(xi, le, uj) + alfa) + Nu(xi, le).'*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)))*(NtetaY(xi, le)*tetaZp(xi, le, upj) + tetaYp(xi, le, upj)*NtetaZ(xi, le)) +...
            ((-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) +...
              Nu(xi, le).'*((tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
                    (tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))))*(NtetaX(xi, le) + NtetaY(xi, le)*tetaZ(xi, le, uj)) +...
             (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(-Nw(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) +...
                                      Nu(xi, le).'*((NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
                                            (NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)))))))) +...
    dm*...
    (2*l2*((-Nw(xi, le).'*(NtetaY(xi, le)*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)) +...
                   tetaYp(xi, le, upj)*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa))) +...
            (Nu(xi, le).' - Nw(xi, le).'*tetaY(xi, le, uj))*(cos(tetaX(xi, le, uj) + alfa)*(NtetaY(xi, le)*tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*NtetaX(xi, le)) +...
                                 sin(tetaX(xi, le, uj) + alfa)*(NtetaZ(xi, le)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*NtetaX(xi, le)))) +...
           (-(NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))*(Nw(xi, le)*tetaYp(xi, le, upj) + wp(xi, le, upj)*NtetaY(xi, le)) +...
            ((NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))*(Nu(xi, le) - Nw(xi, le)*tetaY(xi, le, uj)) +...
             (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)))))) +...
    dm*...
    (l2^2*2*((NtetaY(xi, le).'*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa))*(cos(tetaX(xi, le, uj) + alfa)*(NtetaY(xi, le)*tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*NtetaX(xi, le)) + sin(tetaX(xi, le, uj) + alfa)*(NtetaZ(xi, le)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*NtetaX(xi, le))) +...
             ((NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa)) +...
              (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa))*(NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le) + NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)))) -...
     2*(Nv(xi, le).'*(((Nw(xi, le)*tetaZp(xi, le, upj) + wp(xi, le, upj)*NtetaZ(xi, le))*tetaY(xi, le, uj) + (Nw(xi, le)*tetaYp(xi, le, upj) + wp(xi, le, upj)*NtetaY(xi, le))*tetaZ(xi, le, uj)) +...
              (2*tetaZp(xi, le, upj)*NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*(NtetaZ(xi, le)*tetaXp(xi, le, upj) + tetaZp(xi, le, upj)*NtetaX(xi, le)))*l2 +...
              (2*tetaXp(xi, le, upj)*NtetaX(xi, le)*cos(tetaX(xi, le, uj) + alfa))*l2) +...
        ((Nw(xi, le).'*(tetaZp(xi, le, upj)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*tetaYp(xi, le, upj)) +...
          l2*NtetaZ(xi, le).'*(tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj)) +...
          l2*NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*tetaXp(xi, le, upj))*Nv(xi, le) +...
         vp(xi, le, upj)*(Nw(xi, le).'*(NtetaZ(xi, le)*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le)) +...
             l2*NtetaZ(xi, le).'*(NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le)) +...
             l2*NtetaX(xi, le).'*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le))))) -...
    dm*...
    (2*l2^2*NtetaZ(xi, le).'*((NtetaX(xi, le)*tetaYp(xi, le, upj) + tetaXp(xi, le, upj)*NtetaY(xi, le)) + 2*tetaYp(xi, le, upj)*NtetaY(xi, le)*tetaZ(xi, le, uj)) +...
     2*tetaZ(xi, le, uj)*NtetaZ(xi, le).'*(2*vp(xi, le, upj)*Nv(xi, le) + 2*up(xi, le, upj)*Nu(xi, le)) +...
     2*tetaY(xi, le, uj)*NtetaY(xi, le).'*(2*wp(xi, le, upj)*Nw(xi, le) + 2*up(xi, le, upj)*Nu(xi, le)) +...
     2*l2*(((NtetaZ(xi, le).'*NtetaY(xi, le))*(up(xi, le, upj)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)) + wp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)) +...
            (tetaYp(xi, le, upj)*NtetaZ(xi, le).')*(Nu(xi, le)*(tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa) + tetaY(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)) + Nw(xi, le)*cos(tetaX(xi, le, uj) + alfa))) +...
           ((up(xi, le, upj)*((NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).') +...
                 (NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')) -...
             wp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')*(NtetaX(xi, le) + NtetaY(xi, le)*tetaZ(xi, le, uj)) +...
            (tetaXp(xi, le, upj) + tetaYp(xi, le, upj)*tetaZ(xi, le, uj))*(((NtetaZ(xi, le).'*sin(tetaX(xi, le, uj) + alfa) + tetaZ(xi, le, uj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).') +...
                                         (NtetaY(xi, le).'*cos(tetaX(xi, le, uj) + alfa) - tetaY(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).'))*Nu(xi, le) -...
                                     sin(tetaX(xi, le, uj) + alfa)*(NtetaX(xi, le).'*Nw(xi, le))))) +...
     2*l2*(((-NtetaY(xi, le).'*Nw(xi, le))*(tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)) +...
            (-wp(xi, le, upj)*NtetaY(xi, le).')*(NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa))) +...
           NtetaX(xi, le).'*((Nu(xi, le) - Nw(xi, le)*tetaY(xi, le, uj))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)) +...
                     (up(xi, le, upj) - wp(xi, le, upj)*tetaY(xi, le, uj))*(NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) + NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa)))) +...
     l2^2*2*NtetaX(xi, le).'*((NtetaY(xi, le)*sin(tetaX(xi, le, uj) + alfa) - NtetaZ(xi, le)*cos(tetaX(xi, le, uj) + alfa))*(tetaYp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa) + tetaZp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa)) +...
                      (tetaYp(xi, le, upj)*sin(tetaX(xi, le, uj) + alfa) - tetaZp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa))*(NtetaY(xi, le)*cos(tetaX(xi, le, uj) + alfa) + NtetaZ(xi, le)*sin(tetaX(xi, le, uj) + alfa))) -...
     2*((wp(xi, le, upj)*(NtetaZ(xi, le).'*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le).') +...
         l2*tetaZp(xi, le, upj)*(NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).') +...
         l2*tetaXp(xi, le, upj)*cos(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')*Nv(xi, le) +...
        vp(xi, le, upj)*((NtetaZ(xi, le).'*tetaY(xi, le, uj) + tetaZ(xi, le, uj)*NtetaY(xi, le).')*Nw(xi, le) +...
            l2*(NtetaZ(xi, le).'*cos(tetaX(xi, le, uj) + alfa) - tetaZ(xi, le, uj)*sin(tetaX(xi, le, uj) + alfa)*NtetaX(xi, le).')*NtetaZ(xi, le) +...
            l2*cos(tetaX(xi, le, uj) + alfa)*(NtetaX(xi, le).'*NtetaX(xi, le))))));
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
    (1/2)*(T_m1j_Aux(1, le, M2b, Ip2, I2) + T_m1j_m2_Aux(1, le)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_m1j_Aux(0, le, M2b, Ip2, I2) + T_m1j_m2_Aux(0, le)) +...
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
    (1/2)*(T_m2j_Aux(1, le, uj, Ip2) + T_m2j_m2_Aux(1, le, uj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_m2j_Aux(0, le, uj, Ip2) + T_m2j_m2_Aux(0, le, uj)) +...
    ad_N(nElem)*...
    (1/2)*T_m2j_Aux(1, le, uj, Ip1));
repositorio.T_M2j_c = T_M2j_c;

%(2) Matriz M2juppj_c:
T_M2juppj_c = @(nElem, le, uj, uppj)(ad_1(nElem)*...
    (1/2)*T_m2uppj_M3_Aux(0, le, uj, uppj) +...
    ad_N1(nElem)*...
    (1/2)*(T_m2juppj_Aux(1, le, uj, uppj, Ip2) + T_m2juppj_m2_Aux(1, le, uj, uppj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_m2juppj_Aux(0, le, uj, uppj, Ip2) + T_m2juppj_m2_Aux(0, le, uj, uppj)) +...
    ad_N(nElem)*...
    (1/2)*T_m2juppj_Aux(1, le, uj, uppj, Ip1));
repositorio.T_M2juppj_c = T_M2juppj_c;

%(3) Matriz dM2juppj_dujT_c:
T_dM2juppjdujT_c = @(nElem, le, uj, uppj)(ad_1(nElem)*...
    (1/2)*T_dm2juppjdujT_M3_Aux(0, le, uj, uppj) +...
    ad_N1(nElem)*...
    (1/2)*(T_dm2juppj_dujT_Aux(1, le, uj, uppj, Ip2) + T_dm2juppjdujT_m2_Aux(1, le, uj, uppj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_dm2juppj_dujT_Aux(0, le, uj, uppj, Ip2) + T_dm2juppjdujT_m2_Aux(0, le, uj, uppj)) +...
    ad_N(nElem)*...
    (1/2)*T_dm2juppj_dujT_Aux(1, le, uj, uppj, Ip1));
repositorio.T_dM2juppj_dujT_c = T_dM2juppjdujT_c;

%(7) Matriz Gj_c:
T_Gj_c = @(nElem, le, uj, upj)(ad_1(nElem)*...
    (1/2)*T_gCj_M3_Aux(0, le, uj, upj) +...
    ad_N1(nElem)*...
    (1/2)*(T_gCj_Aux(1, le, uj, upj, Ip2) + T_gCj_m2_Aux(1, le, uj, upj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_gCj_Aux(0, le, uj, upj, Ip2) + T_gCj_m2_Aux(0, le, uj, upj)) +...
    ad_N(nElem)*...
    (1/2)*T_gCj_Aux(1, le, uj, upj, Ip1));
repositorio.T_Gj_c = T_Gj_c;

%(8) Matriz dGj_dujT_c:
T_dGjdujT_c = @(nElem, le, uj, upj)(ad_1(nElem)*...
    (1/2)*T_dgCjdujT_M3_Aux(0, le, uj, upj) +...
    ad_N1(nElem)*...
    (1/2)*(T_dgCjdujT_Aux(1, le, upj, Ip2) + T_dgCjdujT_m2_Aux(1, le, uj, upj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_dgCjdujT_Aux(0, le, upj, Ip2) + T_dgCjdujT_m2_Aux(0, le, uj, upj)) +...
    ad_N(nElem)*...
    (1/2)*T_dgCjdujT_Aux(1, le, upj, Ip1));
repositorio.T_dGj_dujT_c = T_dGjdujT_c;

%(9) Matriz dGj_dupjT_c:
T_dGjdupjT_c = @(nElem, le, uj, upj)(ad_1(nElem)*...
    (1/2)*T_dgCjdupjT_M3_Aux(0, le, uj, upj) +...
    ad_N1(nElem)*...
    (1/2)*(T_dgCjdupjT_Aux(1, le, uj, upj, Ip2) + T_dgCjdupjT_m2_Aux(1, le, uj, upj)) +...
    ad_N1M1(nElem)*...
    (1/2)*(T_dgCjdupjT_Aux(0, le, uj, upj, Ip2) + T_dgCjdupjT_m2_Aux(0, le, uj, upj)) +...
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
    (cjbU_Aux(0, le, bUl) + cjbTx_Aux(0, le, bTx1)) +...
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
%Newmark:
%dt = argumentos.dt_newmark;
gama = argumentos.gama_newmark;
beta = argumentos.beta_newmark;
f_up_nM1 = @(u_nM1, u_n, up_n, upp_n, dt)((gama/(beta*dt))*(u_nM1 - u_n) +...
    (1 - gama/beta)*up_n + (1 - gama/(2*beta))*dt*upp_n);
repositorio.f_up_nM1 = f_up_nM1;
f_upp_nM1 = @(u_nM1, u_n, up_n, upp_n, dt)((1/(beta*dt^2))*(u_nM1 - u_n) -...
    (1/(beta*dt))*up_n + (1 - 1/(2*beta))*upp_n);
repositorio.f_upp_nM1 = f_upp_nM1;

%%
%Matrizes de inércia e rigidez:
[ M1gEF, KgEF ] = fGm1k( argumentos, repositorio );
argumentos.M1gEF = M1gEF;
argumentos.KgEF = KgEF;

%%
%Ponto de equilíbrio no modelo em EF:
Ng = 6*(N + 1);
argumentos.Ng = Ng;
IdNg = eye(Ng);
argumentos.IdNg = IdNg;
%
SigmaRMaU = zeros(Ng,N+1);
for i = 1:N+1
    n = 6*i-5;
    SigmaRMaU(:,i) = IdNg(:,n);
end
SigmaRMauEq = SigmaRMaU(:,2:N+1);
%}
%{
SigmaRMaU = IdNg(:,[(6*2-5):6*(N+1)-5,6*(N+1)]);
SigmaRMauEq = SigmaRMaU;
%}
KrM1 = (SigmaRMauEq.'*KgEF)*SigmaRMauEq;
fcEFeq = @(u)(fGfcEFeq( u, argumentos, repositorio ));
fcEFrM1eq = @(urMa)(SigmaRMauEq.'*fcEFeq(SigmaRMauEq*urMa));
%Derivadas:
dfcEFeqduT = @(u)(fGdfcEFduEq( u, argumentos, repositorio ));
dfcEFrM1eqdurM1T = @(urM1)(SigmaRMauEq.'*dfcEFeqduT(SigmaRMauEq*urM1)*SigmaRMauEq);

SigmaRMaV = zeros(Ng,N+1);
SigmaRMaW = zeros(Ng,N+1);
SigmaRMaTetaX = zeros(Ng,N+1);
for i = 1:N+1
    nV = 6*i-4;
    SigmaRMaV(:,i) = IdNg(:,nV);
    nW = 6*i-2;
    SigmaRMaW(:,i) = IdNg(:,nW);
    nTetaX = 6*i;
    SigmaRMaTetaX(:,i) = IdNg(:,nTetaX);
end
argumentos.SigmaRMaV = SigmaRMaV;
argumentos.SigmaRMaW = SigmaRMaW;
argumentos.SigmaRMaTetaX = SigmaRMaTetaX;

%Ponto de equilíbrio:
%Eta0 = zeros(Ng-(rEstr+1),1);
Eta0 = SigmaRMauEq.'*zeros(Ng,1);
Psi = @(eta)(KrM1*eta + fcEFrM1eq(eta));
dPsidEtaT = @(eta)(KrM1 + dfcEFrM1eqdurM1T(eta));
[ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
uTilrM1eq = EtaConv;
%{
eta0 = uTilrM1eq;
myfun = @(eta)(max(abs(Psi(eta))));
[uTilrM1eq] = fminsearch(myfun,eta0,optimset('TolFun',1e-12));
%Obs.: Processo muito demorado. Sem resultados satisfatórios.
%}
uTilEq = SigmaRMauEq*uTilrM1eq;
%
uTilrM1eqU1 = [0; uTilrM1eq];
uTilEqEF = SigmaRMaU*uTilrM1eqU1;
%}
%{
uTilEqEF = uTilEq;
%}
argumentos.uTilEq = uTilEq;

%%
%Modelo de atrito:
%(1) Forças de interaçao no contato de atrito:
uNM1Eq = uTilEq(6*(N+1)-5);
uC = uNM1Eq;
argumentos.uCef = uC;
kNM1 = argumentos.kNM1;
%kNM1_Hat = argumentos.kNM1_Hat;
cNM1 = argumentos.c_NM1;
cNM1Hat = argumentos.c_NM1_Hat;
R_NM1 = argumentos.R_NM1;
%R_NM1_Hat = argumentos.R_NM1_Hat;

%(1.a)Condição de contato:
IipNM1 = @(u_NM1)(double(u_NM1 >= uC*ones(size(u_NM1))));
repositorio.IipNM1 = IipNM1;

%(1.b)Força normal:
FnNM1 = @(u_NM1)(IipNM1(u_NM1).*(kNM1*(u_NM1 - uC*ones(size(u_NM1)))));
repositorio.FnNM1 = FnNM1;
FnbNM1 = @(ub_NM1)(FnNM1(ub_NM1 + uNM1Eq*ones(size(ub_NM1))));
repositorio.FnbNM1 = FnbNM1;
%F_N_NM1_Hat = @(u_NM1,up_NM1)(I_ip_NM1(u_NM1)*(kNM1_Hat*(u_NM1 - uC) + c_NM1_Hat*up_NM1));
%repositorio.F_N_NM1_Hat = F_N_NM1_Hat;
dFnNM1duNM1 = @(u_NM1)(IipNM1(u_NM1)*(kNM1));
repositorio.dFnNM1duNM1 = dFnNM1duNM1;
%dF_N_NM1_duNM1_Hat = @(u_NM1,up_NM1)(I_ip_NM1(u_NM1)*(kNM1_Hat));
%repositorio.dF_N_NM1_duNM1_Hat = dF_N_NM1_duNM1_Hat;
dFnNM1dupNM1 = @(u_NM1,up_NM1)(IipNM1(u_NM1)*(cNM1));
repositorio.dFnNM1dupNM1 = dFnNM1dupNM1;
dFnNM1dupNM1hat = @(u_NM1,up_NM1)(IipNM1(u_NM1)*(cNM1Hat));
repositorio.dFnNM1dupNM1hat = dFnNM1dupNM1hat;

%(2) Modelo de atrito descontínuo na origem (nó N+1)
%(2.a) Coeficiente de atrito:
mi0 = argumentos.mi0;
%mi0_descont_Hat = argumentos.mi0_descont_Hat;
mi1 = argumentos.mi1;
%mi1_descont_Hat = argumentos.mi10_descont_Hat;
mi2 = argumentos.mi2;
%mi2_descont_Hat = argumentos.mi20_descont_Hat;
beta1 = argumentos.beta1mi;
%beta1_Hat = argumentos.beta1_mi_descont_Hat;
beta2 = argumentos.beta2mi;
%beta2_Hat = argumentos.beta2_mi_descont_Hat;
miDc = @(tetaXp)((mi0 +...
    mi1*(1 - 2/(1 + exp(beta1*abs(tetaXp)))) -...
    mi2*(1 - 2/(1 + exp(beta2*abs(tetaXp)))))*...
    sign(tetaXp));
repositorio.miDc = miDc;
dmidtXpDc = @(tetaXp)(mi1*...
    (2*beta1*exp(beta1*abs(tetaXp))/(1 + exp(beta1*abs(tetaXp)))^2) -...
    mi2*...
    (2*beta2*exp(beta2*abs(tetaXp))/(1 + exp(beta2*abs(tetaXp)))^2));
repositorio.dmidtXpDc = dmidtXpDc;
dmidtXpDcHat = @(tetaXp)(mi1_descont_Hat*...
    (2*beta1_Hat*exp(beta1_Hat*abs(tetaXp))/(1 + exp(beta1_Hat*abs(tetaXp)))^2) -...
    mi2_descont_Hat*...
    (2*beta2_Hat*exp(beta2_Hat*abs(tetaXp))/(1 + exp(beta2_Hat*abs(tetaXp)))^2));
repositorio.dmidtXpDcHat = dmidtXpDcHat;
d2midtXp2Dc = @(tetaXp)((mi1*...
    2*beta1^2*(exp(beta1*abs(tetaXp))/(1 + exp(beta1*abs(tetaXp)))^2)*...
    (1 - 2*exp(beta1*abs(tetaXp))/(1 + exp(beta1*abs(tetaXp)))) -...
    mi2*...
    2*beta2^2*(exp(beta2*abs(tetaXp))/(1 + exp(beta2*abs(tetaXp)))^2)*...
    (1 - 2*exp(beta2*abs(tetaXp))/(1 + exp(beta2*abs(tetaXp)))))*...
    sign(tetaXp));
repositorio.d2midtXp2Dc = d2midtXp2Dc;
d2midtXp2DcHat = @(tetaXp)((mi1_descont_Hat*...
    2*beta1_Hat^2*(exp(beta1_Hat*abs(tetaXp))/(1 + exp(beta1_Hat*abs(tetaXp)))^2)*...
    (1 - 2*exp(beta1_Hat*abs(tetaXp))/(1 + exp(beta1_Hat*abs(tetaXp)))) -...
    mi2_descont_Hat*...
    2*beta2_Hat^2*(exp(beta2_Hat*abs(tetaXp))/(1 + exp(beta2_Hat*abs(tetaXp)))^2)*...
    (1 - 2*exp(beta2_Hat*abs(tetaXp))/(1 + exp(beta2_Hat*abs(tetaXp)))))*...
    sign(tetaXp));
repositorio.d2midtXp2DcHat = d2midtXp2DcHat;
d3midtXp3Dc = @(tetaXp)(mi1*...
    (2*beta1^3*exp(beta1*abs(tetaXp))/(1 + exp(beta1*abs(tetaXp)))^2)*...
    ((1 - 2*exp(beta1*abs(tetaXp))/(1 + exp(beta1*abs(tetaXp))))^2 -...
    (2*exp(beta1*abs(tetaXp))/(1 + exp(beta1*abs(tetaXp))))*...
    (1 - exp(beta1*abs(tetaXp))/(1 + exp(beta1*abs(tetaXp))))) -...
    mi2*...
    (2*beta2^3*exp(beta2*abs(tetaXp))/(1 + exp(beta2*abs(tetaXp)))^2)*...
    ((1 - 2*exp(beta2*abs(tetaXp))/(1 + exp(beta2*abs(tetaXp))))^2 -...
    (2*exp(beta2*abs(tetaXp))/(1 + exp(beta2*abs(tetaXp))))*...
    (1 - exp(beta2*abs(tetaXp))/(1 + exp(beta2*abs(tetaXp))))));
repositorio.d3midtXp3Dc = d3midtXp3Dc;
d3midtXp3DcHat = @(tetaXp)(mi1_descont_Hat*...
    (2*beta1_Hat^3*exp(beta1_Hat*abs(tetaXp))/(1 + exp(beta1_Hat*abs(tetaXp)))^2)*...
    ((1 - 2*exp(beta1_Hat*abs(tetaXp))/(1 + exp(beta1_Hat*abs(tetaXp))))^2 -...
    (2*exp(beta1_Hat*abs(tetaXp))/(1 + exp(beta1_Hat*abs(tetaXp))))*...
    (1 - exp(beta1_Hat*abs(tetaXp))/(1 + exp(beta1_Hat*abs(tetaXp))))) -...
    mi2_descont_Hat*...
    (2*beta2_Hat^3*exp(beta2_Hat*abs(tetaXp))/(1 + exp(beta2_Hat*abs(tetaXp)))^2)*...
    ((1 - 2*exp(beta2_Hat*abs(tetaXp))/(1 + exp(beta2_Hat*abs(tetaXp))))^2 -...
    (2*exp(beta2_Hat*abs(tetaXp))/(1 + exp(beta2_Hat*abs(tetaXp))))*...
    (1 - exp(beta2_Hat*abs(tetaXp))/(1 + exp(beta2_Hat*abs(tetaXp))))));
repositorio.d3midtXp3DcHat = d3midtXp3DcHat;

%(2.b) Torque de atrito:
TxNM1dc = @(nElem,tetaXp_NM1,u_NM1)(ad_N(nElem)*miDc(tetaXp_NM1)*FnNM1(u_NM1)*R_NM1);
repositorio.TxNM1dc = TxNM1dc;


TxNM1dcHat = @(tetaXp_NM1,u_NM1,up_NM1)(mi_descont_Hat(tetaXp_NM1)*F_N_NM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.TxNM1dcHat = TxNM1dcHat;
dTxduNM1dc = @(tetaXp_NM1,u_NM1,up_NM1)(miDc(tetaXp_NM1)*dFnNM1duNM1(u_NM1)*R_NM1);
repositorio.dTxduNM1dc = dTxduNM1dc;
dTxduNM1dcHat = @(tetaXp_NM1,u_NM1,up_NM1)(mi_descont_Hat(tetaXp_NM1)*dF_N_NM1_duNM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.dTxduNM1dcHat = dTxduNM1dcHat;
dTxdtXpNM1dc = @(nElem,tetaXp_NM1,u_NM1)(dmidtXpDc(tetaXp_NM1)*FnNM1(u_NM1)*R_NM1);
repositorio.dTxdtXpNM1dc = dTxdtXpNM1dc;
dTxdtXpNM1dcHat = @(nElem,tetaXp_NM1,u_NM1,up_NM1)(dmidtXpDcHat(tetaXp_NM1)*F_N_NM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.dTx_dtetaXpNM1_descont_Hat = dTxdtXpNM1dcHat;
d2TxdtXpNM1duNM1dc =...
    @(nElem,tetaXp_NM1,u_NM1,up_NM1)(dmidtXpDc(tetaXp_NM1)*dFnNM1duNM1(u_NM1)*R_NM1);
repositorio.d2TxdtXpNM1duNM1dc = d2TxdtXpNM1duNM1dc;
d2TxdtXpNM1duNM1dcHat =...
    @(nElem,tetaXp_NM1,u_NM1,up_NM1)(dmidtXpDcHat(tetaXp_NM1)*dF_N_NM1_duNM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.d2TxdtXpNM1duNM1dcHat = d2TxdtXpNM1duNM1dcHat;
d2TxdtXpNM12Dc =...
    @(nElem,tetaXp_NM1,u_NM1)(d2midtXp2Dc(tetaXp_NM1)*FnNM1(u_NM1)*R_NM1);
repositorio.d2TxdtXpNM12Dc = d2TxdtXpNM12Dc;
d2TxdtXpNM12DcHat =...
    @(nElem,tetaXp_NM1,u_NM1,up_NM1)(d2midtXp2DcHat(tetaXp_NM1)*F_N_NM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.d2TxdtXpNM12dcHat = d2TxdtXpNM12DcHat;

%(3) Modelo de atrito contínuo na origem
%(3.a) Funções adimensionais e auxiliares:
fiAdm = @(nu, lambda_mi_cont, n_mi_cont)(tanh(pi*nu) +...
    lambda_mi_cont*...
    (2*n_mi_cont/(2*n_mi_cont - 1))/...
    (1 + (1/(2*n_mi_cont - 1))*nu^(2*n_mi_cont)));
fiAdmVt = @(nu, lambda_mi_cont, n_mi_cont)(tanh(pi*nu) +...
    lambda_mi_cont*...
    (2*n_mi_cont/(2*n_mi_cont - 1))*ones(size(nu))./...
    (ones(size(nu)) + (1/(2*n_mi_cont - 1))*nu.^(2*n_mi_cont)));
dfiAdmdnu = @(nu, lambda_mi_cont, n_mi_cont)(-pi*(tanh(pi*nu)^2 - 1) -...
    (4*lambda_mi_cont*n_mi_cont^2*nu^(2*n_mi_cont - 1))/...
    ((2*n_mi_cont - 1)^2*(nu^(2*n_mi_cont)/(2*n_mi_cont - 1) + 1)^2));
d2fiAdmdnu2 = @(nu, lambda_mi_cont, n_mi_cont)(2*pi^2*tanh(pi*nu)*(tanh(pi*nu)^2 - 1) -...
    (4*lambda_mi_cont*n_mi_cont^2*nu^(2*n_mi_cont - 2))/...
    ((2*n_mi_cont - 1)*(nu^(2*n_mi_cont)/(2*n_mi_cont - 1) + 1)^2) +...
    (16*lambda_mi_cont*n_mi_cont^3*nu^(4*n_mi_cont - 2))/...
    ((2*n_mi_cont - 1)^3*(nu^(2*n_mi_cont)/(2*n_mi_cont - 1) + 1)^3));
%{
fiAdm = @(eta, alfa1, alfa2)(tanh(eta) + alfa1*eta/(1+alfa2*eta^2));
%}
miAtr = @(tetaXp, mi0, w0, lambda_mi_cont, n_mi_cont)(mi0*...
    fiAdm(tetaXp/w0, lambda_mi_cont, n_mi_cont));
miAtrVt = @(tetaXp, mi0, w0, lambda_mi_cont, n_mi_cont)(mi0*...
    fiAdmVt(tetaXp/w0, lambda_mi_cont, n_mi_cont));
dmiAtrdtXp = @(tetaXp, mi0, w0, lambda_mi_cont, n_mi_cont)((mi0/w0)*...
    dfiAdmdnu(tetaXp/w0, lambda_mi_cont, n_mi_cont));
d2miAtrdtXp2 = @(tetaXp, mi0, w0, lambda_mi_cont, n_mi_cont)((mi0/w0^2)*...
    d2fiAdmdnu2(tetaXp/w0, lambda_mi_cont, n_mi_cont));

%(3.b) Funções de atrito no ponto de contato (nó N+1):
mi0NM1 = argumentos.mic0;
w0NM1 = argumentos.w0;
lambdaNM1miC = argumentos.lambdami;
nNM1miC = argumentos.nmi;

miNM1nElemC = @(nElem, tetaXp_NM1)(ad_N(nElem)*...
    miAtr(tetaXp_NM1, mi0NM1, w0NM1, lambdaNM1miC, nNM1miC));
repositorio.miNM1nElemC = miNM1nElemC;
dmidtXpNM1nElemC = @(nElem, tetaXp_NM1)(ad_N(nElem)*...
    dmiAtrdtXp(tetaXp_NM1, mi0NM1, w0NM1, lambdaNM1miC, nNM1miC));
repositorio.dmidtXpNM1nElemC = dmidtXpNM1nElemC;
d2midtXpNM12nElemC = @(nElem, tetaXp_NM1)(ad_N(nElem)*...
    d2miAtrdtXp2(tetaXp_NM1, mi0NM1, w0NM1, lambdaNM1miC, nNM1miC));
repositorio.d2midtXpNM12nElemC = d2midtXpNM12nElemC;

%(3.c)Torque de atrito:
TcxNM1vt = @(tetaXp_NM1,ub_NM1)(miAtrVt(tetaXp_NM1, mi0NM1, w0NM1, lambdaNM1miC, nNM1miC).*...
    FnbNM1(ub_NM1)*R_NM1);
repositorio.TcxNM1vt = TcxNM1vt;

TxNM1nElem = @(nElem,tetaXp_NM1,u_NM1)(miNM1nElemC(nElem, tetaXp_NM1)*FnNM1(u_NM1)*R_NM1);
repositorio.TxNM1nElemC = TxNM1nElem;
TxNM1nElemHat = @(nElem,tetaXp_NM1,u_NM1,up_NM1)(miNM1nElemC(nElem, tetaXp_NM1)*F_N_NM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.TxNM1nElemChat = TxNM1nElemHat;

dTxduNM1nElem = @(nElem,tetaXp_NM1,u_NM1,up_NM1)(miNM1nElemC(nElem, tetaXp_NM1)*dFnNM1duNM1(u_NM1)*R_NM1);
repositorio.dTxduNM1nElemC = dTxduNM1nElem;
dTxduNM1nElemHat = @(nElem,tetaXp_NM1,u_NM1,up_NM1)(miNM1nElemC(nElem, tetaXp_NM1)*dF_N_NM1_duNM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.dTxduNM1nElemChat = dTxduNM1nElemHat;

dTxdupNM1nElem = @(nElem,tetaXp_NM1,u_NM1,up_NM1)(miNM1nElemC(nElem, tetaXp_NM1)*dFnNM1dupNM1(u_NM1,up_NM1)*R_NM1);
repositorio.dTx_dupNM1_nElem_cont = dTxdupNM1nElem;
dTxdupNM1nElemHat = @(nElem,tetaXp_NM1,u_NM1,up_NM1)(miNM1nElemC(nElem, tetaXp_NM1)*dFnNM1dupNM1hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.dTxdupNM1nElemChat = dTxdupNM1nElemHat;

dTxdtXpNM1nElem = @(nElem,tetaXp_NM1,u_NM1)(dmidtXpNM1nElemC(nElem, tetaXp_NM1)*FnNM1(u_NM1)*R_NM1);
repositorio.dTxdtXpNM1nElemC = dTxdtXpNM1nElem;
dTxdtXpNM1nElemHat = @(nElem,tetaXp_NM1,u_NM1,up_NM1)(dmidtXpNM1nElemC(nElem, tetaXp_NM1)*F_N_NM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.dTxdtXpNM1nElemChat = dTxdtXpNM1nElemHat;

d2TxdtXpNM1duNM1nElem =...
    @(nElem,tetaXp_NM1,u_NM1,up_NM1)(dmidtXpNM1nElemC(nElem, tetaXp_NM1)*dFnNM1duNM1(u_NM1)*R_NM1);
repositorio.d2TxdtetaXpNM1duNM1nElemC = d2TxdtXpNM1duNM1nElem;
d2TxdtXpNM1duNM1nElemHat =...
    @(nElem,tetaXp_NM1,u_NM1,up_NM1)(dmidtXpNM1nElemC(nElem, tetaXp_NM1)*dF_N_NM1_duNM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.d2TxdtetaXpNM1duNM1nElemChat = d2TxdtXpNM1duNM1nElemHat;

d2TxdtXpNM12nElem =...
    @(nElem,tetaXp_NM1,u_NM1)(d2midtXpNM12nElemC(nElem, tetaXp_NM1)*FnNM1(u_NM1)*R_NM1);
repositorio.d2TxdtetaXpNM12nElemC = d2TxdtXpNM12nElem;
d2TxdtXpNM12nElemHat =...
    @(nElem,tetaXp_NM1,u_NM1,up_NM1)(d2midtXpNM12nElemC(nElem, tetaXp_NM1)*F_N_NM1_Hat(u_NM1,up_NM1)*R_NM1_Hat);
repositorio.d2TxdtXpNM12nElemChat = d2TxdtXpNM12nElemHat;

%{
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
%}

%%
%Vetores, matrizes e funçoes globais:
%Equação: M1*upp + C*up + K*u + f(u,up,upp) + kNmk*Gama*lambda = 0

%Aplicação das restrições nas extremidades sem multiplicadores de Lagrange:
KgEFb = KgEF + fGdfuduT( uTilEqEF, argumentos, repositorio );
%SigmaR = [IdNg(:,1),IdNg(:,6),IdNg(:,7:6*N),IdNg(:,6*(N+1)-5),IdNg(:,6*(N+1))];
%SigmaR = [IdNg(:,7:6*N),IdNg(:,6*(N+1)-5),IdNg(:,6*(N+1))];
SigmaR = [IdNg(:,1),IdNg(:,6),IdNg(:,7:6*N1),IdNg(:,6*(N1+1)-5),IdNg(:,6*(N1+1)),IdNg(:,6*(N1+1)+1:6*N),IdNg(:,6*(N+1)-5),IdNg(:,6*(N+1))];
argumentos.SigmaR = SigmaR;
M1restr = SigmaR.'*M1gEF*SigmaR;
argumentos.M1restr = M1restr;
Krestr = SigmaR.'*KgEF*SigmaR;
argumentos.Krestr = Krestr;
KbRestr = SigmaR.'*KgEFb*SigmaR;
argumentos.KbRestr = KbRestr;
%Frequências naturais (rad/s):
freqNatKbw2 = sort(real(eig(M1restr\KbRestr)));
freqNatKbw = sqrt(freqNatKbw2); %rad/s
argumentos.freqNatKbwEF = freqNatKbw;
%{
[Pr,Lmb] = eig(M1restr\KbRestr);
[Upr,Spr,Vpr] = svd(Pr);
figure,plot(diag(Spr)),grid
%}
[ PtX, sigTiltX, Pu, sigTilu, Abvt, Abvl, info ] = fModMMod( M1restr, KbRestr, argumentos );
%Primeiros modos torsionais:
argumentos.PtX = PtX;
%Primeiras frequências torsionais:
argumentos.sigTiltX = sigTiltX;
%Primeiros modos longitudinais:
argumentos.Pu = Pu;
%Primeiras frequências longitudinais:
argumentos.sigTilu = sigTilu;
%{
%Modos intermediários entre os primeiros modos torsionais:
argumentos.Pvw = Pvw;
%Frequências intermediárias entre as primeiras frequências torsionais:
argumentos.sigTilvw = sigTilvw;
%}

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
%{
%Posição do modo Pv em Abvt:
jPv = info.jPv;
argumentos.jPv = jPv;
%Posição do modo Pw em Abvt:
jPw = info.jPw;
argumentos.jPw = jPw;
%}

%%
%Gráficos 3D dos modos de vibraçao:
%{
%Observando modos:
Xef = argumentos.X; %malha de EF
V0 = Xef.' + uTilrM1eqU1;
Auu = SigmaRMaU.';
uTilP = Auu*SigmaR*Abvt;
Avu = SigmaRMaV.';
vTilP = Avu*SigmaR*Abvt;
Awu = SigmaRMaW.';
wTilP = Awu*SigmaR*Abvt;
AtetaXu = SigmaRMaTetaX.';
tetaXTilP = AtetaXu*SigmaR*Abvt;
nAbvt = length(Abvt);
VP = zeros(N+1,3,nAbvt);
for i = 1:nAbvt
    VP(:,:,i) = [V0 + uTilP(:,i), vTilP(:,i), wTilP(:,i)];
end
graf = 0;
if graf>0
    %%
    %
    %Apenas os modos de vibração u1 e u2:
    %(a) Gráficos 3D:
    az = 90;
    el = 180;
    figure
    p = plot3(VP(:,2,jPu1),VP(:,3,jPu1),VP(:,1,jPu1),VP(:,2,jPu2),VP(:,3,jPu2),VP(:,1,jPu2));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos longitudinais - Plano XZ')
    grid
    view(az, el);
    legend('1^{o} modo','2^{o} modo');

    az = 180;
    el = 180;
    figure
    p = plot3(VP(:,2,jPu1),VP(:,3,jPu1),VP(:,1,jPu1),VP(:,2,jPu2),VP(:,3,jPu2),VP(:,1,jPu2));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos longitudinais - Plano XY')
    grid
    view(az, el);
    legend('1^{o} modo','2^{o} modo');
    
    az = 180;
    el = 90;
    figure
    p = plot3(VP(:,2,jPu1),VP(:,3,jPu1),VP(:,1,jPu1),VP(:,2,jPu2),VP(:,3,jPu2),VP(:,1,jPu2));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos longitudinais - Plano YZ')
    grid
    view(az, el);
    legend('1^{o} modo','2^{o} modo');
    
    %(b) Gráficos 2D:
    figure
    p = plot(VP(:,1,jPu1),VP(:,2,jPu1),'b-',VP(:,1,jPu2),VP(:,2,jPu2),'r-',VP(1,1,jPu1),VP(1,2,jPu1),'bo',VP(N1+1,1,jPu1),VP(N1+1,2,jPu1),'bo',VP(N+1,1,jPu1),VP(N+1,2,jPu1),'bo',VP(1,1,jPu2),VP(1,2,jPu2),'ro',VP(N1+1,1,jPu2),VP(N1+1,2,jPu2),'ro',VP(N+1,1,jPu2),VP(N+1,2,jPu2),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    camroll(-90)
    xlabel('x (m)')
    ylabel('v')
    %title('Modos longitudinais')
    grid
    legend('1^{o} modo','2^{o} modo');
    
    figure
    p = plot(VP(:,1,jPu1),VP(:,3,jPu1),'b-',VP(:,1,jPu2),VP(:,3,jPu2),'r-',VP(1,1,jPu1),VP(1,3,jPu1),'bo',VP(N1+1,1,jPu1),VP(N1+1,3,jPu1),'bo',VP(N+1,1,jPu1),VP(N+1,3,jPu1),'bo',VP(1,1,jPu2),VP(1,3,jPu2),'ro',VP(N1+1,1,jPu2),VP(N1+1,3,jPu2),'ro',VP(N+1,1,jPu2),VP(N+1,3,jPu2),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    camroll(-90)
    xlabel('x (m)')
    ylabel('w')
    %title('Modos longitudinais')
    grid
    legend('1^{o} modo','2^{o} modo');
    
    figure
    p = plot(VP(:,2,jPu1),VP(:,3,jPu1),'b-',VP(:,2,jPu2),VP(:,3,jPu2),'r-',VP(1,2,jPu1),VP(1,3,jPu1),'bo',VP(N1+1,2,jPu1),VP(N1+1,3,jPu1),'bo',VP(N+1,2,jPu1),VP(N+1,3,jPu1),'bo',VP(1,2,jPu2),VP(1,3,jPu2),'ro',VP(N1+1,2,jPu2),VP(N1+1,3,jPu2),'ro',VP(N+1,2,jPu2),VP(N+1,3,jPu2),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    axis square
    xlim([-4e-10 4e-10])
    ylim([-7.9e-10 7.9e-10])
    %camroll(-90)
    xlabel('v')
    ylabel('w')
    %title('Modos longitudinais')
    grid
    legend('1^{o} modo','2^{o} modo');
    
    %(c) Gráficos X(N) versus ubTil:
    figure
    p = plot(V0,uTilP(:,jPu1),'b-',V0,uTilP(:,jPu2),'r-',V0(1,1),uTilP(1,jPu1),'bo',V0(1,1),uTilP(1,jPu2),'ro',V0(N1+1,1),uTilP(N1+1,jPu1),'bo',V0(N1+1,1),uTilP(N1+1,jPu2),'ro',V0(N+1,1),uTilP(N+1,jPu1),'bo',V0(N+1,1),uTilP(N+1,jPu2),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    camroll(-90)
    xlabel('x (m)')
    ylabel('ub')
    %title('Modos longitudinais')
    grid
    legend('1^{o} modo','2^{o} modo');
    %
    %%
    %Observando os modos 1 e 4 (modos de torsão):
    %(a) Gráficos 3D:
    az = 90;
    el = 180;
    figure
    p = plot3(VP(:,2,jPtX1),VP(:,3,jPtX1),VP(:,1,jPtX1),VP(:,2,jPtX2),VP(:,3,jPtX2),VP(:,1,jPtX2));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos de torsão - Plano XZ')
    grid
    view(az, el);
    legend('1^{o} modo','2^{o} modo');

    az = 180;
    el = 180;
    figure
    p = plot3(VP(:,2,jPtX1),VP(:,3,jPtX1),VP(:,1,jPtX1),VP(:,2,jPtX2),VP(:,3,jPtX2),VP(:,1,jPtX2));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos de torsão - Plano XY')
    grid
    view(az, el);
    legend('1^{o} modo','2^{o} modo');
    
    az = 180;
    el = 90;
    figure
    p = plot3(VP(:,2,jPtX1),VP(:,3,jPtX1),VP(:,1,jPtX1),VP(:,2,jPtX2),VP(:,3,jPtX2),VP(:,1,jPtX2));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos de torsão - Plano XY')
    grid
    view(az, el);
    legend('1^{o} modo','2^{o} modo');
    
    %(b) Gráficos 2D:
    figure
    p = plot(VP(:,1,jPtX1),VP(:,2,jPtX1),'b-',VP(:,1,jPtX2),VP(:,2,jPtX2),'r-',VP(1,1,jPtX1),VP(1,2,jPtX1),'bo',VP(N1+1,1,jPtX1),VP(N1+1,2,jPtX1),'bo',VP(N+1,1,jPtX1),VP(N+1,2,jPtX1),'bo',VP(1,1,jPtX2),VP(1,2,jPtX2),'ro',VP(N1+1,1,jPtX2),VP(N1+1,2,jPtX2),'ro',VP(N+1,1,jPtX2),VP(N+1,2,jPtX2),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    camroll(-90)
    xlabel('x (m)')
    ylabel('v')
    %title('Modos de torsão - Plano XY')
    grid
    legend('1^{o} modo','2^{o} modo');
    
    figure
    p = plot(VP(:,1,jPtX1),VP(:,3,jPtX1),'b-',VP(:,1,jPtX2),VP(:,3,jPtX2),'r-',VP(1,1,jPtX1),VP(1,3,jPtX1),'bo',VP(N1+1,1,jPtX1),VP(N1+1,3,jPtX1),'bo',VP(N+1,1,jPtX1),VP(N+1,3,jPtX1),'bo',VP(1,1,jPtX2),VP(1,3,jPtX2),'ro',VP(N1+1,1,jPtX2),VP(N1+1,3,jPtX2),'ro',VP(N+1,1,jPtX2),VP(N+1,3,jPtX2),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    camroll(-90)
    xlabel('x (m)')
    ylabel('w')
    %title('Modos de torsão - Plano XZ')
    grid
    legend('1^{o} modo','2^{o} modo');
    
    figure
    p = plot(VP(:,2,jPtX1),VP(:,3,jPtX1),'b-',VP(:,2,jPtX2),VP(:,3,jPtX2),'r-',VP(1,2,jPtX1),VP(1,3,jPtX1),'bo',VP(N1+1,2,jPtX1),VP(N1+1,3,jPtX1),'bo',VP(N+1,2,jPtX1),VP(N+1,3,jPtX1),'bo',VP(1,2,jPtX2),VP(1,3,jPtX2),'ro',VP(N1+1,2,jPtX2),VP(N1+1,3,jPtX2),'ro',VP(N+1,2,jPtX2),VP(N+1,3,jPtX2),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    axis square
    %xlim([-4e-10 4e-10])
    %ylim([-7.9e-10 7.9e-10])
    %camroll(-90)
    xlabel('v')
    ylabel('w')
    %title('Modos de torsão - Plano XZ')
    grid
    legend('1^{o} modo','2^{o} modo');
    
    %(c) Gráficos X(N) versus tetaXTil:
    figure
    p = plot(V0,tetaXTilP(:,jPtX1),'b-',V0,tetaXTilP(:,jPtX2),'r-',V0(1,1),tetaXTilP(1,jPtX1),'bo',V0(N1+1,1),tetaXTilP(N1+1,jPtX1),'bo',V0(N+1,1),tetaXTilP(N+1,jPtX1),'bo',V0(1,1),tetaXTilP(1,jPtX2),'ro',V0(N1+1,1),tetaXTilP(N1+1,jPtX2),'ro',V0(N+1,1),tetaXTilP(N+1,jPtX2),'r-');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    camroll(-90)
    xlabel('x (m)')
    ylabel('teta_{x}')
    %title('Modos de torsão')
    grid
    legend('1^{o} modo','2^{o} modo');
    
    %%
    %Observando os modos 2 e 3 (modos reversos):
    %(a) Gráficos 3D:
    az = 90;
    el = 180;
    figure
    p = plot3(VP(:,2,jPv),VP(:,3,jPv),VP(:,1,jPv),VP(:,2,jPw),VP(:,3,jPw),VP(:,1,jPw));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos reversos - Plano XZ')
    grid
    view(az, el);
    legend('v','w');

    az = 180;
    el = 180;
    figure
    p = plot3(VP(:,2,jPv),VP(:,3,jPv),VP(:,1,jPv),VP(:,2,jPw),VP(:,3,jPw),VP(:,1,jPw));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos reversos - Plano XZ')
    grid
    view(az, el);
    legend('v','w');
    
    az = 180;
    el = 90;
    figure
    p = plot3(VP(:,2,jPv),VP(:,3,jPv),VP(:,1,jPv),VP(:,2,jPw),VP(:,3,jPw),VP(:,1,jPw));
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Modos reversos - Plano XZ')
    grid
    view(az, el);
    legend('v','w');
    
    %(b) Gráficos 2D:
    figure
    p = plot(VP(:,1,jPv),VP(:,2,jPv),'b-',VP(:,1,jPw),VP(:,2,jPw),'r-',VP(1,1,jPv),VP(1,2,jPv),'bo',VP(N1+1,1,jPv),VP(N1+1,2,jPv),'bo',VP(N+1,1,jPv),VP(N+1,2,jPv),'bo',VP(1,1,jPw),VP(1,2,jPw),'ro',VP(N1+1,1,jPw),VP(N1+1,2,jPw),'ro',VP(N+1,1,jPw),VP(N+1,2,jPw),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    camroll(-90)
    xlabel('x (m)')
    ylabel('v')
    %title('Modos reversos - Plano XY')
    grid
    legend('modo "v"','modo "w"');
    
    figure
    p = plot(VP(:,1,jPv),VP(:,3,jPv),'b-',VP(:,1,jPw),VP(:,3,jPw),'r-',VP(1,1,jPv),VP(1,3,jPv),'bo',VP(N1+1,1,jPv),VP(N1+1,3,jPv),'bo',VP(N+1,1,jPv),VP(N+1,3,jPv),'bo',VP(1,1,jPw),VP(1,3,jPw),'ro',VP(N1+1,1,jPw),VP(N1+1,3,jPw),'ro',VP(N+1,1,jPw),VP(N+1,3,jPw),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    camroll(-90)
    xlabel('x (m)')
    ylabel('w')
    %title('Modos reversos - Plano XZ')
    grid
    legend('modo "v"','modo "w"');
    
    figure
    p = plot(VP(:,2,jPv),VP(:,3,jPv),'b-',VP(:,2,jPw),VP(:,3,jPw),'r-',VP(1,2,jPv),VP(1,3,jPv),'bo',VP(N1+1,2,jPv),VP(N1+1,3,jPv),'bo',VP(N+1,2,jPv),VP(N+1,3,jPv),'bo',VP(1,2,jPw),VP(1,3,jPw),'ro',VP(N1+1,2,jPw),VP(N1+1,3,jPw),'ro',VP(N+1,2,jPw),VP(N+1,3,jPw),'ro');
    p(1).LineWidth = 1;
    p(2).LineWidth = 3;
    axis square
    %xlim([-4e-10 4e-10])
    %ylim([-7.9e-10 7.9e-10])
    %camroll(-90)
    xlabel('v')
    ylabel('w')
    %title('Modos reversos - Plano YZ')
    grid
    legend('modo "v"','modo "w"');
    
    %%
    %Observando os modos de 1 a 13:
    az = 90;
    el = 180;
    figure
    plot3(VP(:,2,7),VP(:,3,7),VP(:,1,7),VP(:,2,8),VP(:,3,8),VP(:,1,8),VP(:,2,9),VP(:,3,9),VP(:,1,9),VP(:,2,10),VP(:,3,10),VP(:,1,10),VP(:,2,11),VP(:,3,11),VP(:,1,11),VP(:,2,12),VP(:,3,12),VP(:,1,12),VP(:,2,14),VP(:,3,14),VP(:,1,14),VP(:,2,15),VP(:,3,15),VP(:,1,15));
    xlabel('y')
    ylabel('z')
    zlabel('x')
    %title('Modos 1 a 13 - Plano XZ')
    grid
    view(az, el);
    legend('7o modo','8o modo','9o modo','10o modo','11o modo','12o modo','14o modo','15o modo');

    az = 180;
    el = 180;
    figure
    plot3(VP(:,2,7),VP(:,3,7),VP(:,1,7),VP(:,2,8),VP(:,3,8),VP(:,1,8),VP(:,2,9),VP(:,3,9),VP(:,1,9),VP(:,2,10),VP(:,3,10),VP(:,1,10),VP(:,2,11),VP(:,3,11),VP(:,1,11),VP(:,2,12),VP(:,3,12),VP(:,1,12),VP(:,2,14),VP(:,3,14),VP(:,1,14),VP(:,2,15),VP(:,3,15),VP(:,1,15));
    xlabel('y')
    ylabel('z')
    zlabel('x')
    %title('Modos 1 a 13 - Plano XZ')
    grid
    view(az, el);
    legend('7o modo','8o modo','9o modo','10o modo','11o modo','12o modo','14o modo','15o modo');
    
    az = 180;
    el = 90;
    figure
    plot3(VP(:,2,7),VP(:,3,7),VP(:,1,7),VP(:,2,8),VP(:,3,8),VP(:,1,8),VP(:,2,9),VP(:,3,9),VP(:,1,9),VP(:,2,10),VP(:,3,10),VP(:,1,10),VP(:,2,11),VP(:,3,11),VP(:,1,11),VP(:,2,12),VP(:,3,12),VP(:,1,12),VP(:,2,14),VP(:,3,14),VP(:,1,14),VP(:,2,15),VP(:,3,15),VP(:,1,15));
    xlabel('y')
    ylabel('z')
    zlabel('x')
    %title('Modos 1 a 13 - Plano XZ')
    grid
    view(az, el);
    legend('7o modo','8o modo','9o modo','10o modo','11o modo','12o modo','14o modo','15o modo');
    
    %(b) Gráficos 2D:
    figure
    p = plot(VP(:,1,7),VP(:,2,7),'b-',VP(:,1,8),VP(:,2,8),'r-',VP(:,1,9),VP(:,2,9),'k-',VP(:,1,10),VP(:,2,10),'g-',VP(:,1,11),VP(:,2,11),'b-',VP(:,1,12),VP(:,2,12),'r-',VP(1,1,7),VP(1,2,7),'bo',VP(N1+1,1,7),VP(N1+1,2,7),'bo',VP(N+1,1,7),VP(N+1,2,7),'bo',VP(1,1,8),VP(1,2,8),'ro',VP(N1+1,1,8),VP(N1+1,2,8),'ro',VP(N+1,1,8),VP(N+1,2,8),'ro',VP(1,1,9),VP(1,2,9),'ko',VP(N1+1,1,9),VP(N1+1,2,9),'ko',VP(N+1,1,9),VP(N+1,2,9),'ko',VP(1,1,10),VP(1,2,10),'go',VP(N1+1,1,10),VP(N1+1,2,10),'go',VP(N+1,1,10),VP(N+1,2,10),'go',VP(1,1,11),VP(1,2,11),'ko',VP(N1+1,1,11),VP(N1+1,2,11),'ko',VP(N+1,1,11),VP(N+1,2,11),'ko',VP(1,1,12),VP(1,2,12),'ro',VP(N1+1,1,12),VP(N1+1,2,12),'ro',VP(N+1,1,12),VP(N+1,2,12),'ro');
    p(1).LineWidth = 1;
    p(1).LineStyle = ':';
    %p(1).Marker = '.';
    %p(1).MarkerSize = 6;
    p(2).LineWidth = 3;
    p(2).LineStyle = ':';
    p(4).LineStyle = '--';
    p(5).LineWidth = 3;
    p(5).LineStyle = '--';
    p(6).LineStyle = '-.';
    camroll(-90)
    xlabel('x (m)')
    ylabel('v')
    %title('Modos de torsão - Plano XY')
    grid
    legend('7^{o} modo','8^{o} modo','9^{o} modo','10^{o} modo','11^{o} modo','12^{o} modo');
    
    figure
    p = plot(VP(:,1,7),VP(:,3,7),'b-',VP(:,1,8),VP(:,3,8),'r-',VP(:,1,9),VP(:,3,9),'k-',VP(:,1,10),VP(:,3,10),'g-',VP(:,1,11),VP(:,3,11),'b-',VP(:,1,12),VP(:,3,12),'r-',VP(1,1,7),VP(1,3,7),'bo',VP(N1+1,1,7),VP(N1+1,3,7),'bo',VP(N+1,1,7),VP(N+1,3,7),'bo',VP(1,1,8),VP(1,3,8),'ro',VP(N1+1,1,8),VP(N1+1,3,8),'ro',VP(N+1,1,8),VP(N+1,3,8),'ro',VP(1,1,9),VP(1,3,9),'ko',VP(N1+1,1,9),VP(N1+1,3,9),'ko',VP(N+1,1,9),VP(N+1,3,9),'ko',VP(1,1,10),VP(1,3,10),'go',VP(N1+1,1,10),VP(N1+1,3,10),'go',VP(N+1,1,10),VP(N+1,3,10),'go',VP(1,1,11),VP(1,3,11),'ko',VP(N1+1,1,11),VP(N1+1,3,11),'ko',VP(N+1,1,11),VP(N+1,3,11),'ko',VP(1,1,12),VP(1,3,12),'ro',VP(N1+1,1,12),VP(N1+1,3,12),'ro',VP(N+1,1,12),VP(N+1,3,12),'ro');
    p(1).LineWidth = 1;
    p(1).LineStyle = ':';
    %p(1).Marker = '.';
    %p(1).MarkerSize = 6;
    p(2).LineWidth = 3;
    p(2).LineStyle = ':';
    p(4).LineStyle = '--';
    p(5).LineWidth = 3;
    p(5).LineStyle = '--';
    p(6).LineStyle = '-.';
    camroll(-90)
    xlabel('x (m)')
    ylabel('w')
    %title('Modos de torsão - Plano XZ')
    grid
    legend('7^{o} modo','8^{o} modo','9^{o} modo','10^{o} modo','11^{o} modo','12^{o} modo');
    
    figure
    p = plot(VP(:,2,7),VP(:,3,7),'b-',VP(:,2,8),VP(:,3,8),'r-',VP(:,2,9),VP(:,3,9),'k-',VP(:,2,10),VP(:,3,10),'g-',VP(:,2,11),VP(:,3,11),'b-',VP(:,2,12),VP(:,3,12),'r-',VP(1,2,7),VP(1,3,7),'bo',VP(N1+1,2,7),VP(N1+1,3,7),'bo',VP(N+1,2,7),VP(N+1,3,7),'bo',VP(1,2,8),VP(1,3,8),'ro',VP(N1+1,2,8),VP(N1+1,3,8),'ro',VP(N+1,2,8),VP(N+1,3,8),'ro',VP(1,2,9),VP(1,3,9),'ko',VP(N1+1,2,9),VP(N1+1,3,9),'ko',VP(N+1,2,9),VP(N+1,3,9),'ko',VP(1,2,10),VP(1,3,10),'go',VP(N1+1,2,10),VP(N1+1,3,10),'go',VP(N+1,2,10),VP(N+1,3,10),'go',VP(1,2,11),VP(1,3,11),'ko',VP(N1+1,2,11),VP(N1+1,3,11),'ko',VP(N+1,2,11),VP(N+1,3,11),'ko',VP(1,2,12),VP(1,3,12),'ro',VP(N1+1,2,12),VP(N1+1,3,12),'ro',VP(N+1,2,12),VP(N+1,3,12),'ro');
    p(1).LineWidth = 1;
    p(1).LineStyle = ':';
    %p(1).Marker = '.';
    %p(1).MarkerSize = 6;
    p(2).LineWidth = 3;
    p(2).LineStyle = ':';
    p(4).LineStyle = '--';
    p(5).LineWidth = 3;
    p(5).LineStyle = '--';
    p(6).LineStyle = '-.';
    camroll(-90)
    xlabel('v')
    ylabel('w')
    %title('Modos de torsão - Plano YZ')
    grid
    legend('7^{o} modo','8^{o} modo','9^{o} modo','10^{o} modo','11^{o} modo','12^{o} modo');
    
    %Observando os modos com maiores deslocamento dos rotores:
    epsLim = 0.1;
    limSup = epsLim*ones(1,6*(N+1)-rEstr);
    limInf = -epsLim*ones(1,6*(N+1)-rEstr);
    nRt = 1:(6*(N+1)-rEstr);
    rt1u = Abvt(1,:); rt1tX = Abvt(2,:);
    figure
    plot(nRt,rt1u,'o-',nRt,rt1tX,'x-',nRt,limSup,'k-.',nRt,limInf,'k-.')
    title('Variáveis u e tX no nó 1 em P')
    grid
    legend('rt1u','rt1tX')
    rtN1M1u = Abvt(6*(N1+1)-floor(rEstr/2)-5,:);
    rtN1M1v = Abvt(6*(N1+1)-floor(rEstr/2)-4,:);
    rtN1M1tZ = Abvt(6*(N1+1)-floor(rEstr/2)-3,:);
    rtN1M1w = Abvt(6*(N1+1)-floor(rEstr/2)-2,:);
    rtN1M1tY = Abvt(6*(N1+1)-floor(rEstr/2)-1,:);
    rtN1M1tX = Abvt(6*(N1+1)-floor(rEstr/2),:);
    figure
    plot(nRt,rtN1M1u,'o-',nRt,rtN1M1v,'x-',nRt,rtN1M1w,'*-',nRt,limSup,'k-.',nRt,limInf,'k-.')
    title('Variáveis u, v e w no nó N1+1 em P')
    grid
    legend('rtN1M1u','rtN1M1v','rtN1M1w')
    figure
    plot(nRt,rtN1M1tZ,'o-',nRt,rtN1M1tY,'x-',nRt,rtN1M1tX,'*-',nRt,limSup,'k-.',nRt,limInf,'k-.')
    title('Variáveis tZ, tY e tX no nó N1+1 em P')
    grid
    legend('rtN1M1tZ','rtN1M1tY','rtN1M1tX')
    rtNM1u = Abvt(6*(N+1)-rEstr-1,:); rtNM1tX = Abvt(6*(N+1)-rEstr,:);
    figure
    plot(nRt,rtNM1u,'o-',nRt,rtNM1tX,'x-',nRt,limSup,'k-.',nRt,limInf,'k-.')
    title('Variáveis u e tX no nó N+1 em P')
    grid
    legend('rtNM1u','rtNM1tX')
end
%}

%%
%Frequências naturais (Hz):
freqNatW2Tab = [sort(real(eig(M1restr\Krestr))),sort(real(eig(M1restr\KbRestr)))];
freqNatfTab = (2*pi)^-1*sqrt(freqNatW2Tab);
argumentos.freqNatfTab = freqNatfTab;

%Modelo reduzido a partir da formulação em EF:
opcaoPk = 2;
switch opcaoPk
    case 1
    case 2
        sigTilk = Abvl(1:jPu2,1);
        Pk = Abvt(:,1:jPu2);
    case 3
        sigTilk = [sigTiltX; sigTilu];
        Pk = [PtX, Pu];
    case 4
    case 5
    case 6
        prop = 1; %0 < prop < 1
        Ngr = Ng - rEstr;
        jPteste = floor(prop*Ngr);
        sigTilk = Abvl(1:jPteste,1);
        Pk = Abvt(:,1:jPteste);
    case 7
        sigTilk = Abvl(3:jPu2,1);
        Pk = Abvt(:,3:jPu2);
    otherwise
        disp('other value')
end
%{
if (1 < opcaoPk) && (opcaoPk < 3)
    %Obs.:
    %(1) Apenas para a opçãoPk = 2
    %(2) Preparação do metamodelo para lei de controle
    %(3) Frequências aproximadas para modos de corpo rígido:
    epsCR = 10^-4;
    sigTilk(jPtX0,1) = epsCR;
    sigTilk(jPu0,1) = epsCR;
end
%Obs.: Implementação suspensa
%(1) Pretende-se utilizar o modelo reduzido como representante da planta,
%para testes iniciais da lei de controle via LQR
%}
argumentos.sigTilk = sigTilk;
argumentos.Pk = Pk;
epsZero = 10^-9;
argumentos.epsZero = epsZero;
%{
%Teste:
M1rk = Pk.'*M1restr*Pk;
M1rk = M1rk - M1rk.*double(M1rk < epsZero);
ndiagM1rk = M1rk - diag(diag(M1rk));
NormNdiagM1rk = norm(ndiagM1rk);
Kbrk = Pk.'*KbRestr*Pk;
Kbrk = Kbrk - Kbrk.*double(Kbrk < 10^4*epsZero);
ndiagKbrk = Kbrk - diag(diag(Kbrk));
NormNdiagKbrk = norm(ndiagKbrk);
%ok!
%}

%%
%Matriz de amortecimento:
%(a) Amortecimento estrutural:
nSig = length(sigTiltX);
wRmin = freqNatKbw(jPtX1,1); %rad/s
wRmax = max(freqNatKbw); %rad/s
qsiMax = 1; %Valor máximo dos coeficientes de amortecimento
if nSig<3
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
    wNatTil = freqNatKbw;
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
Crestr = SigmaR.'*CgEF*SigmaR;
argumentos.Crestr = Crestr;
Crk = Pk.'*Crestr*Pk;
nCrk = length(Crk);
if adCconc < 1
    Crk = diag(Crk);
    Crk = diag(Crk);
    for i = 1:nCrk
        if abs(Crk(i,i))<epsZero
            Crk(i,i) = 0;
        end
    end
end
argumentos.Crk = Crk;

Id6NM1 = eye(6*(N+1));
Bef = [Id6NM1(:,6*1),Id6NM1(:,6*1-5)];
Br = SigmaR.'*Bef;
%{
%Equações de estado: Zp = Ab0*Z + Bb0*U
%Obs.:
%(1) Equação de origem: M1r*urpp + Cr*urp + Kr*ur = Br*U
Id6NM1r = eye(6*(N+1) - rEstr);
Ab0 = [zeros(6*(N+1) - rEstr), Id6NM1r; -M1restr\KbRestr, -M1restr\Crestr];
Bb0 = [zeros(6*(N+1) - rEstr, 2); Br];
%Controlabilidade:
CC = ctrb(Ab0,Bb0);
sCC = svd(CC);
idScc = double(sCC>0);
nScc = sum(idScc);
rCC = rank(CC);
%Obs.: 
%(1) rank(CC) = número de valores singulares de CC não nulos
%(2) Se rank(CC) = dim(Ab), então (Ab,Bb) é controlável
%(3) A matriz CC apresentou vários componentes "NaN"
%}

%%
%Equações dee estado (reduzidas): Zp = Abk*Z + Bbk*U
%Obs.:
%(1) Equação de origem: urkpp + Crk*urkp + Krk*urk = Brk*U
Brk = Pk.'*Br;
argumentos.Brk = Brk;
sPk = size(Pk);
Krk = diag(sigTilk.');
argumentos.Krk = Krk;
IdsPkCol = eye(sPk(1,2));
Ab = [zeros(sPk(1,2)), IdsPkCol; 
    -Krk, -Crk];
argumentos.Ab = Ab;
Bb = [zeros(sPk(1,2),2);
    Brk];
argumentos.Bb = Bb;
CC = ctrb(Ab,Bb);
sCC = svd(CC);
idScc = double(sCC>0);
nScc = sum(idScc);
rCC = rank(CC);
%Obs.: 
%(1) rank(CC) = número de valores singulares de CC não nulos
%(2) Se rank(CC) = dim(Ab), então (Ab,Bb) é controlável
%(3) A matriz CC apresentou vários componentes "NaN"

%%
%{
switch opcaoPk
    case 1
        Krk = diag(sigTilk.');
        IdsPkCol = eye(sPk(1,2));
        Ab = [zeros(sPk(1,2)), IdsPkCol; 
            -Krk, -Crk];
        Bb = [zeros(sPk(1,2),2);
            Brk];
    case 2
        Krk = diag(sigTilk.');
        IdsPkCol = eye(sPk(1,2));
        Ab = [zeros(sPk(1,2)), IdsPkCol; 
            -Krk, -Crk];
        Bb = [zeros(sPk(1,2),2);
            Brk];
    case 3
        Krk = diag(sigTilk.');
        IdsPkCol = eye(sPk(1,2));
        Ab = [zeros(sPk(1,2)), IdsPkCol; 
            -Krk, -Crk];
        Bb = [zeros(sPk(1,2),2);
            Brk];
    case 4
        Krk = diag(sigTilk.');
        IdsPkCol = eye(sPk(1,2));
        Ab = [zeros(sPk(1,2)), IdsPkCol; 
            -Krk, -Crk];
        Bb = [zeros(sPk(1,2),2);
            Brk];
    case 5
        Krk = diag(sigTilk.');
        IdsPkCol = eye(sPk(1,2));
        Ab = [zeros(sPk(1,2)), IdsPkCol; 
            -Krk, -Crk];
        Bb = [zeros(sPk(1,2),2);
            Brk];
    case 6
        Krk = diag(sigTilk.');
        IdsPkCol = eye(sPk(1,2));
        Ab = [zeros(sPk(1,2)), IdsPkCol; 
            -Krk, -Crk];
        Bb = [zeros(sPk(1,2),2);
            Brk];
    otherwise
        disp('other value')
end
%}
%{
switch opcaoPk
    case 1
        maxSigTilk = max(sigTilk);
        nSigRef = 0; %Não utilizado para opcaoPk=1
        djAjuste = sPk(1,2)-2;
        epsRk = sigTiltX(3,1)/maxSigTilk;
        A11 = Ab(1:sPk(1,2)+djAjuste,1:sPk(1,2)+djAjuste);
        A12 = Ab(1:sPk(1,2)+djAjuste,sPk(1,2)+djAjuste+1:2*sPk(1,2));
        A21 = Ab(sPk(1,2)+djAjuste+1:2*sPk(1,2),1:sPk(1,2)+djAjuste);
        A22 = Ab(sPk(1,2)+djAjuste+1:2*sPk(1,2),sPk(1,2)+djAjuste+1:2*sPk(1,2));
        B1 = Bb(1:sPk(1,2)+djAjuste,:);
        B2 = Bb(sPk(1,2)+djAjuste+1:2*sPk(1,2),:);
    case 2
        %{
        maxSigTilk = max(sigTilk);
        sigRef = sqrt(sigTilk(3,1)*maxSigTilk);
        idSigRef = double(sigTilk - sigRef < 0);
        nSigRef = sum(idSigRef);
        djAjuste = -2;
        epsRk = sigTilk(nSigRef+djAjuste,1)/maxSigTilk;
        A11 = Ab(1:sPk(1,2)+nSigRef+djAjuste,1:sPk(1,2)+nSigRef+djAjuste);
        A12 = Ab(1:sPk(1,2)+nSigRef+djAjuste,sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        A21 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),1:sPk(1,2)+nSigRef+djAjuste);
        A22 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        B1 = Bb(1:sPk(1,2)+nSigRef+djAjuste,:);
        B2 = Bb(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),:);
        %}
        maxSigTilk = max(sigTilk);
        epsRk = sigTilk(3,1)/maxSigTilk;
        A11 = Ab(1:sPk(1,2)+2,1:sPk(1,2)+2);
        A12 = Ab(1:sPk(1,2)+2,sPk(1,2)+2+1:2*sPk(1,2));
        A21 = Ab(sPk(1,2)+2+1:2*sPk(1,2),1:sPk(1,2)+2);
        A22 = Ab(sPk(1,2)+2+1:2*sPk(1,2),sPk(1,2)+2+1:2*sPk(1,2));
        B1 = Bb(1:sPk(1,2)+2,:);
        B2 = Bb(sPk(1,2)+2+1:2*sPk(1,2),:);
        ad = 1;
    case 3
        maxSigTilk = max(sigTilk);
        nSigRef = 0; %Não utilizado para opcaoPk=3
        djAjuste = sPk(1,2)-2;
        epsRk = sigTiltX(3,1)/maxSigTilk;
        A11 = Ab(1:sPk(1,2)+djAjuste,1:sPk(1,2)+djAjuste);
        A12 = Ab(1:sPk(1,2)+djAjuste,sPk(1,2)+djAjuste+1:2*sPk(1,2));
        A21 = Ab(sPk(1,2)+djAjuste+1:2*sPk(1,2),1:sPk(1,2)+djAjuste);
        A22 = Ab(sPk(1,2)+djAjuste+1:2*sPk(1,2),sPk(1,2)+djAjuste+1:2*sPk(1,2));
        B1 = Bb(1:sPk(1,2)+djAjuste,:);
        B2 = Bb(sPk(1,2)+djAjuste+1:2*sPk(1,2),:);
    case 4
        maxSigTilk = max(sigTilk);
        sigRef = sqrt(sigTilk(3,1)*maxSigTilk);
        idSigRef = double(sigTilk - sigRef < 0);
        nSigRef = sum(idSigRef);
        djAjuste = 0;
        epsRk = sigTilk(nSigRef+djAjuste,1)/maxSigTilk;
        A11 = Ab(1:sPk(1,2)+nSigRef+djAjuste,1:sPk(1,2)+nSigRef+djAjuste);
        A12 = Ab(1:sPk(1,2)+nSigRef+djAjuste,sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        A21 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),1:sPk(1,2)+nSigRef+djAjuste);
        A22 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        B1 = Bb(1:sPk(1,2)+nSigRef+djAjuste,:);
        B2 = Bb(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),:);
    case 5
        maxSigTilk = max(sigTilk);
        sigRef = sqrt(sigTilk(3,1)*maxSigTilk);
        idSigRef = double(sigTilk - sigRef < 0);
        nSigRef = sum(idSigRef);
        djAjuste = 0;
        epsRk = sigTilk(nSigRef+djAjuste,1)/maxSigTilk;
        A11 = Ab(1:sPk(1,2)+nSigRef+djAjuste,1:sPk(1,2)+nSigRef+djAjuste);
        A12 = Ab(1:sPk(1,2)+nSigRef+djAjuste,sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        A21 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),1:sPk(1,2)+nSigRef+djAjuste);
        A22 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        B1 = Bb(1:sPk(1,2)+nSigRef+djAjuste,:);
        B2 = Bb(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),:);
    case 6
        maxSigTilk = max(sigTilk);
        sigRef = sqrt(sigTilk(3,1)*maxSigTilk);
        idSigRef = double(sigTilk - sigRef < 0);
        nSigRef = sum(idSigRef);
        djAjuste = 0;
        epsRk = sigTilk(nSigRef+djAjuste,1)/maxSigTilk;
        A11 = Ab(1:sPk(1,2)+nSigRef+djAjuste,1:sPk(1,2)+nSigRef+djAjuste);
        A12 = Ab(1:sPk(1,2)+nSigRef+djAjuste,sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        A21 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),1:sPk(1,2)+nSigRef+djAjuste);
        A22 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        B1 = Bb(1:sPk(1,2)+nSigRef+djAjuste,:);
        B2 = Bb(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),:);
    case 7
        minSigTilk = min(sigTilk);
        maxSigTilk = max(sigTilk);
        sigRef = sqrt(minSigTilk*maxSigTilk);
        idSigRef = double(sigTilk - sigRef < 0);
        nSigRef = sum(idSigRef);
        djAjuste = -2;
        epsRk = sigTilk(nSigRef+djAjuste,1)/maxSigTilk;
        A11 = Ab(1:sPk(1,2)+nSigRef+djAjuste,1:sPk(1,2)+nSigRef+djAjuste);
        A12 = Ab(1:sPk(1,2)+nSigRef+djAjuste,sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        A21 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),1:sPk(1,2)+nSigRef+djAjuste);
        A22 = Ab(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2));
        B1 = Bb(1:sPk(1,2)+nSigRef+djAjuste,:);
        B2 = Bb(sPk(1,2)+nSigRef+djAjuste+1:2*sPk(1,2),:);
    otherwise
        disp('other value')
end
%{
A0 = A11 - A12*(A22\A21);
B0 = B1 - A12*(A22\B2);
fXps = @(Xs,Us)(A0*Xs + B0*Us);
repositorio.fXps = fXps;
fZs = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
repositorio.fZs = fZs;
fZpf = @(Zf,Uf)(epsRk^-1*(A22*Zf + B2*Uf));
repositorio.fZpf = fZpf;
%}

%%
%Saídas para o modelo sem separação das dinâmicas:
%{
nTetaXNM1 = 6*(N+1)-rEstr;
nuNM1 = 6*(N+1)-rEstr-1;
nvN1M1 = 6*(N1+1)-floor(rEstr/2)-4;
nwN1M1 = 6*(N1+1)-floor(rEstr/2)-2;
%}
%{
nTetaXNM1 = 6*(N+1);
nuNM1 = 6*(N+1)-5;
nvN1M1 = 6*(N1+1)-4;
nwN1M1 = 6*(N1+1)-2;
Pk = SigmaR*Pk;
sPk = size(Pk);
%}
%y1 = tetaXpNM1 = C1*Z:
C1 = [zeros(1,sPk(1,2)), Pk(nTetaXNM1,:)];
%y2 = uNM1 = C2*Z:
C2 = [Pk(nuNM1,:), zeros(1,sPk(1,2))];
%y3 = vN1M1 = C3*Z:
C3 = [Pk(nvN1M1,:), zeros(1,sPk(1,2))];
%y4 = wN1M1 = C4*Z:
C4 = [Pk(nwN1M1,:), zeros(1,sPk(1,2))];
%y5 = tetaXNM1 = C5*Z:
C5 = [Pk(nTetaXNM1,:), zeros(1,sPk(1,2))];
Pk = SigmaR.'*Pk;

%Saídas para o modelo com separação das dinâmicas:
Id2sPkCol = eye(2*sPk(1,2));
Axz = Id2sPkCol(:,1:sPk(1,2)+2);
%
C1Axz = C1*Axz;
C2Axz = C2*Axz;
C5Axz = C5*Axz;
%{
Azz = Id2sPkCol(:,sPk(1,2)+2+1:2*sPk(1,2));
C1Azz = C1*Azz; obsC1Azz = double(C1Azz - 10^-4 > 0); nC1Azz = sum(obsC1Azz);
C2Azz = C2*Azz;
C3Axz = C3*Axz;
C3Azz = C3*Azz;
C4Axz = C4*Axz;
C4Azz = C4*Azz;
C5Azz = C5*Azz;
%ok!
%}
beta0s = C1Axz;
alfa0s = kNM1*C2Axz;
gama0s = C5Axz;

%%
%Estado e entradas permanentes da dinâmica lenta (X0, U0):
wRef = argumentos.wRef;
adNr = 0;
NrRef = adNr*M1*g;
Alqr = [[A0, B0].'*[A0, B0], [beta0s, zeros(1,2)].', [alfa0s, zeros(1,2)].', [gama0s, zeros(1,2)].'; 
        [beta0s, zeros(1,2)],                     0,                     0,                     0;
        [alfa0s, zeros(1,2)],                     0,                     0,                     0;
        [gama0s, zeros(1,2)],                     0,                     0,                     0];
sA0 = size(A0);
sB0 = size(B0);
bLQR = [zeros(sA0(1,1)+2,1); 
                       wRef;
                      NrRef;
                          0];
Xu0lmb = Alqr\bLQR;
X0 = Xu0lmb(1:sA0(1,1),1);
nX0 = length(X0);
for i=1:nX0
    if abs(X0(i))<epsZero
        X0(i) = 0;
    end
end
argumentos.X0 = X0;
U0 = Xu0lmb((sA0(1,1) + 1):(sA0(1,1) + sB0(1,2)),1);
nU0 = length(U0);
for i=1:nU0
    if abs(U0(i))<epsZero
        U0(i) = 0;
    end
end
argumentos.U0 = U0;
%{
lmb0 = Xu0lmb((sA0(1,1) + sB0(1,2) + 1):(sA0(1,1) + sB0(1,2) + 3),1);
%}

%}

%Separação das dinâmicas lenta e rápida:
[ ABvEps, infoS ] = fModOrdMMod( argumentos );
%(1) Ganho do regulador:
Klqr = ABvEps{1,46};
argumentos.Klqr = Klqr;
%(2) Matrizes da dinâmica, com absorvedor de vibração:
Abs = ABvEps{1,13};
argumentos.Abs = Abs;
Bbs = ABvEps{1,14};
argumentos.Bbs = Bbs;
%(3) Matrizes da dinâmica externa:
beta1s = ABvEps{1,16}; 
beta0Bbs = ABvEps{1,17};
alfa2s = ABvEps{1,20}; 
alfa1Bbs = ABvEps{1,21};
%(4) Matrizes dinâmica interna e difeomorfismo:
AmiS = ABvEps{1,22};
argumentos.AmiS = AmiS;
sigDifs = ABvEps{1,23};
argumentos.sigDifs = sigDifs;
%Matriz permuta de linhas Zb = Perl*Z e Z = Perl.'*Zb:
Perl = infoS.Perlqv;
argumentos.Perl = Perl;

%%
%{
%Ponto de equilíbrio (modelos da planta, por parâmetros concentrados):
%Obs.: considerando u1_Eq = 0
Ku = [k1 + k2, -k2; -k2, k2];
Wu = [M2b + dm; M1]*g;
uTilEqMods = Ku\Wu;
u2_Eq = uTilEqMods(1,1);
u3_Eq = uTilEqMods(2,1);
argumentos.uTilEqMods = uTilEqMods;

Id6 = eye(6);
fteta2 = @(uTilAc)(Id6(1,:)*uTilAc);
repositorio.fteta2 = fteta2;
fyL = @(uTilAc)(Id6(2,:)*uTilAc);
repositorio.fyL = fyL;
fzL = @(uTilAc)(Id6(3,:)*uTilAc);
repositorio.fzL = fzL;
fu2 = @(uTilAc)(Id6(4,:)*uTilAc);
repositorio.fu2 = fu2;
fu3 = @(uTilAc)(Id6(5,:)*uTilAc);
repositorio.fu3 = fu3;
fteta3 = @(uTilAc)(Id6(6,:)*uTilAc);
repositorio.fteta3 = fteta3;
uTilAc_Eq = [0; 0; 0; u2_Eq; u3_Eq; 0];
argumentos.uTilAc_Eq = uTilAc_Eq;
%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RUBBING
opcao_sign = argumentos.opcao_sign;
Kip = argumentos.Kip;
Dr2 = argumentos.Dr2;
R2 = Dr2/2;
fg = argumentos.fg;
delta0 = argumentos.delta0;
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
repositorio.fT2cRubc = fT2c;
%(b) Modelo de atrito descontínuo na origem:
fFydc = @(yL,zL,teta2p,yLp,zLp)(-fFn(fuL(yL,zL))*(yL/fuL(yL,zL)) +...
    fFtdc(fuL(yL,zL),fvCfi(yL,zL,teta2p,yLp,zLp))*(zL/fuL(yL,zL)));
repositorio.fFyRubdc = fFydc;
fFzdc = @(yL,zL,teta2p,yLp,zLp)(-fFn(fuL(yL,zL))*(zL/fuL(yL,zL)) -...
    fFtdc(fuL(yL,zL),fvCfi(yL,zL,teta2p,yLp,zLp))*(yL/fuL(yL,zL)));
repositorio.fFzRubdc = fFzdc;
fT2dc = @(yL,zL,teta2p,yLp,zLp)(-miIpdc(fvCfi(yL,zL,teta2p,yLp,zLp))*...
    fFn(fuL(yL,zL))*(R2 + fg));
repositorio.fT2dcRubc = fT2dc;

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
            inc_dsignS*dsignSdvCfi(vCfi)*...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vetor fEF(u,up,upp):
%(1) Modelo contínuo de atrito:
fcEF = @(u,up,upp)(fGfcEF( u, up, upp, argumentos, repositorio ));
repositorio.fcEF = fcEF;
%(2) Modelo descontínuo de atrito:
fdcEF = @(u,up,upp)(fGfdcEF( u, up, upp, argumentos, repositorio ));
repositorio.fdcEF = fdcEF;

%Força peso sobressalente:
[ W ] = fGw( argumentos, repositorio );
[ Fu ] = fGfu( uTilEq, argumentos, repositorio );
Wsb = W - KgEF*uTilEq - Fu;
Wsbr = SigmaR.'*Wsb;

%Definição de variáveis e funções, em torno do ponto de equilíbrio:
%(1) Modelo contínuo de atrito:
fcEFb = @(ub,ubp,ubpp)(KgEF*uTilEq + fcEF(ub + uTilEq,ubp,ubpp) + Wsb);
repositorio.fcEFb = fcEFb;
%(2) Modelo descontínuo de atrito:
fdcEFb = @(ub,ubp,ubpp)(KgEF*uTilEq + fdcEF(ub + uTilEq,ubp,ubpp) + Wsb);
repositorio.fdcEFb = fdcEFb;

%Derivadas:
dudubT = IdNg;
dupdubpT = IdNg;
duppdubppT = IdNg;
%(1) Modelo contínuo de atrito:
dfcEFduT = @(u,up,upp)(fGdfcEFdu( u, up, upp, argumentos, repositorio ));
dfcEFbdubT = @(ub,ubp,ubpp)(dfcEFduT(ub + uTilEq,ubp,ubpp)*dudubT);
repositorio.dfcEFbdubT = dfcEFbdubT;
dfcEFdupT = @(u,up)(fGdfcEFdup( u, up, argumentos, repositorio ));
dfcEFbdubpT = @(ub,ubp)(dfcEFdupT(ub + uTilEq,ubp)*dupdubpT);
repositorio.dfcEFbdubpT = dfcEFbdubpT;
%(2) Modelo descontínuo de atrito:
dfdcEFduT = @(u,up,upp)(fGdfdcEFdu( u, up, upp, argumentos, repositorio ));
dfdcEFbdubT = @(ub,ubp,ubpp)(dfdcEFduT(ub + uTilEq,ubp,ubpp)*dudubT);
repositorio.dfdcEFbdubT = dfdcEFbdubT;
dfdcEFdupT = @(u,up)(fGdfdcEFdup( u, up, argumentos, repositorio ));
dfdcEFbdubpT = @(ub,ubp)(dfdcEFdupT(ub + uTilEq,ubp)*dupdubpT);
repositorio.dfdcEFbdubpT = dfdcEFbdubpT;
%(3) Derivada em relação à aceleração:
dfEFduppT = @(u)(fGdfEFdupp( u, argumentos, repositorio ));
dfEFbdubppT = @(ub)(dfEFduppT(ub + uTilEq)*duppdubppT);
repositorio.dfEFbdubppT = dfEFbdubppT;

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
fcEFbAprox = @(ub,ubp)(KgEF*uTilEq + fcEFaprox(ub + uTilEq,ubp));
repositorio.fcEFbAprox = fcEFbAprox;
fcEFbAtr = @(ub,ubp)(fGfcEFbAtr( ub, ubp, argumentos, repositorio ));
repositorio.fcEFbAtr = fcEFbAtr;
%(2) Modelo descontínuo de atrito:
fdcEFbAprox = @(ub,ubp)(KgEF*uTilEq + fdcEFaprox(ub + uTilEq,ubp));
repositorio.fdcEFbAprox = fdcEFbAprox;
fdcEFbAtr = @(ub,ubp)(fGfdcEFbAtr( ub, ubp, argumentos, repositorio ));
repositorio.fdcEFbAtr = fdcEFbAtr;

%Derivadas:
%(1) Modelo contínuo de atrito:
dfcEFduT = @(u,up,upp)(fGdfcEFdu( u, up, upp, argumentos, repositorio ));
dfcEFbdubT = @(ub,ubp,ubpp)(dfcEFduT(ub + uTilEq,ubp,ubpp)*dudubT);
repositorio.dfcEFbdubT = dfcEFbdubT;
dfcEFdupT = @(u,up)(fGdfcEFdup( u, up, argumentos, repositorio ));
dfcEFbdubpT = @(ub,ubp)(dfcEFdupT(ub + uTilEq,ubp)*dupdubpT);
repositorio.dfcEFbdubpT = dfcEFbdubpT;
%(1a) Aproximação:
dfcEFduTaprox = @(u,up)(fGdfcEFduAprox( u, up, argumentos, repositorio ));
dfcEFbdubTaprox = @(ub,ubp)(dfcEFduTaprox(ub + uTilEq,ubp)*dudubT);
repositorio.dfcEFbdubTaprox = dfcEFbdubTaprox;
dfcEFdupTaprox = @(u,up)(fGdfcEFdupAprox( u, up, argumentos, repositorio ));
dfcEFbdubpTaprox = @(ub,ubp)(dfcEFdupTaprox(ub + uTilEq,ubp)*dupdubpT);
repositorio.dfcEFbdubpTaprox = dfcEFbdubpTaprox;
dfcEFdubTAtr = @(ub,ubp)(fGdfcEFdubAtr( ub, ubp, argumentos, repositorio ));
repositorio.dfcEFdubTAtr = dfcEFdubTAtr;
dfcEFdubpTAtr = @(ub,ubp)(fGdfcEFdubpAtr( ub, ubp, argumentos, repositorio ));
repositorio.dfcEFdubpTAtr = dfcEFdubpTAtr;

%(2) Modelo descontínuo de atrito:
dfdcEFduT = @(u,up,upp)(fGdfdcEFdu( u, up, upp, argumentos, repositorio ));
dfdcEFbdubT = @(ub,ubp,ubpp)(dfdcEFduT(ub + uTilEq,ubp,ubpp)*dudubT);
repositorio.dfdcEFbdubT = dfdcEFbdubT;
dfdcEFdupT = @(u,up)(fGdfdcEFdup( u, up, argumentos, repositorio ));
dfdcEFbdubpT = @(ub,ubp)(dfdcEFdupT(ub + uTilEq,ubp)*dupdubpT);
repositorio.dfdcEFbdubpT = dfdcEFbdubpT;
%(2a) Aproximação:
dfdcEFduTaprox = @(u,up)(fGdfdcEFduAprox( u, up, argumentos, repositorio ));
dfdcEFbdubTaprox = @(ub,ubp)(dfdcEFduTaprox(ub + uTilEq,ubp)*dudubT);
repositorio.dfdcEFbdubTaprox = dfdcEFbdubTaprox;
dfdcEFdupTaprox = @(u,up)(fGdfdcEFdupAprox( u, up, argumentos, repositorio ));
dfdcEFbdubpTaprox = @(ub,ubp)(dfdcEFdupTaprox(ub + uTilEq,ubp)*dupdubpT);
repositorio.dfdcEFbdubpTaprox = dfdcEFbdubpTaprox;
dfdcEFdubTAtr = @(ub,ubp)(fGdfdcEFdubAtr( ub, ubp, argumentos, repositorio ));
repositorio.dfdcEFdubTAtr = dfdcEFdubTAtr;
dfdcEFdubpTAtr = @(ub,ubp)(fGdfdcEFdubpAtr( ub, ubp, argumentos, repositorio ));
repositorio.dfdcEFdubpTAtr = dfdcEFdubpTAtr;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEI DE CONTROLE

%{
%Frequência de corte (para definição da superfície de escorregamento):
nLmb = 0.25;
wCtr2 = w2d;
wCtr4 = freqNatKbw(jPtX2+1,1); %Menor frequência não modelada da planta
lambda_1 = (1-nLmb)*wCtr2 + nLmb*wCtr4;
argumentos.lambda_ctr_y1 = lambda_1;
%}

%Força normal:
Nref = argumentos.Nref;
wRef = argumentos.wRef;
gama_TD = 0.1; %(N/(rad*s^-1)) Parâmetro de ajuste para cada ponto de operação
Psi = @(beta)(beta*exp(wRef/beta) - wRef - Nref/gama_TD - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
beta0 = 10^6; %rad/s
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
while fiBeta(beta0) <= beta0/2
    beta0 = beta0/2;
end
[ beta_TD ] = f_Newton_Raphson( beta0, Psi, dPsidbeta, argumentos );

Nr0b = @(w)(gama_TD*beta_TD*exp(w/beta_TD) - gama_TD*(w + beta_TD));
repositorio.Nr0b = Nr0b;
dNr0bdw = @(w)(gama_TD*(exp(w/beta_TD) - 1));
repositorio.dNr0bdw = dNr0bdw;
d2Nr0bdw2 = @(w)((gama_TD*exp(w/beta_TD))/beta_TD);
repositorio.d2Nr0bdw2 = d2Nr0bdw2;
d3Nr0bdw3 = @(w)((gama_TD*exp(w/beta_TD))/beta_TD^2);
repositorio.d3Nr0bdw3 = d3Nr0bdw3;

%{
tetaP = @(w3)(betaS_hat*(1/k1_hat + 1/k2_hat + 1/k3_hat)*Normalb0(w3));
repositorio.tetaP = tetaP;
dNormalb0dw3 = @(teta3p)((gama_TD*exp(teta3p/beta_TD) - gama_TD)*double(teta3p>=0));
repositorio.dNormalb0dw3 = dNormalb0dw3;
dtetaPdw3 = @(w3)(betaS_hat*(1/k1_hat + 1/k2_hat + 1/k3_hat)*dNormalb0dw3(w3));
repositorio.dtetaPdw3 = dtetaPdw3;
tetapP = @(w3,w3p)(dtetaPdw3(w3)*w3p);
repositorio.tetapP = tetapP;
d2Normalb0dw32 = @(teta3p)((1/beta_TD)*(gama_TD*exp(teta3p/beta_TD))*double(teta3p>=0));
repositorio.d2Normalb0dw32 = d2Normalb0dw32;
d3Normalb0dw33 = @(teta3p)((1/beta_TD^2)*(gama_TD*exp(teta3p/beta_TD))*double(teta3p>=0));
repositorio.d3Normalb0dw33 = d3Normalb0dw33;
d4Normalb0dw34 = @(teta3p)((1/beta_TD^3)*(gama_TD*exp(teta3p/beta_TD))*double(teta3p>=0));
repositorio.d4Normalb0dw34 = d4Normalb0dw34;
d2tetaPdw32 = @(w3)(betaS_hat*(1/k1_hat + 1/k2_hat + 1/k3_hat)*d2Normalb0dw32(w3));
dtetapPdw3 = @(w3,w3p)(d2tetaPdw32(w3)*w3p);
repositorio.dtetapPdw3 = dtetapPdw3;
dtetapPdw3p = @(w3)(dtetaPdw3(w3));
repositorio.dtetapPdw3p = dtetapPdw3p;
%}

%Formas vetorizadas:
Normalb0vt = @(w)((gama_TD*beta_TD*exp(w/beta_TD) -...
    gama_TD*(w + beta_TD*ones(size(w)))));
repositorio.Normalb0vt = Normalb0vt;

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

%(b) Modelo descontínuo:
mi30dc_hat = argumentos.mi0_hat;
mi31dc_hat = argumentos.mi1_hat;
mi32dc_hat = argumentos.mi2_hat;
beta1mi3dc_hat = argumentos.beta1mi_hat;
beta2mi3dc_hat = argumentos.beta2mi_hat;
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

%Derivadas (para integração):
dfiAtrHatdnu = @(nu)(-pi*(tanh(pi*nu)^2 - 1)  +...
    (2*lambdami_hat*nmi)/(2*nmi + nu^(2*nmi) - 1) -...
    (4*lambdami_hat*nmi^2*nu^(2*nmi))/(2*nmi + nu^(2*nmi) - 1)^2);
repositorio.dfiAtrdnu = dfiAtrHatdnu;
dfiAtrHatdnuVt = @(nu)(-pi*(tanh(pi*nu).^2 - 1)  +...
    (2*lambdami_hat*nmi)./(2*nmi + nu.^(2*nmi) - 1) -...
    (4*lambdami_hat*nmi^2*nu.^(2*nmi))./(2*nmi + nu.^(2*nmi) - 1).^2);
repositorio.dfiAtrdnuVt = dfiAtrHatdnuVt;

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
%}

%{
%Torque de atrito:
Id2 = eye(2);
fteta3 = @(upTil)(Id2(2,:)*upTil);
repositorio.fteta3meta = fteta3;
%(a) Modelo contínuo:
Tatr3c0Hat = @(teta3p)(mi3cHat(teta3p)*Normalb0(teta3p)*R3_hat);
Tc0Hat = @(uTilbp)([0; Tatr3c0Hat(fteta3(uTilbp))]);
repositorio.Tc0Hat = Tc0Hat;
%(b) Modelo descontínuo:
Tatr3dc0Hat = @(teta3p)(mi3dcHat(teta3p)*Normalb0(teta3p)*R3_hat);
Tdc0Hat = @(uTilbp)([0; Tatr3dc0Hat(fteta3(uTilbp))]);
repositorio.Tdc0Hat = Tdc0Hat;

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
repositorio.d2fiAtrdnu2 = d2fiAtrHatdnu2;
d2fiAtrHatdnu2Vt = @(nu)(2*pi^2*tanh(pi*nu).*(tanh(pi*nu).^2 - 1) -...
    (8*lambdami_hat*nmi^2*nu.^(2*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^2 +...
    (16*lambdami_hat*nmi^3*nu.^(4*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^3 -...
    (4*lambdami_hat*nmi^2*nu.^(2*nmi - 1)*(2*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^2);
repositorio.d2fiAtrdnu2Vt = d2fiAtrHatdnu2Vt;

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
%Obs.: Expressão obtida a partir d4fiAtrHatdnu4, mas não verificada simb.
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

%{
%%
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
%Funções úteis:
fcZp = @(Z,Uc)(AbHat*Z - FcHat(Z) + Bhat*Uc);
repositorio.fcZp = fcZp;
fdcZp = @(Z,Uc)(AbHat*Z - FdcHat(Z) + Bhat*Uc);
repositorio.fdcZp = fdcZp;

%%
%FORMA NORMAL:
C1 = Id4(4,:);

%(1) Termos auxiliares:
beta0 = C1;
argumentos.beta0 = beta0;
beta1Hat = beta0*AbHat;
argumentos.beta1Hat = beta1Hat;
beta2Hat = beta1Hat*AbHat;
argumentos.beta2Hat = beta2Hat;
beta3Hat = beta2Hat*AbHat;
argumentos.beta3Hat = beta3Hat;

Qsi0cHat = @(Z)(fQ0cHat( Z, argumentos, repositorio ).');
repositorio.Qsi0cHat = Qsi0cHat;
Qsi1cHat = @(Z)(fQ1cHat( Z, argumentos, repositorio ).');
repositorio.Qsi1cHat = Qsi1cHat;

Qsi0dcHat = @(Z)(fQsi0dcHat( Z, argumentos, repositorio ));
repositorio.Qsi0dcHat = Qsi0dcHat;
Qsi1dcHat = @(Z)(fQsi1dcHat( Z, argumentos, repositorio ));
repositorio.Qsi1dcHat = Qsi1dcHat;

%(3) Dinâmica interna: mi = Ami*Z
epsilon = 10^-5;
q = epsilon;
%q = 1;
argumentos.qAmi = q;
%q1 = 0;
q1 = 1;
q3 = (Ib20_hat/K1_hat)*epsilon;
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

%(2) Dinâmica externa (forma normal): yppp = a[d]c(Z) + b*Uc
acHat = @(Z)(beta3Hat*Z - beta2Hat*FcHat(Z) - Qsi1cHat(Z)*(AbHat*Z - FcHat(Z)));
repositorio.acHat = acHat;
adcHat = @(Z)(beta3Hat*Z - beta2Hat*FdcHat(Z) - Qsi1dcHat(Z).'*(AbHat*Z - FdcHat(Z)));
repositorio.adcHat = adcHat;
bHat = beta2Hat*Bhat;
argumentos.bHat = bHat;
%}

%(1) Trajetória desejada: 
%(1.a) y1d = tetaXp1d:
t0_9wRef = argumentos.t0_9wRef;
alfa = argumentos.alfa_trajetoria_desejada;
fi = @(nu)(tanh(alfa*nu) - (alfa*nu)/((alfa^2*nu^2)/3 + 1));
y1d = @(t)(wRef*fi(t/t0_9wRef));
repositorio.y1d = y1d;
dfidnu = @(nu)((6*alfa^3*nu^2)/(alfa^2*nu^2 + 3)^2 -...
    (3*alfa)/(alfa^2*nu^2 + 3) -...
    alfa*(tanh(alfa*nu)^2 - 1));
y1pd = @(t)((wRef/t0_9wRef)*dfidnu(t/t0_9wRef));
repositorio.y1pd = y1pd;
d2fidnu2 = @(nu)(2*alfa^2*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1) +...
    (18*alfa^3*nu)/(alfa^2*nu^2 + 3)^2 -...
    (24*alfa^5*nu^3)/(alfa^2*nu^2 + 3)^3);
y1ppd = @(t)((wRef/t0_9wRef^2)*d2fidnu2(t/t0_9wRef));
repositorio.y1ppd = y1ppd;
d3fidnu3 = @(nu)((18*alfa^3)/(alfa^2*nu^2 + 3)^2 -...
    2*alfa^3*(tanh(alfa*nu)^2 - 1)^2 -...
    4*alfa^3*tanh(alfa*nu)^2*(tanh(alfa*nu)^2 - 1) -...
    (144*alfa^5*nu^2)/(alfa^2*nu^2 + 3)^3 +...
    (144*alfa^7*nu^4)/(alfa^2*nu^2 + 3)^4);
y1pppd = @(t)((wRef/t0_9wRef^3)*d3fidnu3(t/t0_9wRef));
repositorio.y1pppd = y1pppd;
d4fidnu4 = @(nu)(16*alfa^4*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1)^2 -...
    (360*alfa^5*nu)/(alfa^2*nu^2 + 3)^3 +...
    8*alfa^4*tanh(alfa*nu)^3*(tanh(alfa*nu)^2 - 1) +...
    (1440*alfa^7*nu^3)/(alfa^2*nu^2 + 3)^4 -...
    (1152*alfa^9*nu^5)/(alfa^2*nu^2 + 3)^5);
y1ppppd = @(t)((wRef/t0_9wRef^4)*d4fidnu4(t/t0_9wRef));
repositorio.y1ppppd = y1ppppd;
d5fidnu5 = @(nu)((6480*alfa^7*nu^2)/(alfa^2*nu^2 + 3)^4 -...
    16*alfa^5*(tanh(alfa*nu)^2 - 1)^3 -...
    16*alfa^5*tanh(alfa*nu)^4*(tanh(alfa*nu)^2 - 1) -...
    (360*alfa^5)/(alfa^2*nu^2 + 3)^3 -...
    (17280*alfa^9*nu^4)/(alfa^2*nu^2 + 3)^5 +...
    (11520*alfa^11*nu^6)/(alfa^2*nu^2 + 3)^6 -...
    88*alfa^5*tanh(alfa*nu)^2*(tanh(alfa*nu)^2 - 1)^2);
y1pppppd = @(t)((wRef/t0_9wRef^5)*d5fidnu5(t/t0_9wRef));
repositorio.y1pppppd = y1pppppd;
d6fidnu6 = @(nu)((15120*alfa^7*nu)/(alfa^2*nu^2 + 3)^4 +...
    272*alfa^6*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1)^3 +...
    32*alfa^6*tanh(alfa*nu)^5*(tanh(alfa*nu)^2 - 1) -...
    (120960*alfa^9*nu^3)/(alfa^2*nu^2 + 3)^5 +...
    (241920*alfa^11*nu^5)/(alfa^2*nu^2 + 3)^6 -...
    (138240*alfa^13*nu^7)/(alfa^2*nu^2 + 3)^7 +...
    416*alfa^6*tanh(alfa*nu)^3*(tanh(alfa*nu)^2 - 1)^2);
y1ppppppd = @(t)((wRef/t0_9wRef^6)*d6fidnu6(t/t0_9wRef));
repositorio.y1ppppppd = y1ppppppd;

%Vetorização:
fiVt = @(nu)(tanh(alfa*nu) - (alfa*nu)./((alfa^2*nu.^2)/3 + ones(size(nu))));
y1dVt = @(T)(wRef*fiVt(T/t0_9wRef));
repositorio.y1dVt = y1dVt;
dfidnuVt = @(nu)((6*alfa^3*nu.^2)./(alfa^2*nu.^2 + 3*ones(size(nu))).^2 -...
    (3*alfa*ones(size(nu)))./(alfa^2*nu.^2 + 3*ones(size(nu))) -...
    alfa*(tanh(alfa*nu).^2 - ones(size(nu))));
y1pdVt = @(T)((wRef/t0_9wRef)*dfidnuVt(T/t0_9wRef));
repositorio.y1pdVt = y1pdVt;
d2fidnu2Vt = @(nu)(2*alfa^2*tanh(alfa*nu).*(tanh(alfa*nu).^2 - ones(size(nu))) +...
    (18*alfa^3*nu)./(alfa^2*nu.^2 + 3*ones(size(nu))).^2 -...
    (24*alfa^5*nu.^3)./(alfa^2*nu.^2 + 3*ones(size(nu))).^3);
y1ppdVt = @(T)((wRef/t0_9wRef^2)*d2fidnu2Vt(T/t0_9wRef));
repositorio.y1ppdVt = y1ppdVt;
d3fidnu3Vt = @(nu)((18*alfa^3*ones(size(nu)))./(alfa^2*nu.^2 + 3*ones(size(nu))).^2 -...
    2*alfa^3*(tanh(alfa*nu).^2 - ones(size(nu))).^2 -...
    4*alfa^3*tanh(alfa*nu).^2.*(tanh(alfa*nu).^2 - ones(size(nu))) -...
    (144*alfa^5*nu.^2)./(alfa^2*nu.^2 + 3*ones(size(nu))).^3 +...
    (144*alfa^7*nu.^4)./(alfa^2*nu.^2 + 3*ones(size(nu))).^4);
y1pppdVt = @(T)((wRef/t0_9wRef^3)*d3fidnu3Vt(T/t0_9wRef));
repositorio.y1pppdVt = y1pppdVt;
d4fidnu4Vt = @(nu)(16*alfa^4*tanh(alfa*nu).*(tanh(alfa*nu).^2 - 1).^2 -...
    (360*alfa^5*nu)./(alfa^2*nu.^2 + 3).^3 +...
    8*alfa^4*tanh(alfa*nu).^3.*(tanh(alfa*nu).^2 - 1) +...
    (1440*alfa^7*nu.^3)./(alfa^2*nu.^2 + 3).^4 -...
    (1152*alfa^9*nu.^5)./(alfa^2*nu.^2 + 3).^5);
y1ppppdVt = @(t)((wRef/t0_9wRef^4)*d4fidnu4Vt(t/t0_9wRef));
repositorio.y1ppppdVt = y1ppppdVt;
%(1.b) y2d = ub1d:
kHatNM1 = kNM1;
y2d = @(t)(Nr0b(y1d(t))/kHatNM1);
repositorio.y2d = y2d;
y2pd = @(t)((1/kHatNM1)*dNr0bdw(y1d(t))*y1pd(t));
repositorio.y2pd = y2pd;
y2ppd = @(t)((1/kHatNM1)*(d2Nr0bdw2(y1d(t))*(y1pd(t))^2 + dNr0bdw(y1d(t))*y1ppd(t)));
repositorio.y2ppd = y2ppd;
y2pppd = @(t)((1/kHatNM1)*((d3Nr0bdw3(y1d(t))*(y1pd(t))^3 + d2Nr0bdw2(y1d(t))*2*(y1pd(t))*y1ppd(t)) +...
    (d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppd(t) + dNr0bdw(y1d(t))*y1pppd(t))));
repositorio.y2pppd = y2pppd;

ydTils = @(t)([y1d(t); y2d(t); y2pd(t)]);
repositorio.ydTils = ydTils;
ypdTils = @(t)([y1pd(t); y2pd(t); y2ppd(t)]);
repositorio.ypdTils = ypdTils;

%{
%(2) Lei de controle:
daChatdqsi = @(Z)(fdaChatdqsi( Z, argumentos, repositorio ));
daC = @(Z)(dqsiTil.'*daChatdqsi(Z));
daDcHatdqsi = @(Z)(fdaDcHatdqsi( Z, argumentos, repositorio ));
daDc = @(Z)(dqsiTil.'*daDcHatdqsi(Z));
%dbHatdqsi = fdbHatdqsi( argumentos );
%dbHat = dqsiTil.'*dbHatdqsi;

%sigmaDmax = bHat/abs(dbHat);
%argumentos.sigmaDmax = sigmaDmax;

sigmaF = 1;
FcLC = @(Z)(sigmaF*abs(daC(Z)));
repositorio.FcLC = FcLC;
FdcLC = @(Z)(sigmaF*abs(daDc(Z)));
repositorio.FdcLC = FdcLC;
%sigmaD = 1;
%Dlc = sigmaD*(1/bHat)*abs(dbHat);
%argumentos.Dlc = Dlc;

%Obs.: X = [y; yp; ypp; mi]
y = @(X)(Id4(1,:)*X);
yp = @(X)(Id4(2,:)*X);
ypp = @(X)(Id4(3,:)*X);

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
%Kfi = lambda_1/(1 - Dlc);
%argumentos.Kfi = Kfi;
%Vetorização:
yVt = @(XX)(XX*Id4(:,1));
ydVt = @(T)(wRef*fiVt(T/t0_9wRef));
ypVt = @(XX)(XX*Id4(:,2));
ypdVt = @(T)((wRef/t0_9wRef)*dfidnuVt(T/t0_9wRef));
yppVt = @(XX)(XX*Id4(:,3));
y2rVt = @(T,XX)(yppdVt(T) - 2*lambda_1*(ypVt(XX) - ypdVt(T)) -...
    lambda_1^2*(yVt(XX) - ydVt(T)));
sVt = @(T,XX)(yppVt(XX) - y2rVt(T,XX));
repositorio.sVt = sVt;

%Funções adaptadas para a trajetoria desejada:
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

%Camada limite dinâmica:
DfiBd = @(fipBd)(double(fipBd>=0)*1/(1 + Dlc) + double(fipBd<=0)*1/(1 - Dlc));
repositorio.Dfi = DfiBd;
GamaFiBd = @(fipBd)(DfiBd(fipBd)*fipBd);
repositorio.GamaFiBd = GamaFiBd;

%Expressões de ganhos e lei de controle:
Kcb = @(t,Z,X,fiBd,mid,Zd0)(Kc(t,Z,X) - Kcd(t,mid,Zd0) + Kfi*fiBd);
repositorio.Kcb = Kcb;
Kdcb = @(t,Z,X,fiBd,mid,Zd0)(Kdc(t,Z,X) - Kdcd(t,mid,Zd0) + Kfi*fiBd);
sat = @(nu)(tanh(pi*nu));
repositorio.sat = sat;
Ucc = @(t,Z,X,fiBd,mid,Zd0)(bHat^-1*(UcTilc(t,Z,X) - Kcb(t,Z,X,fiBd,mid,Zd0)*sat(s(t,X)/fiBd)));
repositorio.Ucc = Ucc;
Ucdc = @(t,Z,X,fiBd,mid,Zd0)(bHat^-1*(UcTildc(t,Z,X) - Kdcb(t,Z,X,fiBdTil,mid,Zd0)*sat(s(t,X)/fiBd)));
repositorio.Ucdc = Ucdc;

%Expressões para integração:
Uccu = @(t,uTilc,uTilcp,fiBd,mid,Zd0)(Ucc(t,fZ(uTilc,uTilcp),fcXhat(uTilc,uTilcp),fiBd,mid,Zd0));
repositorio.Uccu = Uccu;
Ucdcu = @(t,uTilc,uTilcp,fiBd,mid,Zd0)(Ucdc(t,fZ(uTilc,uTilcp),fdcXhat(uTilc,uTilcp),fiBd,mid,Zd0));
repositorio.Ucdcu = Ucdcu;

%Expressões para a dinâmica interna:
y3rMid = @(t)(ypppd(t));
UcTilcMid = @(t,mid,Zd0)(y3rMid(t) - acHat(fcZdHat(t,mid,Zd0)));
UccMid = @(t,mid,Zd0)(bHat^-1*(UcTilcMid(t,mid,Zd0)));
repositorio.UccMid = UccMid;
UcTildcMid = @(t,mid,Zd0)(y3rMid(t) - adcHat(fdcZdHat(t,mid,Zd0)));
UcdcMid = @(t,mid,Zd0)(bHat^-1*(UcTildcMid(t,mid,Zd0)));
repositorio.UcdcMid = UcdcMid;

%Derivadas (para integração):
dGamadfipBd = @(fipBd)(DfiBd(fipBd));
repositorio.dGamadfipBd = dGamadfipBd;

dydXT = Id4(1,:);
dypdXT = Id4(2,:);
dyppdXT = Id4(3,:);
dy3rdXT = -2*lambda_1*dyppdXT -...
    lambda_1^2*dypdXT;
dUcTildXT = dy3rdXT;
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
dy2rdXT = -2*lambda_1*dypdXT - lambda_1^2*dydXT;
dsdXT = dyppdXT - dy2rdXT;
dsatsFidXT = @(t,X,fiBd)(dsatdnu(s(t,X)/fiBd)*(1/fiBd)*dsdXT);
dUccdXT = @(t,Z,X,fiBd,mid,Z0)(bHat^-1*...
    (dUcTildXT -...
    (dKcbdXT(t,Z,X,fiBd,mid,Z0)*sat(s(t,X)/fiBd) + Kcb(t,Z,X,fiBd,mid,Z0)*dsatsFidXT(t,X,fiBd))));
DUccDZT = @(t,Z,X,fiBd,mid,Z0)(dUccdZT(t,Z,X,fiBd) + dUccdXT(t,Z,X,fiBd,mid,Z0)*dPsiDifcHatdZT(Z));
repositorio.DUccDZT = DUccDZT;

%dKbdFiBd = Kfi;
dsatsFidFiBd = @(t,X,fiBd)(dsatdnu(s(t,X)/fiBd)*(-s(t,X)/fiBd^2));
dUccdFiBd = @(t,Z,X,fiBd,mid,Z0)(-bHat^-1*...
    (dKbdFiBd*sat(s(t,X)/fiBd) + Kcb(t,Z,X,fiBd,mid,Z0)*dsatsFidFiBd(t,X,fiBd)));
repositorio.dUccdFiBd = dUccdFiBd;
dUcdcdFiBd = @(t,Z,X,fiBd,mid,Z0)(-bHat^-1*(dKbdFiBd*sat(s(t,X)/fiBd) +...
    Kdcb(t,Z,X,fiBd,mid,Z0)*dsatsFidFiBd(t,X,fiBd)));
repositorio.dUcdcdFiBd = dUcdcdFiBd;
dUcdcudFiBd = @(t,uTilc,uTilcp,fiBd,mid,Z0)(dUcdcdFiBd(t,fZ(uTilc,uTilcp),fdcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.dUcdcudFiBd = dUcdcudFiBd;
dXdmi = Ad;
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
dUccudmid = @(t,uTilc,uTilcp,fiBd,mid,Z0)(dUccdmid(t,fcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.dUccudmid = dUccudmid;
dUcdcudmid = @(t,uTilc,uTilcp,fiBd,mid,Z0)(dUcdcdmid(t,fdcXhat(uTilc,uTilcp),fiBd,mid,Z0));
repositorio.dUcdcudmid = dUcdcudmid;

%Expressões para a dinâmica interna:
dUcTilcMiddmid = @(t,mid,Z0)(-daChatdZT(fcZdHat(t,mid,Z0))*dInvPsiDifcdXTHat(fXd(t,mid),Z0)*dXdmi);
dUccMiddmid = @(t,mid,Z0)(bHat^-1*(dUcTilcMiddmid(t,mid,Z0)));
repositorio.dUccMiddmid = dUccMiddmid;
dUcTildcMiddmid = @(t,mid,Z0)(-daDcHatdZT(fdcZdHat(t,mid,Z0))*dInvPsiDifdcdXTHat(fXd(t,mid),Z0)*dXdmi);
dUcdcMiddmid = @(t,mid,Z0)(bHat^-1*(dUcTildcMiddmid(t,mid,Z0)));
repositorio.dUcdcMiddmid = dUcdcMiddmid;

%Outras expressões úteis:
fcdZdHatdmid = @(t,mid,Zd0)(dInvPsiDifcdXTHat(fXd(t,mid),Zd0)*dXdmi);
fdcdZdHatdmid = @(t,mid,Zd0)(dInvPsiDifdcdXTHat(fXd(t,mid),Zd0)*dXdmi);
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

%Derivada parcial de Ucc e Ucdc em relação a t:
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

%%
%Lei de controle por modos delizantes para o sistema Zbp = Abs*Zb + Bbs*Uc:
%FORMA NORMAL, dinâmica externa:
a1s = @(Zb)(beta1s*Zb);
b1s = beta0Bbs(1,1);
a2s = @(Zb)(alfa2s*Zb);
b2s = alfa1Bbs(1,2);

aTils = @(Zb)([a1s(Zb); a2s(Zb)]);
repositorio.aTils = aTils;
Es = diag([b1s, b2s]);
argumentos.Es = Es;

sigmaA1s = 0.1;
da1s = @(Zb)(sigmaA1s*a1s(Zb));
sigmaA2s = 0.1;
da2s = @(Zb)(sigmaA2s*a2s(Zb));
sigmaB1s = 0.1;
db1s = sigmaB1s*b1s;
sigmaB2s = 0.1;
db2s = sigmaB2s*b2s;

datils = @(Zb)([da1s(Zb); da2s(Zb)]);
FlcTils = @(Zb)(abs(datils(Zb)));
repositorio.FlcTils = FlcTils;
D11s = abs(db1s)/b1s;
D22s = abs(db2s)/b2s;
DlcS = diag([D11s, D22s]);
argumentos.DlcS = DlcS;

Id2jPu2 = eye(2*jPu2);
y1s = @(Xb)(Id2jPu2(1,:)*Xb);
y2s = @(Xb)(Id2jPu2(2,:)*Xb);
y2pS = @(Xb)(Id2jPu2(3,:)*Xb);

%Frequência de corte:
sigjPtX2 = Abvl(jPtX2);
sigjPtX2M1 = Abvl(jPtX2+1);
alfaLmbTx = 0.25;
lambda1 = (1 - alfaLmbTx)*sigjPtX2 + alfaLmbTx*sigjPtX2M1;
sigjPu2 = Abvl(jPu2);
sigjPu2M1 = Abvl(jPu2+1);
alfaLmbU = 0.25;
lambda2 = (1 - alfaLmbU)*sigjPu2 + alfaLmbU*sigjPu2M1;

y1r0s = @(t)(y1d(t));
y2r1s = @(t,Xb)(y2pd(t) - lambda2*(y2s(Xb) - y2d(t)));
yr1Tils = @(t,Xb)([y1r0s(t); y2r1s(t,Xb)]);
repositorio.yr1Tils = yr1Tils;

%Variáveis de escorregamento:
s1s = @(t,Xb)(y1s(Xb) - y1r0s(t));
s2s = @(t,Xb)(y2pS(Xb) - y2r1s(t));
sTils = @(t,Xb)([s1s(t,Xb); s2s(t,Xb)]);
repositorio.sTils = sTils;

y1r1s = @(t)(y1pd(t));
y2r2s = @(t,Xb)(y2ppd(t) - lambda2*(y2pS(Xb) - y2pd(t)));
yr3Tils = @(t,Xb)([y1r1s(t); y2r2s(t,Xb)]);
y2r2dS = @(t)(y2ppd(t));
y3rdTils = @(t)([y1r1s(t); y2r2dS(t)]);
y1pr1s = @(t)(y1ppd(t));
y2pr2dS = @(t)(y2pppd(t));
yr3dTilpS = @(t)([y1pr1s(t); y2pr2dS(t)]);

%Vetor de estado desejado:
nSigDif0s = length(ydTils(0));
argumentos.nSigDif0s = nSigDif0s;
nMids = 2*jPu2 - nSigDif0s;
argumentos.nMids = nMids;
Ads = Id2jPu2(:,nSigDif0s+1:2*jPu2);
Bds = Id2jPu2(:,1:nSigDif0s);
fXds = @(t,mid)(Ads*mid + Bds*ydTils(t)); %Obs.: Xd = [ydTil(t); mid]
fZbds = @(t,mid)(sigDifs\fXds(t,mid));
repositorio.fZbds = fZbds;

UcHats = @(t,Zb,Xb)(yr3Tils(t,Xb) - aTils(Zb));
UcdHats = @(t,mid)(y3rdTils(t) - aTils(fZbds(t,mid)));

%Ganhos:
Id2 = eye(2);
etaLc1s = 1;
etaLc2s = 1;
etaLcTils = [etaLc1s; etaLc2s];
%(1) Expressões diretas:
KlcS = @(t,Zb,Xb)((Id2 - DlcS)\(FlcTils(Zb) + etaLcTils + DlcS*abs(UcHats(t,Zb,Xb))));
repositorio.KlcS = KlcS;
Klcds = @(t,mid)((Id2 - DlcS)\(FlcTils(fZbds(t,mid)) + etaLcTils + DlcS*abs(UcdHats(t,mid))));
repositorio.Klcds = Klcds;
KfiS = diag([lambda1/(1 - D11s), lambda2/(1 - D22s)]);
argumentos.KfiS = KfiS;
KbLcS = @(t,Zb,Xb,mid,fiBdTil)(KlcS(t,Zb,Xb) - Klcds(t,mid) + KfiS*fiBdTil);
repositorio.KbLcS = KbLcS;
%(2) Derivadas:
da1sdZbT = beta1s;
da2sdZbT = alfa2s;
dda1sdZbT = sigmaA1s*da1sdZbT;
dda2sdZbT = sigmaA2s*da2sdZbT;
ddaTilsdZbT = [dda1sdZbT; dda2sdZbT];
dFlcTilsdZbT = @(Zb)(diag(sign(datils(Zb).'))*ddaTilsdZbT);
dZbdXbT = sigDifs\Id2jPu2;
dXbddt = @(t)(Bds*ypdTils(t));
dZbddt = @(t)(dZbdXbT*dXbddt(t));
dFlcTilsdt = @(t,mid)(dFlcTilsdZbT(fZbds(t,mid))*dZbddt(t));
daTilsdZbT = [da1sdZbT; da2sdZbT];
dUcdHatsdt = @(t)(yr3dTilpS(t) - daTilsdZbT*dZbddt(t));
dabsUcdHatsdt = @(t,mid)(diag(sign(UcdHats(t,mid).'))*dUcdHatsdt(t));
dKlcdsdt = @(t,mid)((Id2 - DlcS)\(dFlcTilsdt(t,mid) + DlcS*dabsUcdHatsdt(t,mid)));
dXbdmidT = Ads;
dZbddmidT = dZbdXbT*dXbdmidT;
dFlcTilsdmidT = @(t,mid)(dFlcTilsdZbT(fZbds(t,mid))*dZbddmidT);
dUcdHatsdmidT = -daTilsdZbT*dZbddmidT;
%Obs.:
%UcdHat = @(t,mid)(yr3dTil(t) - aTils(fZbd(t,mid)));
dabsUcdHatsdmidT = @(t,mid)(diag(sign(UcdHats(t,mid).'))*dUcdHatsdmidT);
dKlcdsdmidT = @(t,mid)((Id2 - DlcS)\(dFlcTilsdmidT(t,mid) + DlcS*dabsUcdHatsdmidT(t,mid)));
repositorio.dKlcddmidT = dKlcdsdmidT;
DKlcdsDt = @(t,mid,mipd)(dKlcdsdt(t,mid) + dKlcdsdmidT(t,mid)*mipd);
repositorio.DKlcdsDt = DKlcdsDt;

%Dinâmica interna: 
%Obs.: mip = GamaMis*mi + SigmaMis*yTil
GamaMis = AmiS*Abs*(sigDifs\Ads);
argumentos.GamaMis = GamaMis;
SigmaMis = AmiS*Abs*(sigDifs\Bds);
argumentos.SigmaMis = SigmaMis;
fmipds = @(t,mid)(GamaMis*mid + SigmaMis*ydTils(t));
repositorio.fmipds = fmipds;
fmippds = @(t,mipd)(GamaMis*mipd + SigmaMis*ypdTils(t));
repositorio.fmippds = fmippds;

%Camada limite dinâmica:
%Obs.: DfiBd(fipBdTil)*fipBdTil + Kfi*fiBdTil - Klcd(t,mid) = 0
%(1) Expressões diretas:
D1fiS = diag([1/(1+D11s), 1/(1+D22s)]); 
D2fiS = diag([1/(1+D11s), 1/(1-D22s)]);
D3fiS = diag([1/(1-D11s), 1/(1+D22s)]);
D4fiS = diag([1/(1-D11s), 1/(1-D22s)]);
DfiBdS = @(fipBdTil)(...
    D1fiS*double((Id2(1,:)*fipBdTil >= 0) && (Id2(2,:)*fipBdTil >= 0)) +...
    D2fiS*double((Id2(1,:)*fipBdTil >= 0) && (Id2(2,:)*fipBdTil < 0)) +...
    D3fiS*double((Id2(1,:)*fipBdTil < 0) && (Id2(2,:)*fipBdTil >= 0)) +...
    D4fiS*double((Id2(1,:)*fipBdTil < 0) && (Id2(2,:)*fipBdTil < 0)));
repositorio.DfiBdS = DfiBdS;
GamaFiBdS = @(fipBdTil)(DfiBdS(fipBdTil)*fipBdTil);
repositorio.GamaFiBdS = GamaFiBdS;
%(2) Derivadas:
dGamaFiSdfipBdT = @(fipBdTil)(DfiBdS(fipBdTil));
repositorio.dGamaFiSdfipBdT = dGamaFiSdfipBdT;
ffippBds = @(t,mid,mipd,fipBd)(DfiBdS(fipBd)\(-KfiS*fipBd + DKlcdsDt(t,mid,mipd)));
repositorio.ffippBds = ffippBds;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lei de controle por modos delizantes para o sistema Zp = A*Z + Bb*Uc:

%{
ABvEps(1,33) = {beta0NM1}; ABvEps(1,34) = {beta1NM1}; 
ABvEps(1,35) = {beta2NM1}; ABvEps(1,36) = {beta3NM1}; 
ABvEps(1,37) = {beta4NM1}; ABvEps(1,28) = {beta5NM1};
ABvEps(1,39) = {beta4NM1Bb};
ABvEps(1,40) = {alfa0NM1}; ABvEps(1,41) = {alfa1NM1}; 
ABvEps(1,42) = {alfa2NM1};
ABvEps(1,43) = {alfa1NM1Bb};
%(11) Matrizes dinâmica interna e difeomorfismo (com LQR):
ABvEps(1,44) = {AmiNM1}; ABvEps(1,45) = {sigDifNM1};
%}

beta5NM1 = ABvEps{1,28};
beta4NM1Bb = ABvEps{1,39};
alfa2NM1 = ABvEps{1,42};
alfa1NM1Bb = ABvEps{1,43};
sigDifNM1 = ABvEps{1,45};
AmiNM1 = ABvEps{1,44};
argumentos.AmiNM1 = AmiNM1;

%FORMA NORMAL, dinâmica externa:
a1fn = @(Z)(beta5NM1*Z);
b1fn = beta4NM1Bb(1,1);
a2fn = @(Z)(alfa2NM1*Z);
b2fn = alfa1NM1Bb(1,2);

aTil = @(Z)([a1fn(Z); a2fn(Z)]);
repositorio.aTil = aTil;
Efn = diag([b1fn, b2fn]);
argumentos.Efn = Efn;

sigmaA1 = 0.1;
da1fn = @(Z)(sigmaA1*a1fn(Z));
sigmaA2 = 0.1;
da2fn = @(Z)(sigmaA2*a2fn(Z));
sigmaB1 = 0.1;
db1fn = sigmaB1*b1fn;
sigmaB2 = 0.1;
db2fn = sigmaB2*b2fn;

datil = @(Z)([da1fn(Z); da2fn(Z)]);
FlcTilFN = @(Z)(abs(datil(Z)));
repositorio.FlcTilFN = FlcTilFN;
D11fn = abs(db1fn)/b1fn;
D22fn = abs(db2fn)/b2fn;
DlcFN = diag([D11fn, D22fn]);
argumentos.DlcFN = DlcFN;

y1 = @(X)(Id2jPu2(1,:)*X);
y1p = @(X)(Id2jPu2(2,:)*X);
y1pp = @(X)(Id2jPu2(3,:)*X);
y1ppp = @(X)(Id2jPu2(4,:)*X);
y1pppp = @(X)(Id2jPu2(5,:)*X);
y2 = @(X)(Id2jPu2(6,:)*X);
y2p = @(X)(Id2jPu2(7,:)*X);

ydTil = @(t)([y1d(t); y1pd(t); y1ppd(t); y1pppd(t); y1ppppd(t); y2d(t); y2pd(t)]);
repositorio.ydTil = ydTil;
ypdTil = @(t)([y1pd(t); y1ppd(t); y1pppd(t); y1ppppd(t); y1pppppd(t); y2pd(t); y2ppd(t)]);
repositorio.ypdTil = ypdTil;

%Frequência de corte:
sigjPtX2 = Abvl(jPtX2);
sigjPtX2M1 = Abvl(jPtX2+1);
alfaLmbTx = 0.25;
lambda1 = (1 - alfaLmbTx)*sigjPtX2 + alfaLmbTx*sigjPtX2M1;
sigjPu2 = Abvl(jPu2);
sigjPu2M1 = Abvl(jPu2+1);
alfaLmbU = 0.25;
lambda2 = (1 - alfaLmbU)*sigjPu2 + alfaLmbU*sigjPu2M1;

y4r1 = @(t,X)(y1ppppd(t) - 4*lambda1*(y1ppp(X) - y1pppd(t)) -...
    6*lambda1^2*(y1pp(X) - y1ppd(t)) - 4*lambda1^3*(y1p(X) - y1pd(t)) -...
    lambda1^4*(y1(X) - y1d(t)));
y1r2 = @(t,X)(y2pd(t) - lambda2*(y2(X) - y2d(t)));
y5rTil = @(t,X)([y4r1(t,X); y1r2(t,X)]);
repositorio.y5rTil = y5rTil;

%Variáveis de escorregamento:
s1 = @(t,X)(y1pppp(X) - y4r1(t,X));
s2 = @(t,X)(y2p(X) - y1r2(t,X));
sTil = @(t,X)([s1(t,X); s2(t,X)]);
repositorio.sTil = sTil;

y5r1 = @(t,X)(y1pppppd(t) - 4*lambda1*(y1pppp(X) - y1ppppd(t)) -...
    6*lambda1^2*(y1ppp(X) - y1pppd(t)) - 4*lambda1^3*(y1pp(X) - y1ppd(t)) -...
    lambda1^4*(y1p(X) - y1pd(t)));
y2r2 = @(t,X)(y2ppd(t) - lambda2*(y2p(X) - y2pd(t)));
y7rTil = @(t,X)([y5r1(t,X); y2r2(t,X)]);
y5r1d = @(t)(y1pppppd(t));
y2r2d = @(t)(y2ppd(t));
y7rdTil = @(t)([y5r1d(t); y2r2d(t)]);

y5pr1d = @(t)(y1ppppppd(t));
y2pr2d = @(t)(y2pppd(t));
y7rdTilp = @(t)([y5pr1d(t); y2pr2d(t)]);

%Vetor de estado desejado:
nSigDif0 = length(ydTil(0));
argumentos.nSigDif0 = nSigDif0;
nMid = 2*jPu2 - nSigDif0;
argumentos.nMid = nMid;
Ad = Id2jPu2(:,nSigDif0+1:2*jPu2);
Bd = Id2jPu2(:,1:nSigDif0);
fXd = @(t,mid)(Ad*mid + Bd*ydTil(t)); %Obs.: Xd = [ydTil(t); mid]
fZd = @(t,mid)(sigDifNM1\fXd(t,mid));
repositorio.fZd = fZd;

UcHat = @(t,Z,X)(y7rTil(t,X) - aTil(Z));
UcdHat = @(t,mid)(y7rdTil(t) - aTil(fZd(t,mid)));

%Ganhos:
Id2 = eye(2);
etaLc1 = 1;
etaLc2 = 1;
etaLcTil = [etaLc1; etaLc2];
%(1) Expressões diretas:
Klc = @(t,Z,X)((Id2 - DlcFN)\(FlcTilFN(Z) + etaLcTil + DlcFN*abs(UcHat(t,Z,X))));
repositorio.Klc = Klc;
Klcd = @(t,mid)((Id2 - DlcFN)\(FlcTilFN(fZd(t,mid)) + etaLcTil + DlcFN*abs(UcdHat(t,mid))));
repositorio.Klcd = Klcd;
Kfi = diag([lambda1/(1 - D11fn), lambda2/(1 - D22fn)]);
argumentos.Kfi = Kfi;
KbLc = @(t,Z,X,mid,fiBdTil)(Klc(t,Z,X) - Klcd(t,mid) + Kfi*fiBdTil);
repositorio.KbLc = KbLc;
%(2) Derivadas:
da1fndZT = beta5NM1;
da2fndZT = alfa2NM1;
dda1fndZT = sigmaA1*da1fndZT;
dda2fndZT = sigmaA2*da2fndZT;
ddaTildZT = [dda1fndZT; dda2fndZT];
dFlcTildZT = @(Z)(diag(sign(datil(Z).'))*ddaTildZT);
dZdXT = sigDifNM1\Id2jPu2;
dXddt = @(t)(Bd*ypdTil(t));
dZddt = @(t)(dZdXT*dXddt(t));
dFlcTildt = @(t,mid)(dFlcTildZT(fZd(t,mid))*dZddt(t));
daTildZT = [da1fndZT; da2fndZT];
dUcdHatdt = @(t)(y7rdTilp(t) - daTildZT*dZddt(t));
dabsUcdHatdt = @(t,mid)(diag(sign(UcdHat(t,mid).'))*dUcdHatdt(t));
dKlcddt = @(t,mid)((Id2 - DlcFN)\(dFlcTildt(t,mid) + DlcFN*dabsUcdHatdt(t,mid)));
dXdmidT = Ad;
dZddmidT = dZdXT*dXdmidT;
dFlcTildmidT = @(t,mid)(dFlcTildZT(fZd(t,mid))*dZddmidT);
dUcdHatdmidT = -daTildZT*dZddmidT;
%Obs.:
%UcdHat = @(t,mid)(yr3dTil(t) - aTils(fZbd(t,mid)));
dabsUcdHatdmidT = @(t,mid)(diag(sign(UcdHat(t,mid).'))*dUcdHatdmidT);
dKlcddmidT = @(t,mid)((Id2 - DlcFN)\(dFlcTildmidT(t,mid) + DlcFN*dabsUcdHatdmidT(t,mid)));
repositorio.dKlcddmidT = dKlcddmidT;
DKlcdDt = @(t,mid,mipd)(dKlcddt(t,mid) + dKlcddmidT(t,mid)*mipd);
repositorio.DKlcdDt = DKlcdDt;

%Dinâmica interna: 
%Obs.: mip = GamaMis*mi + SigmaMis*yTil
GamaMi = AmiNM1*Ab*(sigDifNM1\Ad);
argumentos.GamaMi = GamaMi;
SigmaMi = AmiNM1*Ab*(sigDifNM1\Bd);
argumentos.SigmaMi = SigmaMi;
fmipd = @(t,mid)(GamaMi*mid + SigmaMi*ydTil(t));
repositorio.fmipd = fmipd;
fmippd = @(t,mipd)(GamaMi*mipd + SigmaMi*ypdTil(t));
repositorio.fmippd = fmippd;

%Camada limite dinâmica:
%Obs.: DfiBd(fipBdTil)*fipBdTil + Kfi*fiBdTil - Klcd(t,mid) = 0
%(1) Expressões diretas:
D1fiFN = diag([1/(1+D11fn), 1/(1+D22fn)]); 
D2fiFN = diag([1/(1+D11fn), 1/(1-D22fn)]);
D3fiFN = diag([1/(1-D11fn), 1/(1+D22fn)]);
D4fiFN = diag([1/(1-D11fn), 1/(1-D22fn)]);
DfiBdFN = @(fipBdTil)(...
    D1fiFN*double((Id2(1,:)*fipBdTil >= 0) && (Id2(2,:)*fipBdTil >= 0)) +...
    D2fiFN*double((Id2(1,:)*fipBdTil >= 0) && (Id2(2,:)*fipBdTil < 0)) +...
    D3fiFN*double((Id2(1,:)*fipBdTil < 0) && (Id2(2,:)*fipBdTil >= 0)) +...
    D4fiFN*double((Id2(1,:)*fipBdTil < 0) && (Id2(2,:)*fipBdTil < 0)));
repositorio.DfiBdFN = DfiBdFN;
GamaFiBdFN = @(fipBdTil)(DfiBdFN(fipBdTil)*fipBdTil);
repositorio.GamaFiBdFN = GamaFiBdFN;
%(2) Derivadas:
dGamaFiFNdfipBdT = @(fipBdTil)(DfiBdFN(fipBdTil));
repositorio.dGamaFiFNdfipBdT = dGamaFiFNdfipBdT;
ffippBdFN = @(t,mid,mipd,fipBd)(DfiBdFN(fipBd)\(-Kfi*fipBd + DKlcdDt(t,mid,mipd)));
repositorio.ffippBdFN = ffippBdFN;

%%
simulacao = argumentos.simulacao;
switch simulacao
    case 1
        M = M1restr;
        argumentos.M = M;
        K = KbRestr;
        argumentos.K = K;
        C = Crestr;
        argumentos.C = C;
        B = Br;
        argumentos.B = B;
        W0 = Wsbr;
        argumentos.W0 = W0;
        %Atrito contínuo na origem:
        fTilc = @(ubr,ubpr)(SigmaR.'*fcEFbAtr(SigmaR*ubr,SigmaR*ubpr));
        repositorio.fTilc = fTilc;
        dfTilcdubT = @(ubr,ubpr)(SigmaR.'*dfcEFdubTAtr(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        repositorio.dfTilcdubT = dfTilcdubT;
        dfTilcdubpT = @(ubr,ubpr)(SigmaR.'*dfcEFdubpTAtr(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        repositorio.dfTilcdubpT = dfTilcdubpT;
        %Atrito descontínuo na origem:
        fTildc = @(ubr,ubpr)(SigmaR.'*fdcEFbAtr(SigmaR*ubr,SigmaR*ubpr));
        repositorio.fTildc = fTildc;
        dfTildcdubT = @(ubr,ubpr)(SigmaR.'*dfdcEFdubTAtr(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        repositorio.dfTildcdubT = dfTildcdubT;
        dfTildcdubpT = @(ubr,ubpr)(SigmaR.'*dfdcEFdubpTAtr(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        repositorio.dfTildcdubpT = dfTildcdubpT;
        %Entradas:
        %(a) Movimento de torção:
        R3 = argumentos.R3;
        nu = 1;
        %mic0 = 0.6;
        mic0 = 0;
        mi = mic0*fiAdm(nu, lambdaNM1miC, nNM1miC);
        Ttx0 = mi*M1*g*R3;
        fTtxRef = @(t)(1.5*Ttx0);
        %(b) Movimento longitudinal:
        %Obs.: 
        %(1) Função descontínua na origem
        %(2) Sinal contrário ao da velocidade u1bp
        %{
        alfaF = 1.2;
        dF = alfaF*(M3pl + M2b + dm + M1)*g;
        fFuRef = @(u1bp)(-dF*double(u1bp>0) + dF*double(u1bp<0));
        fUc = @(t,u1bp)([fTtxRef(t); fFuRef(u1bp)]);
        %}
        fFuRef = @(t)(0);
        fUc = @(t)([fTtxRef(t); fFuRef(t)]);
        repositorio.fUc = fUc;
        %Obs.: 
        %indG = 1: inclusão da função G em fc[d]EFaprox
        %indG = 0: retirada da função G em fc[d]EFaprox
    case 2
        M = eye(size(Krk));
        argumentos.M = M;
        K = Krk;
        argumentos.K = K;
        C = Crk;
        argumentos.C = C;
        B = Brk;
        argumentos.B = B;
        W0 = Pk.'*Wsbr;
        argumentos.W0 = W0;
        %Atrito contínuo na origem:
        fTilc = @(ubr,ubpr)(SigmaR.'*fcEFbAtr(SigmaR*ubr,SigmaR*ubpr));
        fTilck = @(ubrk,ubprk)(Pk.'*fTilc(Pk*ubrk,Pk*ubprk));
        repositorio.fTilc = fTilck;
        dfTilcdubT = @(ubr,ubpr)(SigmaR.'*dfcEFdubTAtr(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        dfTilckdukT = @(ubrk,ubprk)(Pk.'*dfTilcdubT(Pk*ubrk,Pk*ubprk)*Pk);
        repositorio.dfTilcduT = dfTilckdukT;
        dfTilcdubpT = @(ubr,ubpr)(SigmaR.'*dfcEFdubpTAtr(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        dfTilckdupkT = @(ubrk,ubprk)(Pk.'*dfTilcdubpT(Pk*ubrk,Pk*ubprk)*Pk);
        repositorio.dfTilcdupT = dfTilckdupkT;
        %Atrito descontínuo na origem:
        fTildc = @(ubr,ubpr)(SigmaR.'*fdcEFbAtr(SigmaR*ubr,SigmaR*ubpr));
        fTildck = @(ubrk,ubprk)(Pk.'*fTildc(Pk*ubrk,Pk*ubprk));
        repositorio.fTildc = fTildck;
        dfTildcdubT = @(ubr,ubpr)(SigmaR.'*dfdcEFdubTAtr(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        dfTildckdukT = @(ubrk,ubprk)(Pk.'*dfTildcdubT(Pk*ubrk,Pk*ubprk)*Pk);
        repositorio.dfTildcduT = dfTildckdukT;
        dfTildcdubpT = @(ubr,ubpr)(SigmaR.'*dfdcEFdubpTAtr(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        dfTildckdupkT = @(ubrk,ubprk)(Pk.'*dfTildcdubpT(Pk*ubrk,Pk*ubprk)*Pk);
        repositorio.dfTildcdupT = dfTildckdupkT;
        %Entradas:
        %(a) Movimento de torção:
        R3 = argumentos.R3;
        nu = 1;
        %mic0 = 0.6;
        mic0 = 0;
        mi = mic0*fiAdm(nu, lambdaNM1miC, nNM1miC);
        Ttx0 = mi*M1*g*R3;
        fTtxRef = @(t)(1.5*Ttx0);
        %(b) Movimento longitudinal:
        %Obs.: 
        %(1) Função descontínua na origem
        %(2) Sinal contrário ao da velocidade u1bp
        %{
        alfaF = 1.2;
        dF = alfaF*(M3pl + M2b + dm + M1)*g;
        fFuRef = @(u1bp)(-dF*double(u1bp>0) + dF*double(u1bp<0));
        fUc = @(t,u1bp)([fTtxRef(t); fFuRef(u1bp)]);
        %}
        fFuRef = @(t)(0);
        fUc = @(t)([fTtxRef(t); fFuRef(t)]);
        repositorio.fUc = fUc;
        %Obs.: 
        %indG = 1: inclusão da função G em fc[d]EFaprox
        %indG = 0: retirada da função G em fc[d]EFaprox
    case 3
        M = M1restr;
        argumentos.M = M;
        K = Krestr;
        argumentos.K = K;
        C = Crestr;
        argumentos.C = C;
        fTilc = @(ubr,ubpr,ubppr)(SigmaR.'*fcEFb(SigmaR*ubr,SigmaR*ubpr,SigmaR*ubppr));
        repositorio.fTilc = fTilc;
        B = Br;
        argumentos.B = B;
        W0 = Wsbr;
        argumentos.W0 = W0;
        %Atrito contínuo na origem:
        dfTilcdubT = @(ubr,ubpr,ubppr)(SigmaR.'*dfcEFbdubT(SigmaR*ubr,SigmaR*ubpr,SigmaR*ubppr)*SigmaR);
        repositorio.dfTilcdubT = dfTilcdubT;
        dfTilcdubpT = @(ubr,ubpr)(SigmaR.'*dfcEFbdubpT(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        repositorio.dfTilcdubpT = dfTilcdubpT;
        %Atrito descontínuo na origem:
        dfTildcdubT = @(ubr,ubpr,ubppr)(SigmaR.'*dfdcEFbdubT(SigmaR*ubr,SigmaR*ubpr,SigmaR*ubppr)*SigmaR);
        repositorio.dfTildcdubT = dfTildcdubT;
        dfTildcdubpT = @(ubr,ubpr)(SigmaR.'*dfdcEFbdubpT(SigmaR*ubr,SigmaR*ubpr)*SigmaR);
        repositorio.dfTildcdubpT = dfTildcdubpT;
        %Atrito contínuo/descontínuo na origem:
        dfTildubppT = @(ubr)(SigmaR.'*dfEFbdubppT(SigmaR*ubr)*SigmaR);
        repositorio.dfTildubppT = dfTildubppT;
        %Entradas:
        %(a) Movimento de torção:
        R3 = argumentos.R3;
        nu = 1;
        %mic0 = 0.6;
        mic0 = 0;
        mi = mic0*fiAdm(nu, lambdaNM1miC, nNM1miC);
        Ttx0 = mi*M1*g*R3;
        fTtxRef = @(t)(1.5*Ttx0);
        %(b) Movimento longitudinal:
        %Obs.: 
        %(1) Função descontínua na origem
        %(2) Sinal contrário ao da velocidade u1bp
        %{
        alfaF = 1.2;
        dF = alfaF*(M3pl + M2b + dm + M1)*g;
        fFuRef = @(u1bp)(-dF*double(u1bp>0) + dF*double(u1bp<0));
        fUc = @(t,u1bp)([fTtxRef(t); fFuRef(u1bp)]);
        %}
        fFuRef = @(t)(0);
        fUc = @(t)([fTtxRef(t); fFuRef(t)]);
        repositorio.fUc = fUc;
    case 4
        M = M1gEF;
        argumentos.M = M;
        K = KgEF;
        argumentos.K = K;
        C = CgEF;
        argumentos.C = C;
        fTilc = @(ub,ubp,ubpp)(fcEFb(ub,ubp,ubpp));
        repositorio.fTilc = fTilc;
        B = Bef;
        argumentos.B = B;
        W0 = Wsb;
        argumentos.W0 = W0;
        %Atrito contínuo na origem:
        dfTilcdubT = @(ub,ubp,ubpp)(dfcEFbdubT(ub,ubp,ubpp));
        repositorio.dfTilcdubT = dfTilcdubT;
        dfTilcdubpT = @(ub,ubp)(dfcEFbdubpT(ub,ubp));
        repositorio.dfTilcdubpT = dfTilcdubpT;
        %Atrito descontínuo na origem:
        dfTildcdubT = @(ub,ubp,ubpp)(dfdcEFbdubT(ub,ubp,ubpp));
        repositorio.dfTildcdubT = dfTildcdubT;
        dfTildcdubpT = @(ub,ubp)(dfdcEFbdubpT(ub,ubp));
        repositorio.dfTildcdubpT = dfTildcdubpT;
        %Atrito contínuo/descontínuo na origem:
        dfTildubppT = @(ub)(dfEFbdubppT(ub));
        repositorio.dfTildubppT = dfTildubppT;
        %Entradas:
        %(a) Movimento de torção:
        R3 = argumentos.R3;
        nu = 1;
        %mic0 = 0.6;
        mic0 = 0;
        mi = mic0*fiAdm(nu, lambdaNM1miC, nNM1miC);
        Ttx0 = mi*M1*g*R3;
        fTtxRef = @(t)(1.5*Ttx0);
        %(b) Movimento longitudinal:
        %Obs.: 
        %(1) Função descontínua na origem
        %(2) Sinal contrário ao da velocidade u1bp
        %{
        alfaF = 1.2;
        dF = alfaF*(M3pl + M2b + dm + M1)*g;
        fFuRef = @(u1bp)(-dF*double(u1bp>0) + dF*double(u1bp<0));
        fUc = @(t,u1bp)([fTtxRef(t); fFuRef(u1bp)]);
        %}
        fFuRef = @(t)(0);
        fUc = @(t)([fTtxRef(t); fFuRef(t)]);
        repositorio.fUc = fUc;
    case 5
    case 6
    case 7
    case 8
    case 9
    otherwise
        disp('other value')
end
%%
%Limitação de variáveis cíclicas:
fLVC = @(q,condicao,IdLVC)(q - 2*pi*IdLVC*double(condicao));
repositorio.fLVC = fLVC;
%}
end