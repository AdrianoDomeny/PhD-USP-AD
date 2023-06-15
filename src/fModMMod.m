function [ PtetaX, sigTiltetaX, Pu, sigTilu, Pvw, sigTilvw, Abvt, Abvl, info ] =...
    fModMMod( MefRestr, KefRestr, argumentos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Gráficos dos vetores unitários:
graf = 0;
%Obs.: se graf>0, então serão plotados os gráficos. Caso contrário, não.

pond = 'p';
%Obs.:
%(1) pond = 'p': normalização ponderada dos autovetores
%(2) pond = 'np': normalização não ponderada

%Garantia de que as matrizes de inércia e rigidez são simétricas:
KefRestr = (1/2)*(KefRestr + KefRestr.');
MefRestr = (1/2)*(MefRestr + MefRestr.');
%{
dK = (1/2)*(KefRestr - KefRestr.');
dM = (1/2)*(MefRestr - MefRestr.');
ndK = norm(dK);
ndM = norm(dM);
%Obs.:
%(1) dK e dM são nulas para N = 100
%}
switch pond
    case 'p'
        [ Avt, Avl ] = eig(KefRestr,MefRestr);
    case 'np'
        [ Avt, Avl ] = eig(MefRestr\KefRestr);
    otherwise
        disp('other value')
end
Avl = diag(Avl);
[ Abvl, Ivl] = sort(Avl);
nAvl = length(Avl);
Id = eye(nAvl);
Sigma = zeros(nAvl);
for i = 1:nAvl
    Sigma(i,:) = Id(Ivl(i),:);
end
Abvt = Avt*Sigma.';

%Modos de corpo rígido (existência):
sigTol = 10^-2;
IdnModCR = double(Abvl<sigTol*ones(nAvl,1));
%Obs.: IdnModCR indica a posição dos modos de corpo rígido em Abvt
nModCR = sum(IdnModCR);
%Obs.: nModCR indica a número de modos de corpo rígido do sistema
if nModCR>1
    Abvl(1,1) = 0;
    Abvl(2,1) = 0;
end

%Modos de vibração em sistema com parâmetros concentrados:
Ip1 = argumentos.Ip1;
M1 = argumentos.M1rt;
Ip2b = argumentos.Ip2b;
M2b = argumentos.M2b;
dm = argumentos.dm;
l2 = argumentos.l2;
M2 = M2b + dm;
Ibp2 = Ip2b + dm*l2^2;
Ibp3 = argumentos.Ip3bMt;
M3 = argumentos.M3_pl;
K10 = argumentos.K10;
K20 = argumentos.K20;
k10 = argumentos.k10;
k20 = argumentos.k20;
%(1) Sistemas com três graus de liberdade, e modo de corpo rígido em tetaX:
I12 = (Ip1^-1 + Ibp2^-1)^-1;
I23 = (Ibp2^-1 + Ibp3^-1)^-1;
I123 = (1/(Ibp3*Ip1) + 1/(Ibp2*Ip1) + 1/(Ibp2*Ibp3))^-1;
b3Tx = K20/I12 + K10/I23;
c3Tx = K10*K20/I123;
delta3Tx = b3Tx^2 - 4*c3Tx;
%(1a) Frequências naturais:
lmb32tX = (1/2)*(b3Tx - sqrt(delta3Tx));
lmb33tX = (1/2)*(b3Tx + sqrt(delta3Tx));
%(1b) Modos de vibração:
bdTx = K10/Ibp3 + (K10+K20)/Ibp2;
cdTx = K10*K20/(Ibp2*Ibp3);
dC3tX = @(lmb)(lmb^2 - bdTx*lmb + cdTx);
P32tX = [K10/Ibp3, K10/Ibp3 - lmb32tX, (Ibp2/K20)*dC3tX(lmb32tX)].';
P33tX = [K10/Ibp3, K10/Ibp3 - lmb33tX, (Ibp2/K20)*dC3tX(lmb33tX)].';
%(2) Sistemas com dois graus de liberdade, sem modo de CR em tetaX:
b2Tx = K10/Ibp2 + K20/I12;
c2Tx = K10*K20/(Ibp2*Ip1);
delta2Tx = b2Tx^2 - 4*c2Tx;
%(2a) Frequências naturais:
lmb21 = (1/2)*(b2Tx - sqrt(delta2Tx));
lmb22 = (1/2)*(b2Tx + sqrt(delta2Tx));
%(2b) Modos de vibração:
dC2tX = @(lmb)((K10 + K20)/Ibp2 - lmb);
P21tX = [(K20/Ibp2), dC2tX(lmb21)].';
P22tX = [(K20/Ibp2), dC2tX(lmb22)].';

%(3) Sistemas com três graus de liberdade, e modo de corpo rígido em u:
M12 = (M1^-1 + M2^-1)^-1;
M23 = (M2^-1 + M3^-1)^-1;
M123 = (1/(M3*M1) + 1/(M2*M1) + 1/(M2*M3))^-1;
b3u = k20/M12 + k10/M23;
c3u = k10*k20/M123;
delta3u = b3u^2 - 4*c3u;
%(1a) Frequências naturais:
lmb32u = (1/2)*(b3u - sqrt(delta3u));
lmb33u = (1/2)*(b3u + sqrt(delta3u));
%(1b) Modos de vibração:
bdu = k10/M3 + (k10+k20)/M2;
cdu = k10*k20/(M2*M3);
dC3u = @(lmb)(lmb^2 - bdu*lmb + cdu);
P32u = [k10/M3, k10/M3 - lmb32u, (M2/k20)*dC3u(lmb32u)].';
P33u = [k10/M3, k10/M3 - lmb33u, (M2/k20)*dC3u(lmb33u)].';
%(2) Sistemas com dois graus de liberdade, sem modo de CR em tetaX:
b2u = k10/M2 + k20/M12;
c2u = k10*k20/(M2*M1);
delta2u = b2u^2 - 4*c2u;
%(2a) Frequências naturais:
lmb21 = (1/2)*(b2u - sqrt(delta2u));
lmb22 = (1/2)*(b2u + sqrt(delta2u));
%(2b) Modos de vibração:
dC2u = @(lmb)((k10 + k20)/M2 - lmb);
P21u = [(k20/M2), dC2u(lmb21)].';
P22u = [(k20/M2), dC2u(lmb22)].';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch pond
    case 'p'
        AbvtNorm2 = Abvt.'*MefRestr*Abvt;
        %Obs.: AvtNorm2:nXn
        AbvtNorm = sqrt(abs(diag(AbvtNorm2)));
        %Obs.: AvtNorm:nX1
        Abvt = Abvt*diag(AbvtNorm.^-1);
        %Obs.: As colunas de Abvt estão normalizadas segundo MefRestr
    case 'np'
        AbvtNorm2 = Avt.'*Avt;
        %Obs.: AvtNorm2:nXn
        AbvtNorm = sqrt(abs(diag(AbvtNorm2)));
        %Obs.: AvtNorm:nX1
        Abvt = Abvt*diag(AbvtNorm.^-1);
        %Obs.: As colunas de Abvt estão normalizadas segundo Id
    otherwise
        disp('other value')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
L = argumentos.L;
L1 = argumentos.L1;
N = argumentos.N;
N1 = argumentos.N1;
X = argumentos.X;
SigmaR = argumentos.SigmaReF;

%Modos torcionais:
%(1) Modo de corpo rígido:
vTilTetaX0 = zeros(6*(N+1),1);
vTx0 = zeros(N+1,1);
%(2) Modo 1:
vTilTetaX1 = zeros(6*(N+1),1);
vTx1 = zeros(N+1,1);
%(3) Modo 2:
vTilTetaX2 = zeros(6*(N+1),1);
vTx2 = zeros(N+1,1);
for i = 1:N+1
    if nModCR>0
        vTilTetaX1(6*1) = P32tX(1,1);
        vTilTetaX1(6*(N1+1)) = P32tX(2,1);
        vTilTetaX1(6*(N+1)) = P32tX(3,1);
        vTilTetaX2(6*1) = P33tX(1,1);
        vTilTetaX2(6*(N1+1)) = P33tX(2,1);
        vTilTetaX2(6*(N+1)) = P33tX(3,1);
        if i<N1+2
            n = 6*i;
            vTilTetaX0(n) = 1;
            vTx0(i) = vTilTetaX0(n);
            vTilTetaX1(n) = vTilTetaX1(6*1) + ((vTilTetaX1(6*(N1+1)) - vTilTetaX1(6*1))/L1)*(X(i) - X(1));
            vTx1(i) = vTilTetaX1(n);
            vTilTetaX2(n) = vTilTetaX2(6*1) + ((vTilTetaX2(6*(N1+1)) - vTilTetaX2(6*1))/L1)*(X(i) - X(1));
            vTx2(i) = vTilTetaX2(n);
        else
            n = 6*i;
            vTilTetaX0(n) = 1;
            vTx0(i) = vTilTetaX0(n);
            vTilTetaX1(n) = vTilTetaX1(6*(N1+1)) +...
                ((vTilTetaX1(6*(N+1)) - vTilTetaX1(6*(N1+1)))/(L - L1))*(X(i) - X(N1+1));
            vTx1(i) = vTilTetaX1(n);
            vTilTetaX2(n) = vTilTetaX2(6*(N1+1)) +...
                ((vTilTetaX2(6*(N+1)) - vTilTetaX2(6*(N1+1)))/(L - L1))*(X(i) - X(N1+1));
            vTx2(i) = vTilTetaX2(n);
        end
    else
        vTilTetaX1(6*1) = 0;
        vTilTetaX1(6*(N1+1)) = P21tX(1,1);
        vTilTetaX1(6*(N+1)) = P21tX(2,1);
        vTilTetaX2(6*1) = 0;
        vTilTetaX2(6*(N1+1)) = P22tX(1,1);
        vTilTetaX2(6*(N+1)) = P22tX(2,1);
        if i<N1+2
            n = 6*i;
            vTilTetaX1(n) = vTilTetaX1(6*1) +...
                ((vTilTetaX1(6*(N1+1)) - vTilTetaX1(6*1))/L1)*(X(i) - X(1));
            vTx1(i) = vTilTetaX1(n);
            vTilTetaX2(n) = vTilTetaX2(6*1) +...
                ((vTilTetaX2(6*(N1+1)) - vTilTetaX2(6*1))/L1)*(X(i) - X(1));
            vTx2(i) = vTilTetaX2(n);
        else
            n = 6*i;
            vTilTetaX1(n) = vTilTetaX1(6*(N1+1)) +...
                ((vTilTetaX1(6*(N+1)) - vTilTetaX1(6*(N1+1)))/(L - L1))*(X(i) - X(N1+1));
            vTx1(i) = vTilTetaX1(n);
            vTilTetaX2(n) = vTilTetaX2(6*(N1+1)) +...
                ((vTilTetaX2(6*(N+1)) - vTilTetaX2(6*(N1+1)))/(L - L1))*(X(i) - X(N1+1));
            vTx2(i) = vTilTetaX2(n);
        end
    end
end

%Modos longitudinais:
%(1) Modo de corpo rígido:
vTilu0 = zeros(6*(N+1),1);
vU0 = zeros(N+1,1);
%(2) Modo 1:
vTilu1 = zeros(6*(N+1),1);
vU1 = zeros(N+1,1);
%(3) Modo 2:
vTilu2 = zeros(6*(N+1),1);
vU2 = zeros(N+1,1);
for i = 1:N+1
    if nModCR>0
        vTilu1(6*1-5) = P32u(1,1);
        vTilu1(6*(N1+1)-5) = P32u(2,1);
        vTilu1(6*(N+1)-5) = P32u(3,1);
        vTilu2(6*1-5) = P33u(1,1);
        vTilu2(6*(N1+1)-5) = P33u(2,1);
        vTilu2(6*(N+1)-5) = P33u(3,1);
        if i<N1+2
            n = 6*i-5;
            vTilu0(n) = 1;
            vU0(i) = vTilu0(n);
            vTilu1(n) = vTilu1(6*1-5) + ((vTilu1(6*(N1+1)-5) - vTilu1(6*1-5))/L1)*(X(i) - X(1));
            vU1(i) = vTilu1(n);
            vTilu2(n) = vTilu2(6*1-5) + ((vTilu2(6*(N1+1)-5) - vTilu2(6*1-5))/L1)*(X(i) - X(1));
            vU2(i) = vTilu2(n);
        else
            n = 6*i-5;
            vTilu0(n) = 1;
            vU0(i) = vTilu0(n);
            vTilu1(n) = vTilu1(6*(N1+1)-5) + ((vTilu1(6*(N+1)-5) - vTilu1(6*(N1+1)-5))/(L - L1))*(X(i) - X(N1+1));
            vU1(i) = vTilu1(n);
            vTilu2(n) = vTilu2(6*(N1+1)-5) + ((vTilu2(6*(N+1)-5) - vTilu2(6*(N1+1)-5))/(L - L1))*(X(i) - X(N1+1));
            vU2(i) = vTilu2(n);
        end
    else
        vTilu1(6*1-5) = 0;
        vTilu1(6*(N1+1)-5) = P21u(1,1);
        vTilu1(6*(N+1)-5) = P21u(2,1);
        vTilu2(6*1-5) = 0;
        vTilu2(6*(N1+1)-5) = P22u(1,1);
        vTilu2(6*(N+1)-5) = P22u(2,1);
        if i<N1+2
            n = 6*i-5;
            vTilu1(n) = vTilu1(6*1-5) + ((vTilu1(6*(N1+1)-5) - vTilu1(6*1-5))/L1)*(X(i) - X(1));
            vU1(i) = vTilu1(n);
            vTilu2(n) = vTilu2(6*1-5) + ((vTilu2(6*(N1+1)-5) - vTilu2(6*1-5))/L1)*(X(i) - X(1));
            vU2(i) = vTilu2(n);
        else
            n = 6*i-5;
            vTilu1(n) = vTilu1(6*(N1+1)-5) + ((vTilu1(6*(N+1)-5) - vTilu1(6*(N1+1)-5))/(L - L1))*(X(i) - X(N1+1));
            vU1(i) = vTilu1(n);
            vTilu2(n) = vTilu2(6*(N1+1)-5) + ((vTilu2(6*(N+1)-5) - vTilu2(6*(N1+1)-5))/(L - L1))*(X(i) - X(N1+1));
            vU2(i) = vTilu2(n);
        end
    end
end

%
%Modos intermediários, com frequências próximas aos modos torsionais:
%(1) Modo v:
vTilv = zeros(6*(N+1),1);
vV = zeros(N+1,1);
%(2) Modo w:
vTilw = zeros(6*(N+1),1);
vW = zeros(N+1,1);
for i = 1:N+1
    vTilv(6*1-4) = 0;
    vTilv(6*(N1+1)-4) = 1;
    vTilv(6*(N+1)-4) = 0;
    vTilw(6*1-2) = 0;
    vTilw(6*(N1+1)-2) = 1;
    vTilw(6*(N+1)-2) = 0;
    if i<N1+2
        nV = 6*i-4;
        vTilv(nV) = vTilv(6*1-4) + ((vTilv(6*(N1+1)-4) - vTilv(6*1-4))/L1)*(X(i) - X(1));
        vV(i) = vTilv(nV);
        nW = 6*i-2;
        vTilw(nW) = vTilw(6*1-2) + ((vTilw(6*(N1+1)-2) - vTilw(6*1-2))/L1)*(X(i) - X(1));
        vW(i) = vTilw(nW);
    else
        nV = 6*i-4;
        vTilv(nV) = vTilv(6*(N1+1)-4) +...
            ((vTilv(6*(N+1)-4) - vTilv(6*(N1+1)-4))/(L - L1))*(X(i) - X(N1+1));
        vV(i) = vTilv(nV);
        nW = 6*i-2;
        vTilw(nW) = vTilw(6*(N1+1)-2) +...
            ((vTilw(6*(N+1)-2) - vTilw(6*(N1+1)-2))/(L - L1))*(X(i) - X(N1+1));
        vW(i) = vTilw(nW);
    end
end
%}

%%
%Definição da norma para os autovetores:
%Obs.:
%(1) AsimPd é matriz quadrada simétrica positiva definida
%(2) Se v:nX1, então AsimPd:nXn
normV = @(v,AsimPd)(sqrt(v.'*AsimPd*v));

%Modos torcionais:
switch pond
    case 'p'
        if nModCR>0
            vTilTx0 = SigmaR.'*vTilTetaX0;
            normVtX0 = normV(vTilTx0,MefRestr);
            vTilTx0 = vTilTx0/normVtX0;
            [ PtX0, jPtX0, MacTx0 ] = fModPxj( Abvt, vTilTx0, MefRestr );
            info.jPtX0 = jPtX0;
            info.tTx0 = MacTx0;
            vTilTx1 = SigmaR.'*vTilTetaX1;
            normVtX1 = normV(vTilTx1,MefRestr);
            vTilTx1 = vTilTx1/normVtX1;
            [ PtX1, jPtX1, MacTx1 ] = fModPxj( Abvt, vTilTx1, MefRestr );
            vTilTx2 = SigmaR.'*vTilTetaX2;
            normVtX2 = normV(vTilTx2,MefRestr);
            vTilTx2 = vTilTx2/normVtX2;
            [ PtX2, jPtX2, MacTx2 ] = fModPxj( Abvt, vTilTx2, MefRestr );
        else
            vTilTx1 = SigmaR.'*vTilTetaX1;
            normVtX1 = normV(vTilTx1,MefRestr);
            vTilTx1 = vTilTx1/normVtX1;
            [ PtX1, jPtX1, MacTx1 ] = fModPxj( Abvt, vTilTx1, MefRestr );
            vTilTx2 = SigmaR.'*vTilTetaX2;
            normVtX2 = normV(vTilTx2,MefRestr);
            vTilTx2 = vTilTx2/normVtX2;
            [ PtX2, jPtX2, MacTx2 ] = fModPxj( Abvt, vTilTx2, MefRestr );
        end
    case 'np'
        if nModCR>0
            vTilTx0 = SigmaR.'*vTilTetaX0;
            normVtX0 = normV(vTilTx0,Id);
            vTilTx0 = vTilTx0/normVtX0;
            [ PtX0, jPtX0, MacTx0 ] = fModPxj( Abvt, vTilTx0, Id );
            info.jPtX0 = jPtX0;
            info.tTx0 = MacTx0;
            vTilTx1 = SigmaR.'*vTilTetaX1;
            normVtX1 = normV(vTilTx1,Id);
            vTilTx1 = vTilTx1/normVtX1;
            [ PtX1, jPtX1, MacTx1 ] = fModPxj( Abvt, vTilTx1, Id );
            vTilTx2 = SigmaR.'*vTilTetaX2;
            normVtX2 = normV(vTilTx2,Id);
            vTilTx2 = vTilTx2/normVtX2;
            [ PtX2, jPtX2, MacTx2 ] = fModPxj( Abvt, vTilTx2, Id );
        else
            vTilTx1 = SigmaR.'*vTilTetaX1;
            normVtX1 = normV(vTilTx1,Id);
            vTilTx1 = vTilTx1/normVtX1;
            [ PtX1, jPtX1, MacTx1 ] = fModPxj( Abvt, vTilTx1, Id );
            vTilTx2 = SigmaR.'*vTilTetaX2;
            normVtX2 = normV(vTilTx2,Id);
            vTilTx2 = vTilTx2/normVtX2;
            [ PtX2, jPtX2, MacTx2 ] = fModPxj( Abvt, vTilTx2, Id );
        end
    otherwise
        disp('other value')
end
if nModCR>0
    vTx0 = vTx0/normVtX0;
    sigTx0 = Abvl(jPtX0);
end
vTx1 = vTx1/normVtX1;
vTx2 = vTx2/normVtX2;
sigTx1 = Abvl(jPtX1);
sigTx2 = Abvl(jPtX2);
if graf>0
    if nModCR>0
        figure
        plot(X,vTx0,'o-',X,vTx1,'.-',X,vTx2,'x-')
        xlabel('X (m)')
        ylabel('teta_{X}')
        %title('Modos de torção')
        legend('Modo de corpo rígido','1o modo','2o modo')
        grid
    else
        figure
        plot(X,vTx1,'o-',X,vTx2,'o-')
        xlabel('X')
        ylabel('vTx')
        legend('vTx1','vTx2')
        grid
    end
end

%Matrizes com modos e frequências (w^2) torcionais:
if nModCR>0
    PtetaX = [PtX0,PtX1,PtX2];
    PtetaX(:,1) = vTilTx0;
    Abvt(:,jPtX0) = vTilTx0;
    sigTiltetaX = [sigTx0; sigTx1; sigTx2];
    sigTiltetaX(1,1) = 0;
else
    PtetaX = [PtX1,PtX2];
    sigTiltetaX = [sigTx1; sigTx2];
end

%%
%Modos longitudinais:
switch pond
    case 'p'
        if nModCR>0
            vTilU0 = SigmaR.'*vTilu0;
            normVu0 = normV(vTilU0,MefRestr);
            vTilU0 = vTilU0/normVu0;
            [ Pu0, jPu0, MacU0 ] = fModPxj( Abvt, vTilU0, MefRestr );
            info.jPu0 = jPu0;
            info.tU0 = MacU0;
            vTilU1 = SigmaR.'*vTilu1;
            normVu1 = normV(vTilU1,MefRestr);
            vTilU1 = vTilU1/normVu1;
            [ Pu1, jPu1, MacU1 ] = fModPxj( Abvt, vTilU1, MefRestr );
            vTilU2 = SigmaR.'*vTilu2;
            normVu2 = normV(vTilU2,MefRestr);
            vTilU2 = vTilU2/normVu2;
            [ Pu2, jPu2, MacU2 ] = fModPxj( Abvt, vTilU2, MefRestr );
        else
            vTilU1 = SigmaR.'*vTilu1;
            normVu1 = normV(vTilU1,MefRestr);
            vTilU1 = vTilU1/normVu1;
            [ Pu1, jPu1, MacU1 ] = fModPxj( Abvt, vTilU1, MefRestr );
            vTilU2 = SigmaR.'*vTilu2;
            normVu2 = normV(vTilU2,MefRestr);
            vTilU2 = vTilU2/normVu2;
            [ Pu2, jPu2, MacU2 ] = fModPxj( Abvt, vTilU2, MefRestr );
        end
    case 'np'
        if nModCR>0
            vTilU0 = SigmaR.'*vTilu0;
            normVu0 = normV(vTilU0,Id);
            vTilU0 = vTilU0/normVu0;
            [ Pu0, jPu0, MacU0 ] = fModPxj( Abvt, vTilU0, Id );
            info.jPu0 = jPu0;
            info.tU0 = MacU0;
            vTilU1 = SigmaR.'*vTilu1;
            normVu1 = normV(vTilU1,Id);
            vTilU1 = vTilU1/normVu1;
            [ Pu1, jPu1, MacU1 ] = fModPxj( Abvt, vTilU1, Id );
            vTilU2 = SigmaR.'*vTilu2;
            normVu2 = normV(vTilU2,Id);
            vTilU2 = vTilU2/normVu2;
            [ Pu2, jPtX2, MacTx2 ] = fModPxj( Abvt, vTilU2, Id );
        else
            vTilU1 = SigmaR.'*vTilu1;
            normVu1 = normV(vTilU1,Id);
            vTilU1 = vTilU1/normVu1;
            [ Pu1, jPu1, MacU1 ] = fModPxj( Abvt, vTilU1, Id );
            vTilU2 = SigmaR.'*vTilu2;
            normVu2 = normV(vTilU2,Id);
            vTilU2 = vTilU2/normVu2;
            [ Pu2, jPtX2, MacTx2 ] = fModPxj( Abvt, vTilU2, Id );
        end
    otherwise
        disp('other value')
end
if nModCR>0
    vU0 = vU0/normVu0;
    sigU0 = Abvl(jPu0);
end
vU1 = vU1/normVu1;
vU2 = vU2/normVu2;
sigU1 = Abvl(jPu1);
sigU2 = Abvl(jPu2);
if graf>0
    if nModCR>0
        figure
        plot(X,vU0,'o-',X,vU1,'.-',X,vU2,'x-')
        xlabel('X (m)')
        ylabel('u')
        %title('Modos de torção')
        legend('Modo de corpo rígido','1o modo','2o modo')
        grid
    else
        figure
        plot(X,vU1,'o-',X,vU2,'o-')
        xlabel('X')
        ylabel('vU')
        legend('vU1','vU2')
        grid
    end
end

if nModCR>0
    Pu = [Pu0,Pu1,Pu2];
    Pu(:,1) = vTilU0;
    Abvt(:,jPu0) = vTilU0;
    sigTilu = [sigU0; sigU1; sigU2];
    sigTilu(1,1) = 0;
else
    Pu = [Pu1,Pu2];
    sigTilu = [sigU1; sigU2];
end

%
%%
%Modos intermediários, com frequências próximas aos modos torsionais:
switch pond
    case 'p'
        vTilv = SigmaR.'*vTilv;
        normVv = normV(vTilv,MefRestr);
        vTilv = vTilv/normVv;
        [ Pv, jPv, MacV ] = fModPxj( Abvt, vTilv, MefRestr );
        vTilw = SigmaR.'*vTilw;
        normVw = normV(vTilw,MefRestr);
        vTilw = vTilw/normVw;
        [ Pw, jPw, MacW ] = fModPxj( Abvt, vTilw, MefRestr );
    case 'np'
        vTilv = SigmaR.'*vTilv;
        normVv = normV(vTilv,Id);
        vTilv = vTilv/normVv;
        [ Pv, jPv, MacV ] = fModPxj( Abvt, vTilv, Id );
        vTilw = SigmaR.'*vTilw;
        normVw = normV(vTilw,MefRestr);
        vTilw = vTilw/normVw;
        [ Pw, jPw, MacW ] = fModPxj( Abvt, vTilw, Id );
    otherwise
        disp('other value')
end
vV = vV/normVv;
vW = vW/normVw;
sigV = Abvl(jPv);
sigW = Abvl(jPw);
if graf>0
    figure
    plot(X,vV,'o-')
    xlabel('X (m)')
    ylabel('v')
    grid
    figure
    plot(X,vW,'o-')
    xlabel('X (m)')
    ylabel('w')
    grid
end

%Matrizes com modos e frequências (w^2) intermediárias aos modos torsionais:
Pvw = [Pv,Pw];
sigTilvw = [sigV; sigW];
%}

%{
%%
%Modo espiral:
vTilEsp = @(wP,fi0)(fModEsp( wP, fi0, argumentos ));
vTilEspr = @(wP,fi0)(SigmaR.'*vTilEsp(wP,fi0));
%{
QsiOtm = @(vOtm)(-vTilEspr([1,0]*vOtm,[0,1]*vOtm).'*(Abvt*Abvt.')*vTilEspr([1,0]*vOtm,[0,1]*vOtm));
vOtm0 = [1;0];
vOtm = fminsearch(QsiOtm,vOtm0);
wPotm = [1,0]*vOtm;
fi0otm = [0,1]*vOtm;
%Obs.: vOtm = [wP; fi0]
%}
wP0 = 2*pi/L;
%nWP = 1; MacEsp 0.99
nWP = 2;
wP = nWP*wP0;
fi0 = 0;
vTilEsprOtm = vTilEspr(wP,fi0);
switch pond
    case 'p'
        normVesp = normV(vTilEsprOtm,MefRestr);
        vTilEsprOtm = vTilEsprOtm/normVesp;
        [ Pesp, jPesp, MacEsp ] = fModPxj( Abvt, vTilEsprOtm, MefRestr );
    case 'np'
        normVesp = normV(vTilEsprOtm,Id);
        vTilEsprOtm = vTilEsprOtm/normVesp;
        [ Pesp, jPesp, MacEsp ] = fModPxj( Abvt, vTilEsprOtm, Id );
    otherwise
        disp('other value')
end
SigmaRMaV = argumentos.SigmaRMaV;
SigmaRMaW = argumentos.SigmaRMaW;
vTilvtEsp = SigmaR*vTilEsprOtm;
vEsp = SigmaRMaV.'*vTilvtEsp;
vEsp = vEsp/normVesp;
wEsp = SigmaRMaW.'*vTilvtEsp;
wEsp = wEsp/normVesp;
sigEsp = Abvl(jPesp);
if graf>0
    figure
    plot(X,vEsp,'o-')
    xlabel('X')
    ylabel('vV')
    grid
    figure
    plot(X,wEsp,'o-')
    xlabel('X')
    ylabel('vW')
    grid
end
N1 = argumentos.N1;
vAbvt = SigmaRMaV.'*SigmaR*Abvt;
wAbvt = SigmaRMaW.'*SigmaR*Abvt;
nEsp = floor(N1/2) - 3;
X6 = X(1,nEsp+1:nEsp+5);
v6 = vAbvt(1:N1,:);
w6 = wAbvt(1:N1,:);
LmbPa = zeros(N+1,6);
iF = 1;
ModEsp = zeros(1,5);
for i = 1:N+1
    ApA = [v6(:,i).^2, 2*v6(:,i).*w6(:,i), w6(:,i).^2, v6(:,i), w6(:,i)];
    bPA = -ones(N1,1);
    QsiPA = @(p)((1/2)*(ApA*p-bPA).'*(ApA*p-bPA));
    p0 = ones(5,1);
    [pA,fVal] = fminsearch(QsiPA,p0);
    a = pA(1,1);
    b = pA(2,1);
    c = pA(3,1);
    d = pA(4,1);
    e = pA(5,1);
    A = [a, b; b, c];
    [Pap,lmbAp] = eig(A);
    lmbAp = diag(lmbAp);
    pq = [d, e]*Pap;
    p = pq(1,1);
    q = pq(1,2);
    F = p^2/(4*lmbAp(1,1)) + q^2/(4*lmbAp(2,1)) - 1;
    LmbPa(i,1:6) = [lmbAp.',p,q,F,fVal];
    if (lmbAp(1,1)>0 && lmbAp(2,1)>0) && F>0
        ModEsp(iF,:) = [i,lmbAp.',F,fVal];
        iF = iF + 1;
    end
end
%}

%%
%Outros dados:
info.jPtX1 = jPtX1;
info.tTx1 = MacTx1;
info.jPtX2 = jPtX2;
info.tTx2 = MacTx2;
info.jPu1 = jPu1;
info.tU1 = MacU1;
info.jPu2 = jPu2;
info.tU2 = MacU2;
%
%Modos intermeriários entre os dois modos torsionais:
info.jPv = jPv;
info.tV = MacV;
info.jPw = jPw;
info.tW = MacW;
%}
%{
%Modo espiral:
info.jPesp = jPesp;
info.tEsp = MacEsp;
%}
end