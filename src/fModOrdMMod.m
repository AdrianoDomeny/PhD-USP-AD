function [ ABvEps, info ] = fModOrdMMod( argumentos )
%UNTITLED Summary of this function goes here
%   Ordem de grandeza das frequências(^2)
Ab = argumentos.Ab;
Bb = argumentos.Bb;
%Ng = argumentos.Ng;
%rEstr = argumentos.rEstr;
%Krk = argumentos.Krk;
%Crk = argumentos.Crk;
%Brk = argumentos.Brk;
%sigTilk = argumentos.sigTilk;
Pk = argumentos.Pk;
%Rs = argumentos.Rs;
%Ns = argumentos.Ns;
%Rf = argumentos.Rf;
%Nf = argumentos.Nf;
%jPtX0 = argumentos.jPtX0;
%jPu0 = argumentos.jPu0;
%jPtX1 = argumentos.jPtX1;
%jPtX2 = argumentos.jPtX2;
%jPu1 = argumentos.jPu1;
jPu2 = argumentos.jPu2;
%epsZero = argumentos.epsZero;
%nSigTilk = length(sigTilk);
%idSig = double(abs(sigTilk) < epsZero*10^4);
%ndSig0 = sum(idSig);

%%
nAb = length(Ab);
IdAb = eye(nAb);
sPk = size(Pk);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FORMA NORMAL SEM LQR (saídas ref. acoplamento)
%Obs.: 
%(1) Z é o vetor original de estado
%(2) Saídas: tetaXp1 e ub1 (ATENÇÃO!)

%Saídas: y1 = tetaXp1 e y2 = ub1 (ATENÇÃO!)
C11 = [zeros(1,sPk(1,2)),Pk(2,:)];
C21 = [Pk(1,:),zeros(1,sPk(1,2))];

%Dinâmica externa:
%Obs.: yTil = aTils(X) + bTilb*Us
beta01 = C11;
beta11 = beta01*Ab;
beta01Bb = beta01*Bb;
%{
beta2 = beta1*Ab;
beta1Bb = beta1*Bb;
beta3 = beta2*Ab;
beta2Bb = beta2*Bb;
beta4 = beta3*Ab;
beta3Bb = beta3*Bb;
beta5 = beta4*Ab;
beta4Bb = beta4*Bb;
beta6 = beta5*Ab;
beta5Bb = beta5*Bb;
beta7 = beta6*Ab;
beta6Bb = beta6*Bb;
beta8 = beta7*Ab;
beta7Bb = beta7*Bb;
beta9 = beta8*Ab;
beta8Bb = beta8*Bb;
beta10 = beta9*Ab;
beta9Bb = beta9*Bb;
beta11 = beta10*Ab;
beta10Bb = beta10*Bb;
%}
alfa01 = C21;
alfa11 = alfa01*Ab;
alfa21 = alfa11*Ab;
alfa11Bb = alfa11*Bb;
%{
alfa0Bbs = alfa0*Bbs;
alfa3 = alfa2*Abs;
alfa2Bbs = alfa2*Bbs;
alfa4 = alfa3*Abs;
alfa3Bbs = alfa3*Bbs;
ad = [beta0Bbs, beta1Bbs, beta2Bbs, beta3Bbs, beta4Bbs, beta5Bbs, beta6Bbs, beta7Bbs, beta8Bbs, beta9Bbs, beta10Bbs, alfa0Bbs, alfa1Bbs, alfa2Bbs, alfa3Bbs];
%}

%%
%Dinâmica interna:
%Obs.:
%(1) mi = Ami*Z
%(2) mip = GamaMi*mi + SigmaMi*yTil + Ami*Bb*U
[maxBb1,jM1] = max(abs(Bb(:,1)));
[maxBb2,jM2] = max(abs(Bb(:,2)));
info.maxBb11 = maxBb1;
info.maxBb21 = maxBb2;
bn11 = Bb(jM1,1);
bn22 = Bb(jM2,2);
bn21 = Bb(jM2,1);
bn12 = Bb(jM1,2);
Bn = [bn11, bn21; bn12, bn22];
invBn = Bn\eye(2);
ibn11 = invBn(1,1);
ibn22 = invBn(2,2);
ibn21 = invBn(2,1);
ibn12 = invBn(1,2);
AmiTotal = zeros(nAb-2,nAb);
j = 1;
for i=1:nAb
    if i ~= jM1 && i ~= jM2
        AmiTotal(j,:) = IdAb(i,:) -...
            IdAb(jM1,:)*(Bb(i,1)*ibn11 + Bb(i,2)*ibn21) -...
            IdAb(jM2,:)*(Bb(i,1)*ibn12 + Bb(i,2)*ibn22);
        j = j + 1;
    end
end
argumentos.AmiTotal = AmiTotal;

%%
%Difeomorfismo:
%Obs.: X = sigDif*Z
sigDif01 = [beta01; alfa01; alfa11];
[ Ami1 ] = fAG( sigDif01, AmiTotal, Ab );
sigDif1 = [sigDif01; Ami1];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FORMA NORMAL SEM LQR (saídas ref. rotor principal)
%Obs.: 
%(1) Z é o vetor original de estado
%(2) Saídas: tetaXpNM1 e Nr = kNM1*ubNM1

%Saídas: y1 = tetaXpNM1 e y2 = kNM1*ubNM1 (ATENÇÃO!)
SigmaR = argumentos.SigmaR;
PkNg = SigmaR*Pk;
Ng = argumentos.Ng;
kNM1 = argumentos.kNM1;
C1NM1 = [zeros(1,sPk(1,2)),PkNg(Ng,:)];
C2NM1 = [PkNg(Ng - 5,:),zeros(1,sPk(1,2))];

%%
%Dinâmica externa:
%Obs.: yTil = aTils(X) + bTilb*Us
beta0NM1 = C1NM1;
beta1NM1 = beta0NM1*Ab;
beta2NM1 = beta1NM1*Ab;
beta3NM1 = beta2NM1*Ab;
beta4NM1 = beta3NM1*Ab;
beta5NM1 = beta4NM1*Ab;
beta4NM1Bb = beta4NM1*Bb;
%
beta0NM1Bb = beta0NM1*Bb;
beta1NM1Bb = beta1NM1*Bb;
beta2NM1Bb = beta2NM1*Bb;
beta3NM1Bb = beta3NM1*Bb;
beta6NM1 = beta5NM1*Ab;
beta5NM1Bb = beta5NM1*Bb;
beta7NM1 = beta6NM1*Ab;
beta6NM1Bb = beta6NM1*Bb;
beta8NM1 = beta7NM1*Ab;
beta7NM1Bb = beta7NM1*Bb;
beta9NM1 = beta8NM1*Ab;
beta8NM1Bb = beta8NM1*Bb;
beta10NM1 = beta9NM1*Ab;
beta9NM1Bb = beta9NM1*Bb;
beta11NM1 = beta10NM1*Ab;
beta10NM1Bb = beta10NM1*Bb;
%}
alfa0NM1 = kNM1*C2NM1;
alfa1NM1 = alfa0NM1*Ab;
alfa2NM1 = alfa1NM1*Ab;
alfa1NM1Bb = alfa1NM1*Bb;
%
alfa0NM1Bb = alfa0NM1*Bb;
alfa3NM1 = alfa2NM1*Ab;
alfa2NM1Bb = alfa2NM1*Bb;
alfa4NM1 = alfa3NM1*Ab;
alfa3NM1Bb = alfa3NM1*Bb;
ad = [beta0NM1Bb, beta1NM1Bb, beta2NM1Bb, beta3NM1Bb, beta4NM1Bb, beta5NM1Bb, beta6NM1Bb, beta7NM1Bb, beta8NM1Bb, beta9NM1Bb, beta10NM1Bb, alfa0NM1Bb, alfa1NM1Bb, alfa2NM1Bb, alfa3NM1Bb];
%}

%%
%Dinâmica interna:
%Obs.:
%(1) mi = Ami*Z
%(2) mip = GamaMi*mi + SigmaMi*yTil + Ami*Bb*Uc

%%
%Difeomorfismo:
%Obs.: X = sigDif*Z
sigDif0NM1 = [beta0NM1; beta1NM1; beta2NM1; beta3NM1; beta4NM1; alfa0NM1; alfa1NM1];
%
disp('Ami: inicio de verificação')
[ AmiNM1 ] = fAG( sigDif0NM1, AmiTotal, Ab );
disp('Ami: fim de verificação')
%}
%Obs.:
%{
%(0) Tempo de verificação (no PC do Lucas M.):
%    (a) com impressão de "i": 6min50s
%    (b) sem impressão de "i": 3min49s
%(1) Nesse caso, são mais de 1.3*10^6 possibilidades a serem testadas
%(2) Esse cálculo já foi feito e consta no anexo C da tese escrita
%(3) Resultado mostrou que no mínimo, Fmi = 1 (uma instabilidade)
%(4) Esse fato motivou a utilização do LQR
%(5) Apenas para constar, segue um resultado com a supressão da linha lmi:
%}
%{
lmiNM1 = 13; %Linha retirada de AmiTotal, com "melhor" resultado geral
AmiNM1 = [AmiTotal(1:lmiNM1-1,:); AmiTotal(lmiNM1+1:nAb-2,:)];
%}
sigDifNM1 = [sigDif0NM1; AmiNM1];
%
condSigDifNM1 = cond(sigDifNM1);
sAmiNM1 = size(AmiNM1);
lAmiNM1 = sAmiNM1(1,1);
sSigDif0NM1 = size(sigDif0NM1);
lSigDif0NM1 = sSigDif0NM1(1,1);
Ad = IdAb(:,lSigDif0NM1+1:lSigDif0NM1+lAmiNM1);
GamaMi = AmiNM1*Ab*(sigDifNM1\Ad);
eGamaMi = eig(GamaMi);
ReaLeGamaMi = real(eGamaMi);
idRealeGamaMi = double(ReaLeGamaMi>0);
nRealeAbs = sum(idRealeGamaMi);
maxRealeAbs = max(abs(ReaLeGamaMi));
minRealeAbs = min(abs(ReaLeGamaMi));
ad = 1;
%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Regulador com quebra de dinâmicas de vibração:
Qs01 = eye(4);
Rs01 = diag([10^-3, 10^-3]);
[K01,Ss01,eS01] = lqr(A01,B01,Qs01,Rs01,Ns);
info.Ss01 = Ss01;
info.eS01 = eS01;
KlqrQv = [K01 + K21*(A221\(A211 - B21*K01)), K21];
%Obs.: Uc = -KlqrQv*Zbb
%Autalores do sistema:
%(1) Planta (metamodelo) e Controlador (metamodelo):
eLQRqv = eig(AbbLQR - BbbLQR*KlqrQv);
%{
ReaLeLQR = real(eLQR);
idRealeLQR = double(ReaLeLQR>0);
nRealeLQR = sum(idRealeLQR);
maxRealeLQR = max(abs(ReaLeLQR));
minRealeLQR = min(abs(ReaLeLQR));
log10dOG = log10(maxRealeLQR) - log10(minRealeLQR);
%ok!
%Obs.:
%(1) Todos os autovalores com parte real negativa
%(2) Módulo da menor parte real ~10^-4
%09/06/2022
%}
info.eLQRqv = eLQRqv;
%(2) Planta (modelo sem aprox de CR) e Controlador (metamodelo):
ePLQRqv = eig(Ab - Bb*KlqrQv*Perlqv);
%{
ReaLePLQR = real(ePLQRqv);
idRealePLQR = double(ReaLePLQR>0);
nRealePLQR = sum(idRealePLQR);
maxRealePLQR = max(abs(ReaLePLQR));
minRealePLQR = min(abs(ReaLePLQR));
log10dOG = log10(maxRealePLQR) - log10(minRealePLQR);
%ok!
%Obs.:
%(1) Todos os autovalores com parte real negativa
%(2) Módulo da menor parte real: ~1.37*10^-4
%26/06/2022
%}
info.ePLQRqv = ePLQRqv;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FORMA NORMAL COM LQR (saídas ref. rotor principal)
%Obs.: 
%(1) Z é o vetor original de estado
%(2) Saídas: tetaXpNM1 e Nr = kNM1*ubNM1
%(3) Tentativa de absorção dos modos reversos
eps0 = 10^-9;
Qrev = eps0*eye(2*jPu2);
Qrev(jPv,jPv) = 1;
Qrev(jPw,jPw) = 1;
Qrev(jPu2 + jPv,jPu2 + jPv) = eps0;
Qrev(jPu2 + jPw,jPu2 + jPw) = eps0;
Rrev = diag([10^-3, 10^-3]);
[Krev,Ss01,eS01] = lqr(Ab,Bb,Qrev,Rrev);
info.Ss01 = Ss01;
info.eS01 = eS01;
AbRev = Ab - Bb*Krev;

%Dinâmica externa:
%Obs.: yTil = aTils(X) + bTilb*Us
beta0NM1 = C1NM1;
beta1NM1 = beta0NM1*AbRev;
beta2NM1 = beta1NM1*AbRev;
beta3NM1 = beta2NM1*AbRev;
beta4NM1 = beta3NM1*AbRev;
beta5NM1 = beta4NM1*AbRev;
beta4NM1Bb = beta4NM1*Bb;
%{
beta0NM1Bb = beta0NM1*Bb;
beta1NM1Bb = beta1NM1*Bb;
beta2NM1Bb = beta2NM1*Bb;
beta3NM1Bb = beta3NM1*Bb;
beta6NM1 = beta5NM1*Ab;
beta5NM1Bb = beta5NM1*Bb;
beta7NM1 = beta6NM1*Ab;
beta6NM1Bb = beta6NM1*Bb;
beta8NM1 = beta7NM1*Ab;
beta7NM1Bb = beta7NM1*Bb;
beta9NM1 = beta8NM1*Ab;
beta8NM1Bb = beta8NM1*Bb;
beta10NM1 = beta9NM1*Ab;
beta9NM1Bb = beta9NM1*Bb;
beta11NM1 = beta10NM1*Ab;
beta10NM1Bb = beta10NM1*Bb;
%}
alfa0NM1 = kNM1*C2NM1;
alfa1NM1 = alfa0NM1*AbRev;
alfa2NM1 = alfa1NM1*AbRev;
alfa1NM1Bb = alfa1NM1*Bb;
%{
alfa0NM1Bb = alfa0NM1*Bb;
alfa3NM1 = alfa2NM1*Ab;
alfa2NM1Bb = alfa2NM1*Bb;
alfa4NM1 = alfa3NM1*Ab;
alfa3NM1Bb = alfa3NM1*Bb;
%ad = [beta0NM1Bb, beta1NM1Bb, beta2NM1Bb, beta3NM1Bb, beta4NM1Bb, beta5NM1Bb, beta6NM1Bb, beta7NM1Bb, beta8NM1Bb, beta9NM1Bb, beta10NM1Bb, alfa0NM1Bb, alfa1NM1Bb, alfa2NM1Bb, alfa3NM1Bb];
ad = [beta0NM1Bb, beta1NM1Bb, beta2NM1Bb, beta3NM1Bb, beta4NM1Bb, alfa0NM1Bb, alfa1NM1Bb, alfa2NM1Bb];
%}

%%
%Difeomorfismo:
%Obs.: X = sigDif*Z
sigDif0NM1 = [beta0NM1; beta1NM1; beta2NM1; beta3NM1; beta4NM1; alfa0NM1; alfa1NM1];
%{
disp('Ami: inicio de verificação')
[ AmiNM1 ] = fAG( sigDif0NM1, AmiTotal, Ab );
disp('Ami: fim de verificação')
%}
sigDifNM1 = [sigDif0NM1; AmiNM1];
%{
condSigDifNM1 = cond(sigDifNM1);
sAmiNM1 = size(AmiNM1);
lAmiNM1 = sAmiNM1(1,1);
sSigDif0NM1 = size(sigDif0NM1);
lSigDif0NM1 = sSigDif0NM1(1,1);
Ad = IdAb(:,lSigDif0NM1+1:lSigDif0NM1+lAmiNM1);
GamaMi = AmiNM1*Ab*(sigDifNM1\Ad);
eGamaMi = eig(GamaMi);
ReaLeGamaMi = real(eGamaMi);
idRealeGamaMi = double(ReaLeGamaMi>0);
nRealeAbs = sum(idRealeGamaMi);
maxRealeAbs = max(abs(ReaLeGamaMi));
minRealeAbs = min(abs(ReaLeGamaMi));
ad = 1;
%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Armazenamento:
%(1) Ganho do regulador:
ABvEps(1,1) = {Klqr};
%(2) Matrizes da dinâmica, com absorvedor de vibração:
ABvEps(1,2) = {Abs}; ABvEps(1,3) = {Bbs};
%(3) Matrizes da dinâmica externa (com LQR):
ABvEps(1,4) = {beta0s}; ABvEps(1,5) = {beta1s}; 
ABvEps(1,6) = {beta0Bbs};
ABvEps(1,7) = {alfa0s}; ABvEps(1,8) = {alfa1s}; ABvEps(1,9) = {alfa2s}; 
ABvEps(1,10) = {alfa1Bbs};
%(4) Matrizes dinâmica interna e difeomorfismo (com LQR):
ABvEps(1,11) = {AmiS}; ABvEps(1,12) = {sigDifs};
%(5) Matrizes da dinâmica, com absorvedor de vibração:
ABvEps(1,13) = {Abbs}; ABvEps(1,14) = {Bbbs};
%(6) Matrizes da dinâmica externa (com LQR):
ABvEps(1,15) = {beta0s2}; ABvEps(1,16) = {beta1s2}; 
ABvEps(1,17) = {beta0Bbs2};
ABvEps(1,18) = {alfa0s2}; ABvEps(1,19) = {alfa1s2}; ABvEps(1,20) = {alfa2s2}; 
ABvEps(1,21) = {alfa1Bbs2};
%(7) Matrizes dinâmica interna e difeomorfismo (com LQR):
ABvEps(1,22) = {AmiS2}; ABvEps(1,23) = {sigDifs2}; 
%(8) Matrizes da dinâmica externa (sem LQR):
ABvEps(1,24) = {beta01}; ABvEps(1,25) = {beta11};
ABvEps(1,26) = {beta01Bb};
ABvEps(1,27) = {alfa01}; ABvEps(1,28) = {alfa11}; ABvEps(1,29) = {alfa21}; 
ABvEps(1,30) = {alfa11Bb};
%(9) Matrizes dinâmica interna e difeomorfismo (sem LQR):
ABvEps(1,31) = {Ami1}; ABvEps(1,32) = {sigDif1};
%(10) Matrizes da dinâmica externa (sem LQR):
ABvEps(1,33) = {beta0NM1}; ABvEps(1,34) = {beta1NM1}; 
ABvEps(1,35) = {beta2NM1}; ABvEps(1,36) = {beta3NM1}; 
ABvEps(1,37) = {beta4NM1}; ABvEps(1,28) = {beta5NM1};
ABvEps(1,39) = {beta4NM1Bb};
ABvEps(1,40) = {alfa0NM1}; ABvEps(1,41) = {alfa1NM1}; 
ABvEps(1,42) = {alfa2NM1};
ABvEps(1,43) = {alfa1NM1Bb};
%(11) Matrizes dinâmica interna e difeomorfismo (com LQR):
ABvEps(1,44) = {AmiNM1}; ABvEps(1,45) = {sigDifNM1};
%(12) Regulador, com quebra da dinâmica de vibração: fzS1
ABvEps(1,46) = {KlqrQv}; ABvEps(1,47) = {fzS1}; ABvEps(1,48) = {fzS2}; 
end