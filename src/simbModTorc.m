clear
clc

%%
%Sistema linear: M*uppTil + C*upTil + K*uTil = B*T
syms I1 I2 I3 K1 K2 alfa beta c1 c3
beta = 0;
M = diag([I3, I2, I1]);
K = [K1, -K1, 0; -K1, K1+K2, -K2; 0, -K2, K2];
Cestr = beta*K;
Conc = diag([c1, 0, c3]);
C = Cestr + Conc;
Id = eye(3);
B = Id(:,1);

%%
%Equações de estado: Zp = A*Z + Bb*T
A = [zeros(3), Id; -M\K, -M\C];
Bb = [zeros(3,1); M\B];

%%
%Matriz total para escolha da dinâmica interna: AmiTotal; AmiTotal*Bb = 0
n = length(M);
[bMax,jbMax] = max(abs(B));
jMax = 3 + jbMax;
bjMax = bMax;
ibjMax = bjMax^-1;
AmiTotal = zeros(2*n-1,2*n);
IdA = eye(2*n);
j = 1;
for i=1:2*n
    if i ~= jMax
        AmiTotal(j,:) = IdA(i,:) - IdA(jMax,:)*Bb(i,1)*ibjMax;
        j = j + 1;
    end
end
AmiTotalBb = AmiTotal*Bb;

%%
%Saída: y = tetap1
C3 = [zeros(1,3), Id(1,:)];
beta03 = C3;
beta13 = beta03*A;
beta03Bb = beta03*Bb;
sigDif03 = beta03;
sigDif03Bb = sigDif03*Bb;
%Ordem relativa:
sSigDif03 = size(sigDif03);
r03 = sSigDif03(1,1);
%Ordem da dinâmica interna:
nMi03 = 2*n - r03;
%Dinâmica interna: miTil = Ami*Z

%%
%Saída: y = tetap2
C2 = [zeros(1,3), Id(2,:)];
beta02 = C2;
beta12 = beta02*A;
beta02Bb = beta02*Bb;
beta22 = beta12*A;
beta12Bb = beta12*Bb;
%Obs.: 
%(1) beta12Bb é nulo apenas na presença de amortecimento estrutural
beta32 = beta22*A;
beta22Bb = beta22*Bb;
sigDif02 = [beta02; beta12; beta22];
sigDif02Bb = sigDif02*Bb;
%Ordem relativa:
sSigDif02 = size(sigDif02);
r02 = sSigDif02(1,1);
%Ordem da dinâmica interna:
nMi02 = 2*n - r02;

%%
%Saída: y = tetap3
C1 = [zeros(1,3), Id(3,:)];
beta01 = C1;
beta11 = beta01*A;
beta01Bb = beta01*Bb;
beta21 = beta11*A;
beta11Bb = beta11*Bb;
beta31 = beta21*A;
beta21Bb = beta21*Bb;
beta41 = beta31*A;
beta31Bb = beta31*Bb;
beta51 = beta41*A;
beta41Bb = beta41*Bb;
sigDif01 = [beta01; beta11; beta21; beta31; beta41];
sigDif01Bb = sigDif01*Bb;
%Ordem relativa:
sSigDif01 = size(sigDif01);
r01 = sSigDif01(1,1);
%Ordem da dinâmica interna:
nMi01 = 2*n - r01;




