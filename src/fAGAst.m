function [ AmiOp, xOp, x , CUSTO, Rsig, FMI, CSIG, REALMax ] = fAGAst( sigDif0, AmiAstTotal, Ab, M )
%Otimização por Algorítmos Genéticos:
%   D: número de linhas em cada indivíduo
%   NP: número de indivíduos numa população
%   CR: probabilidade de cruzamento, CR pertence a [0,1]
%   F: fator de mutação, F pertence a {0,1,2}
%   LB e UB: delimitadores dos indivíduos
sSig = size(sigDif0);
sAmi = size(AmiAstTotal);
LB = 1;
UB = sAmi(1,1);
D = sAmi(1,2) - sSig(1,1);
v = LB:UB;
x = nchoosek(v,D);
sC = size(x);
%fAGc = @(v)(fAGcusto(v, argumentos));
CUSTO = zeros(sC(1,1),1);
Rsig = zeros(sC(1,1),1);
FMI = zeros(sC(1,1),1);
CSIG = zeros(sC(1,1),1);
REALMax = zeros(sC(1,1),1);
minCUSTO = 10^200;
for i = 1:sC(1,1)
    [ custo, rSig, Fmi, cSig, RealMax ] = fAGastVer( x(i,:), sigDif0, AmiAstTotal, Ab, M );
    CUSTO(i,1) = custo;
    Rsig(i,1) = rSig;
    FMI(i,1) = Fmi;
    nFmi = 1;
    %{
    if Fmi<nFmi+1
        iOp = i;
        minCUSTO = custo;
    end
    %}
    %
    if Fmi<nFmi+1 && custo<=minCUSTO
        iOp = i;
        minCUSTO = custo;
    end
    %}
    CSIG(i,1) = cSig;
    REALMax(i,1) = RealMax;
    i
end
xOp = x(iOp,:);
AmiOp = AmiAstTotal(xOp,:);
end