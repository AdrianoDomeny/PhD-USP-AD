function [ custo, rSig, Fmi, cSig, RealMax ] = fAGverify( v, sigDif0, AmiTotal, Ab )
%Função CUSTO para otimização por algorítmos genéticos:
%   v: vetor linha a ser avaliado
Ami = AmiTotal(v,:);
sigDif = [sigDif0; Ami];
sAb = size(Ab);
sSig0 = size(sigDif0);
rSig = sAb(1,1) - rank(sigDif);
IdAb = eye(sAb(1,1));
Ad = IdAb(:,(sSig0(1,1)+1):sAb(1,1)); 
if rSig < 1
    GamaMi = Ami*Ab*(sigDif\Ad);
    eGamaMi = eig(GamaMi);
    [maxeG,jM] = max(real(eGamaMi));
    if maxeG>0 && maxeG<10^-17
        eGamaMi(jM,1) = 0;
    end
    ReaLeGamaMi = real(eGamaMi);
    RealMax = max(ReaLeGamaMi);
    idReaLeGamaMi = double(ReaLeGamaMi>0);
    nReaLeGamaMi = sum(idReaLeGamaMi);
    Fmi = nReaLeGamaMi;
    cSig = abs(log10(abs(cond(sigDif))));
else
    RealMax = 1000;
    Fmi = sAb(1,1) - sSig0(1,1);
    cSig = 30;
end
k1 = 10^2;
k2 = 10^2;
custo = k1*rSig + k2*Fmi + cSig;
end