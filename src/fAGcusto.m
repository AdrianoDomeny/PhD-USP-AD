function [ custo ] = fAGcusto( v, argumentos )
%Função CUSTO para otimização por algorítmos genéticos:
%   v: vetor linha a ser avaliado
sigDif0 = argumentos.sigDif0;
AmiTotal = argumentos.AmiTotal;
Ab = argumentos.Ab;
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
    ReaLeGamaMi = real(eGamaMi);
    idReaLeGamaMi = double(ReaLeGamaMi>0);
    nReaLeGamaMi = sum(idReaLeGamaMi);
    Fmi = nReaLeGamaMi;
    cSig = abs(log10(abs(cond(sigDif))));
else
    Fmi = sAb(1,1) - sSig0(1,1);
    cSig = 30;
end
k1 = 10^2;
k2 = 10^2;
custo = k1*rSig + k2*Fmi + cSig;
end

