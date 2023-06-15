function [ ABvEps, info ] = fModOrdGrand( argumentos )
%UNTITLED Summary of this function goes here
%   Ordem de grandeza das frequências(^2)
Ab = argumentos.Ab;
Bb = argumentos.Bb;
%Ng = argumentos.Ng;
%rEstr = argumentos.rEstr;
%Krk = argumentos.Krk;
Crk = argumentos.Crk;
%Brk = argumentos.Brk;
sigTilk = argumentos.sigTilk;
Pk = argumentos.Pk;
%Rs = argumentos.Rs;
Ns = argumentos.Ns;
%Rf = argumentos.Rf;
Nf = argumentos.Nf;
jPtX0 = argumentos.jPtX0;
jPu0 = argumentos.jPu0;
jPtX1 = argumentos.jPtX1;
jPtX2 = argumentos.jPtX2;
jPv = argumentos.jPv;
jPw = argumentos.jPw;
jPu1 = argumentos.jPu1;
jPu2 = argumentos.jPu2;
epsZero = argumentos.epsZero;
nSigTilk = length(sigTilk);
idSig = double(abs(sigTilk) < epsZero*10^4);
ndSig0 = sum(idSig);

%%
nAb = length(Ab);

%{
dKbrk = sigTilk(ndSig0+1:nSigTilk,1);
OrdGrand = floor(log10(dKbrk)).*double(log10(dKbrk) - floor(log10(dKbrk)) < 0.5) +...
    ceil(log10(dKbrk)).*double(log10(dKbrk) - floor(log10(dKbrk)) >= 0.5);
OGmax = max(OrdGrand);
Nordg = zeros(OGmax+1,1);
sep = 6;
switch sep
    case 1
        %Primeira separação das dinâmicas:
        epS = 1;
        A11 = Ab(1:nSigTilk+ndSig0,1:nSigTilk+ndSig0);
        A12 = Ab(1:nSigTilk+ndSig0,nSigTilk+ndSig0+1:2*nSigTilk);
        A21 = epS*Ab(nSigTilk+ndSig0+1:2*nSigTilk,1:nSigTilk+ndSig0);
        A22 = epS*Ab(nSigTilk+ndSig0+1:2*nSigTilk,nSigTilk+ndSig0+1:2*nSigTilk);
        B1 = Bb(1:nSigTilk+ndSig0,:);
        B2 = epS*Bb(nSigTilk+ndSig0+1:2*nSigTilk,:);
        A0 = A11 - A12*(A22\A21);
        B0 = B1 - A12*(A22\B2);
        varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
        %Ganhos parciais:
        %(1) Dinâmica lenta:
        sA0 = size(A0);
        Qs = eye(sA0);
        [Ks,Ss,eS] = lqr(A0,B0,Qs,Rs,Ns);
        info.Ss = Ss;
        info.eS = eS;
        %Armazenamento:
        ABvEps(1,1) = {A0}; ABvEps(1,2) = {B0}; 
        ABvEps(1,3) = {A21}; ABvEps(1,4) = {A22}; ABvEps(1,5) = {B2};
        ABvEps(1,6) = {varDif}; ABvEps(1,7) = {epS};
        ABvEps(1,8) = {Ks}; 
        %Demais separações:
        OG = 0;
        j = 2;
        nSig = 0;
        Aaux = A22;
        Baux = B2;
        for i = 1:OGmax
            idOG = double((OG-1 < OrdGrand) & (OrdGrand < OG+1));
            nOG = sum(idOG);
            Nordg(i,1) = nOG;
            if nOG > 0
                nAaux = length(Aaux);
                epS = sigTilk(ndSig0+nSig+1,1)/sigTilk(ndSig0+nSig+nOG+1,1);
                A11 = Aaux(1:nOG,1:nOG);
                A12 = Aaux(1:nOG,nOG+1:nAaux);
                A21 = epS*Aaux(nOG+1:nAaux,1:nOG);
                A22 = epS*Aaux(nOG+1:nAaux,nOG+1:nAaux);
                B1 = Baux(1:nOG,:);
                B2 = epS*Baux(nOG+1:nAaux,:);
                A0 = A11 - A12*(A22\A21);
                B0 = B1 - A12*(A22\B2);
                varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
                %Ganhos parciais:
                %(1) Dinâmica lenta:
                sA0 = size(A0);
                Qs = eye(sA0);
                [Ks] = lqr(A0,B0,Qs,Rs,Ns);
                %Armazenamento:
                ABvEps(j,1) = {A0}; ABvEps(j,2) = {B0}; 
                ABvEps(j,3) = {A21}; ABvEps(j,4) = {A22}; ABvEps(j,5) = {B2};
                ABvEps(j,6) = {varDif}; ABvEps(j,7) = {epS};
                ABvEps(j,8) = {Ks}; 
                %Atualização:
                Aaux = A22;
                Baux = B2;
                j = j + 1;
                nSig = nSig + nOG;
            end
            OG = OG + 1;
        end
        %(2) Dinâmica rápida:
        sA22 = size(A22);
        Qf = eye(sA22);
        [Kf,Sf,eF] = lqr(A22,B2,Qf,Rf,Nf);
        info.Sf = Sf;
        info.eF = eF;
        %Ganho final de controle:
        for i = j-1:-1:1
            A21 = ABvEps{i,3};
            A22 = ABvEps{i,4};
            B2 = ABvEps{i,5};
            Ks = ABvEps{i,8};
            [ Kf ] = fModKfn( Ks, Kf, A22, A21, B2  );
            %ad = 1;
        end
        Klqr = Kf;
        eLQR = eig(Ab + Bb*Klqr);
        %{
        ReaLeLQR = real(eLQR);
        idRealeLQR = double(ReaLeLQR>0);
        nRealeLQR = sum(idRealeLQR);
        %}
        info.Klqr = Klqr;
        info.eLQR = eLQR;
    case 2
        %Primeira separação das dinâmicas:
        epS = 1;
        A11 = Ab(1:nSigTilk+ndSig0,1:nSigTilk+ndSig0);
        A12 = Ab(1:nSigTilk+ndSig0,nSigTilk+ndSig0+1:2*nSigTilk);
        A21 = epS*Ab(nSigTilk+ndSig0+1:2*nSigTilk,1:nSigTilk+ndSig0);
        A22 = epS*Ab(nSigTilk+ndSig0+1:2*nSigTilk,nSigTilk+ndSig0+1:2*nSigTilk);
        B1 = Bb(1:nSigTilk+ndSig0,:);
        B2 = epS*Bb(nSigTilk+ndSig0+1:2*nSigTilk,:);
        A0 = A11 - A12*(A22\A21);
        B0 = B1 - A12*(A22\B2);
        varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
        %Ganhos parciais:
        %(1) Dinâmica lenta:
        sA0 = size(A0);
        Qs = eye(sA0);
        [Ks,Ss,eS] = lqr(A0,B0,Qs,Rs,Ns);
        info.Ss = Ss;
        info.eS = eS;
        %Armazenamento:
        ABvEps(1,1) = {A0}; ABvEps(1,2) = {B0}; 
        ABvEps(1,3) = {A21}; ABvEps(1,4) = {A22}; ABvEps(1,5) = {B2};
        ABvEps(1,6) = {varDif}; ABvEps(1,7) = {epS};
        ABvEps(1,8) = {Ks}; 
        %Demais separações:
        OG = 0;
        j = 2;
        nSig = 0;
        Aaux = A22;
        Baux = B2;
        for i = 1:1
            idOG = double((OG-1 < OrdGrand) & (OrdGrand < OG+5));
            nOG = sum(idOG);
            Nordg(i,1) = nOG;
            if nOG > 0
                nAaux = length(Aaux);
                epS = sigTilk(ndSig0+nSig+1,1)/sigTilk(ndSig0+nSig+nOG+1,1);
                A11 = Aaux(1:nOG,1:nOG);
                A12 = Aaux(1:nOG,nOG+1:nAaux);
                A21 = epS*Aaux(nOG+1:nAaux,1:nOG);
                A22 = epS*Aaux(nOG+1:nAaux,nOG+1:nAaux);
                B1 = Baux(1:nOG,:);
                B2 = epS*Baux(nOG+1:nAaux,:);
                A0 = A11 - A12*(A22\A21);
                B0 = B1 - A12*(A22\B2);
                varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
                %Ganhos parciais:
                %(1) Dinâmica lenta:
                sA0 = size(A0);
                Qs = eye(sA0);
                [Ks] = lqr(A0,B0,Qs,Rs,Ns);
                %Armazenamento:
                ABvEps(j,1) = {A0}; ABvEps(j,2) = {B0}; 
                ABvEps(j,3) = {A21}; ABvEps(j,4) = {A22}; ABvEps(j,5) = {B2};
                ABvEps(j,6) = {varDif}; ABvEps(j,7) = {epS};
                ABvEps(j,8) = {Ks}; 
                %Atualização:
                Aaux = A22;
                Baux = B2;
                j = j + 1;
                nSig = nSig + nOG;
            end
            OG = OG + 5;
        end
        %(2) Dinâmica rápida:
        sA22 = size(A22);
        Qf = eye(sA22);
        [Kf,Sf,eF] = lqr(A22,B2,Qf,Rf,Nf);
        info.Sf = Sf;
        info.eF = eF;
        %Ganho final de controle:
        for i = 2:-1:1
            A21 = ABvEps{i,3};
            A22 = ABvEps{i,4};
            B2 = ABvEps{i,5};
            Ks = ABvEps{i,8};
            [ Kf ] = fModKfn( Ks, Kf, A22, A21, B2  );
            %ad = 1;
        end
        Klqr = Kf;
        eLQR = eig(Ab - Bb*Klqr);
        %{
        ReaLeLQR = real(eLQR);
        idRealeLQR = double(ReaLeLQR>0);
        nRealeLQR = sum(idRealeLQR);
        %}
        info.Klqr = Klqr;
        info.eLQR = eLQR;
    case 3
        %%
        jPtX0 = argumentos.jPtX0;
        jPu0 = argumentos.jPu0;
        jPtX1 = argumentos.jPtX1;
        jPu1 = argumentos.jPu1;
        jPtX2 = argumentos.jPtX2;
        jPu2 = argumentos.jPu2;
        jPv = argumentos.jPv;
        jPw = argumentos.jPw;
        %Aproximação dos modos de corpo rígido:
        eps0 = 10^-10;
        %eps0 = 10^-6;
        sigTilk(jPtX0,1) = eps0;
        sigTilk(jPu0,1) = eps0;
        %Redefinição de Krk, pela aproximação das freq. de CR:
        Krk = diag(sigTilk.');
        %Separação entre modos de vibração e de corpo rígido:
        sKrk = size(Krk);
        IdKrk = eye(sKrk);
        Ark0 = IdKrk(:,1:2);
        Arkv = IdKrk(:,3:sKrk(1,2));
        Ark0CrkArk0 = Ark0.'*(Crk*Ark0);
        Ark0CrkArkv = Ark0.'*(Crk*Arkv);
        ArkvCrkArk0 = Arkv.'*(Crk*Ark0);
        ArkvCrkArkv = Arkv.'*(Crk*Arkv);
        Ark0KrkArk0 = Ark0.'*(Krk*Ark0);
        Ark0KrkArkv = Ark0.'*(Krk*Arkv);
        ArkvKrkArk0 = Arkv.'*(Krk*Ark0);
        ArkvKrkArkv = Arkv.'*(Krk*Arkv);
        Ark0Brk = Ark0.'*Brk;
        ArkvBrk = Arkv.'*Brk;
        %{
        Ark0CrkArkv = Ark0CrkArkv - Ark0CrkArkv.*double(abs(Ark0CrkArkv) < 10^2*epsZero);
        ArkvCrkArk0 = ArkvCrkArk0 - ArkvCrkArk0.*double(abs(ArkvCrkArk0) < 10^2*epsZero);
        ArkvCrkArkv = ArkvCrkArkv - ArkvCrkArkv.*double(abs(ArkvCrkArkv) < 10^2*epsZero);
        ArkvBrk = ArkvBrk - ArkvBrk.*double(abs(ArkvBrk) < 10^2*epsZero);
        %}
        Crk0 = Ark0CrkArk0;
        Crk0v = Ark0CrkArkv;
        Crkv0 = ArkvCrkArk0;
        Crkv = ArkvCrkArkv;
        Kbrk0 = Ark0KrkArk0;
        Kbrk0v = Ark0KrkArkv;
        Kbrkv0 = ArkvKrkArk0;
        Kbrkv = ArkvKrkArkv;
        Brk0 = Ark0Brk;
        Brkv = ArkvBrk;
        % Rearranjo da equação original da dinâmica em espaço de estados:
        sPk = size(Pk);
        Ab0 = [zeros(ndSig0), eye(ndSig0);
            -Kbrk0, -Crk0];
        Ab0v = [zeros(ndSig0,sPk(1,2)-ndSig0), zeros(ndSig0,sPk(1,2)-ndSig0);
            -Kbrk0v, -Crk0v];
        Abv0 = [zeros(sPk(1,2)-ndSig0,ndSig0), zeros(sPk(1,2)-ndSig0,ndSig0);
            -Kbrkv0, -Crkv0];
        Abv = [zeros(sPk(1,2)-ndSig0), eye(sPk(1,2)-ndSig0);
            -Kbrkv, -Crkv];
        Bbrk0 = [zeros(ndSig0,2); 
            Brk0];
        Bbrkv = [zeros(sPk(1,2)-ndSig0,2); 
            Brkv];
        Abb = [Ab0, Ab0v;
            Abv0, Abv];
        Bbb = [Bbrk0;
            Bbrkv];
        %Separação das dinâmicas:
        epS = sigTilk(1,1)/sigTilk(ndSig0+1,1);
        A11 = Abb(1:2*ndSig0,1:2*ndSig0);
        A12 = Abb(1:2*ndSig0,2*ndSig0+1:2*sPk(1,2));
        A21 = epS*Abb(2*ndSig0+1:2*sPk(1,2),1:2*ndSig0);
        A22 = epS*Abb(2*ndSig0+1:2*sPk(1,2),2*ndSig0+1:2*sPk(1,2));
        B1 = Bbb(1:2*ndSig0,:);
        B2 = epS*Bbb(2*ndSig0+1:2*sPk(1,2),:);
        A0 = A11 - A12*(A22\A21);
        B0 = B1 - A12*(A22\B2);
        varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
        %Ganhos parciais:
        %(1) Dinâmica lenta:
        sA0 = size(A0);
        Qs = eye(sA0);
        Rs = diag([10^-3, 10^-3]); %Testar diferentes valores
        [K0,Ss,eS] = lqr(A0,B0,Qs,Rs,Ns);
        info.Ss = Ss;
        info.eS = eS;
        %(2) Dinâmica rápida:
        sA22 = size(A22);
        Qf = eye(sA22);
        %Pesos em uTilrk0:
        Qf(jPtX1-ndSig0,jPtX1-ndSig0) = 10^3;
        Qf(jPtX2-ndSig0,jPtX2-ndSig0) = 10^3;
        Qf(jPv-ndSig0,jPv-ndSig0) = 10^3;
        Qf(jPw-ndSig0,jPw-ndSig0) = 10^3;
        Qf(jPu1-ndSig0,jPu1-ndSig0) = 10^3;
        Qf(jPu2-ndSig0,jPu2-ndSig0) = 10^3;
        %Pesos em vTilrk0:
        nrkv = sPk(1,2) - ndSig0;
        Qf(nrkv+jPtX1-ndSig0,nrkv+jPtX1-ndSig0) = 10^3;
        Qf(nrkv+jPtX2-ndSig0,nrkv+jPtX2-ndSig0) = 10^3;
        Qf(nrkv+jPv-ndSig0,nrkv+jPv-ndSig0) = 10^3;
        Qf(nrkv+jPw-ndSig0,nrkv+jPw-ndSig0) = 10^3;
        Qf(nrkv+jPu1-ndSig0,nrkv+jPu1-ndSig0) = 10^3;
        Qf(nrkv+jPu2-ndSig0,nrkv+jPu2-ndSig0) = 10^3;
        Rf = diag([10^-3, 10^-3]); %Testar diferentes valores
        [K2,Sf,eF] = lqr(A22,B2,Qf,Rf,Nf);
        info.Sf = Sf;
        info.eF = eF;
        
        %Ganho total:
        K1 = K0 + K2*(A22\(A21 - B2*K0));
        Klqr = [K1, K2];
        eLQR = eig(Abb - Bbb*Klqr);
        %{
        ReaLeLQR = real(eLQR);
        idRealeLQR = double(ReaLeLQR>0);
        nRealeLQR = sum(idRealeLQR);
        maxRealeLQR = max(abs(ReaLeLQR));
        minRealeLQR = min(abs(ReaLeLQR));
        log10dOG = log10(maxRealeLQR) - log10(minRealeLQR);
        %ok!
        %}
        info.eLQR = eLQR;
        %Resultado: 
        %(1) 2x Avl com parte real positiva no regulador
        
        %Novo sistema, com aplicação de regulador sobre a dinâmica rápida:
        %Obs.: Zp = Abs*Z + Bbs*Us
        %Redefinição de Ab, pela aproximação das freq. de CR:
        Idz = eye(2*(nSigTilk-ndSig0));
        Abs = Abb - Bbb*K2*[A22\A21, Idz];
        IdnSig0 = eye(ndSig0);
        Bbs = Bbb*(IdnSig0 - K2*(A22\B2));
        %Teste:
        eAbs = eig(Abs);
        info.eAbs = eAbs;
        %{
        ReaLeAbs = real(eAbs);
        idRealeAbs = double(ReaLeAbs>0);
        nRealeAbs = sum(idRealeAbs);
        maxRealeAbs = max(abs(ReaLeAbs));
        minRealeAbs = min(abs(ReaLeAbs));
        log10dOG = log10(maxRealeAbs) - log10(minRealeAbs);
        ad = 1;
        %ok!
        %}
        %Resultado: 
        %(1) Todos os avl de Abs com parte real negativa 
        
        %Matriz permuta de linhas Zb = Perl*Z e Z = Perl.'*Zb:
        sAbs = size(Abs);
        IdAbs = eye(sAbs);
        Perl = [IdAbs(1:ndSig0,:); 
            IdAbs(sPk(1,2)+1:sPk(1,2)+ndSig0,:); 
            IdAbs(ndSig0+1:sPk(1,2),:); 
            IdAbs(sPk(1,2)+ndSig0+1:2*sPk(1,2),:)];
        
        %Saídas:
        %{
        Ng = argumentos.Ng;
        rEstr = argumentos.rEstr;
        jtXNM1 = Ng-rEstr;
        C1 = [zeros(1,sPk(1,2)),Pk(jtXNM1,:)]*Perl.';
        juNM1 = Ng-rEstr-1;
        C2 = [Pk(juNM1,:),zeros(1,sPk(1,2))]*Perl.';
        kNM1 = argumentos.kNM1;
        %}
        %{
        C1 = [zeros(1,sPk(1,2)),Pk(2,:)]*Perl.';
        C2 = [Pk(1,:),zeros(1,sPk(1,2))]*Perl.';
        kNM1 = 1;
        %}
        
        %Análise linear do sistema "Planta + Regulador":
        testeSS = 1;
        if testeSS>0
            Cb = eye(sAbs)*Perl;
            sys = ss(Abs,Bbs,Cb,[]);
            t = 0:0.0001:50;
            [y,t] = impulse(sys,t);
            figure
            subplot(2,3,1)
            plot(t,y(:,jPtX1,1))
            grid
            title('modo tX1')
            xlabel('Tempo (s)')
            subplot(2,3,2)
            plot(t,y(:,jPtX2,1))
            grid
            title('modo tX2')
            xlabel('Tempo (s)')
            subplot(2,3,3)
            plot(t,y(:,jPv,1))
            grid
            title('modo v')
            xlabel('Tempo (s)')
            subplot(2,3,4)
            plot(t,y(:,jPw,1))
            grid
            title('modo w')
            xlabel('Tempo (s)')
            subplot(2,3,5)
            plot(t,y(:,jPu1,1))
            grid
            title('modo u1')
            xlabel('Tempo (s)')
            subplot(2,3,6)
            plot(t,y(:,jPu2,1))
            grid
            title('modo u2')
            xlabel('Tempo (s)')
            %{
            figure
            plot(t,y(:,sPk(1,2)-1,1))
            grid
            figure
            plot(t,y(:,sPk(1,2),1))
            grid
            %}
            %{
            figure
            plot(t,y(:,1,2))
            grid
            figure
            plot(t,y(:,2,2))
            grid
            %}
            %{
            figure
            plot(t,y(:,sPk(1,2)-1,2))
            grid
            figure
            plot(t,y(:,sPk(1,2),2))
            grid
            %}
        end
        
        beta0 = C1;
        beta1 = beta0*Abs;
        beta0Bbs = beta0*Bbs;
        beta2 = beta1*Abs;
        beta1Bbs = beta1*Bbs;
        beta3 = beta2*Abs;
        beta2Bbs = beta2*Bbs;
        beta4 = beta3*Abs;
        beta3Bbs = beta3*Bbs;
        beta5 = beta4*Abs;
        beta4Bbs = beta4*Bbs;
        beta6 = beta5*Abs;
        beta5Bbs = beta5*Bbs;
        beta7 = beta6*Abs;
        beta6Bbs = beta6*Bbs;
        beta8 = beta7*Abs;
        beta7Bbs = beta7*Bbs;
        beta9 = beta8*Abs;
        beta8Bbs = beta8*Bbs;
        beta10 = beta9*Abs;
        beta9Bbs = beta9*Bbs;
        beta11 = beta10*Abs;
        beta10Bbs = beta10*Bbs;
        
        alfa0 = kNM1*C2;
        alfa1 = alfa0*Abs;
        alfa0Bbs = alfa0*Bbs;
        alfa2 = alfa1*Abs;
        alfa1Bbs = alfa1*Bbs;
        alfa3 = alfa2*Abs;
        alfa2Bbs = alfa2*Bbs;
        alfa4 = alfa3*Abs;
        alfa3Bbs = alfa3*Bbs;
        ad = [beta0Bbs, beta1Bbs, beta2Bbs, beta3Bbs, beta4Bbs, beta5Bbs, beta6Bbs, beta7Bbs, beta8Bbs, beta9Bbs, beta10Bbs, alfa0Bbs, alfa1Bbs, alfa2Bbs, alfa3Bbs];
        
        %Armazenamento:
        ABvEps(1,1) = {A0}; ABvEps(1,2) = {B0}; 
        ABvEps(1,3) = {A21}; ABvEps(1,4) = {A22}; ABvEps(1,5) = {B2};
        ABvEps(1,6) = {varDif}; ABvEps(1,7) = {epS};
        ABvEps(1,8) = {K0}; ABvEps(1,9) = {K2}; ABvEps(1,10) = {Klqr};
    case 4
        sPk = size(Pk);
        %Aproximação dos modos de corpo rígido:
        %eps0 = 10^-10;
        eps0 = 10^-6;
        sigTilk(1:ndSig0,1) = eps0*ones(ndSig0,1);
        %Redefinição de Krk, pela aproximação das freq. de CR:
        Krk = diag(sigTilk.');
        %Redefinição de Ab, pela aproximação das freq. de CR:
        Ab = [zeros(sPk(1,2)), eye(sPk(1,2)); 
            -Krk, -Crk];
        %Primeira separação das dinâmicas, segundo os modos de vibração:
        nSigTilk = length(sigTilk);
        epS = sigTilk(1,1)/sigTilk(ndSig0+1,1);
        %epS = sigTilk(1,1)/sigTilk(nSigTilk,1);
        A11 = Ab(1:nSigTilk+ndSig0,1:nSigTilk+ndSig0);
        A12 = Ab(1:nSigTilk+ndSig0,nSigTilk+ndSig0+1:2*nSigTilk);
        A21 = epS*Ab(nSigTilk+ndSig0+1:2*nSigTilk,1:nSigTilk+ndSig0);
        A22 = epS*Ab(nSigTilk+ndSig0+1:2*nSigTilk,nSigTilk+ndSig0+1:2*nSigTilk);
        B1 = Bb(1:nSigTilk+ndSig0,:);
        B2 = epS*Bb(nSigTilk+ndSig0+1:2*nSigTilk,:);
        %B2 = B2 - B2.*double(abs(B2)<epsZero);
        A0 = A11 - A12*(A22\A21);
        B0 = B1 - A12*(A22\B2);
        varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
        
        %Ganhos parciais:
        %(1) Dinâmica lenta:
        sA0 = size(A0);
        Qs = eye(sA0);
        %{
        Qs = zeros(sA0);
        Qs(1:sPk(1,2),1:sPk(1,2)) = Krk;
        Qs(sPk(1,2)+1:sA0(1,1),sPk(1,2)+1:sA0(1,1)) = Krk(1:ndSig0,1:ndSig0);
        %}
        Rs = diag([10^-3, 10^-3]);
        [K0,Ss,eS] = lqr(A0,B0,Qs,Rs,Ns);
        info.Ss = Ss;
        info.eS = eS;
        %(2) Dinâmica rápida:
        sA22 = size(A22);
        Qf = eye(sA22);
        Qf(1:4,1:4) = 10^3*diag(ones(1,4));
        Qf(jPu1-ndSig0,jPu1-ndSig0) = 10^3;
        Qf(jPu2-ndSig0,jPu2-ndSig0) = 10^3;
        %Qf = Krk(ndSig0+1:sPk(1,2),ndSig0+1:sPk(1,2));
        Rf = diag([10^-3, 10^-3]); %Testar diferentes valores
        [K2,Sf,eF] = lqr(A22,B2,Qf,Rf,Nf);
        info.Sf = Sf;
        info.eF = eF;
        
        %Ganho total para projeto do regulador:
        K1 = K0 + K2*(A22\(A21 - B2*K0));
        Klqr = [K1, K2];
        eLQR = eig(Ab - Bb*Klqr);
        %{
        ReaLeLQR = real(eLQR);
        idRealeLQR = double(ReaLeLQR>0);
        nRealeLQR = sum(idRealeLQR);
        maxRealeLQR = max(abs(ReaLeLQR));
        minRealeLQR = min(abs(ReaLeLQR));
        log10dOG = log10(maxRealeLQR) - log10(minRealeLQR);
        %ok!
        %}
        info.Klqr = Klqr;
        info.eLQR = eLQR;
        
        %Novo sistema, com aplicação de regulador sobre a dinâmica rápida:
        %Obs.: Zp = Abs*Z + Bbs*Us
        Idz = eye(nSigTilk-ndSig0);
        Abs = Ab - Bb*K2*[A22\A21, Idz];
        IdnSig0 = eye(ndSig0);
        Bbs = Bb*(IdnSig0 - K2*(A22\B2));
        %Teste:
        eAbs = eig(Abs);
        info.eAbs = eAbs;
        %{
        ReaLeAbs = real(eAbs);
        idRealeAbs = double(ReaLeAbs>0);
        nRealeAbs = sum(idRealeAbs);
        maxRealeAbs = max(abs(ReaLeAbs));
        minRealeAbs = min(abs(ReaLeAbs));
        log10dOG = log10(maxRealeAbs) - log10(minRealeAbs);
        ad = 1;
        %ok!
        %}
        
        %Saídas:
        %{
        sPk = size(Pk);
        jPtXNM1 = Ng-rEstr;
        C1 = [zeros(1,sPk(1,2)),Pk(jPtXNM1,:)];
        jPuNM1 = Ng-rEstr-1;
        C2 = [Pk(jPuNM1,:),zeros(1,sPk(1,2))];
        %}
        %{
        sPk = size(Pk);
        jPtX1 = 2;
        C1 = [zeros(1,sPk(1,2)),Pk(jPtX1,:)];
        jPu1 = 1;
        C2 = [Pk(jPu1,:),zeros(1,sPk(1,2))];
        %}
        %{
        sAb = size(Ab);
        IdAb = eye(sAb);
        Axz = IdAb(:,1:nSigTilk+ndSig0);
        Azz = IdAb(:,nSigTilk+ndSig0+1:2*nSigTilk);
        C1xz = C1*Axz;
        C1zz = C1*Azz;
        C2xz = C2*Axz;
        C2zz = C2*Azz;
        %}
        
        %Análise linear do sistema "Planta + Regulador":
        testeSS = 1;
        if testeSS>0
            nAb = length(Ab);
            Cb = eye(nAb);
            sysP = ss(Ab,Bb,Cb,[]);
            sysK = ss([],[],[],Klqr);
            %{
            sBb = size(Bb);
            feedin = 1:sBb(1,2);
            feedout = 1:nAb;
            sys = feedback(sysP,sysK,feedin,feedout,-1);
            %}
            sys = feedback(sysP,sysK);
            t = 0:0.0001:50;
            [y,t] = impulse(sys,t);
            figure
            subplot(2,3,1)
            plot(t,y(:,3,1))
            grid
            title('modo 1')
            xlabel('Tempo (s)')
            subplot(2,3,2)
            plot(t,y(:,2,1))
            grid
            title('modo 2')
            xlabel('Tempo (s)')
            subplot(2,3,3)
            plot(t,y(:,3,1))
            grid
            title('modo 3')
            xlabel('Tempo (s)')
            subplot(2,3,4)
            plot(t,y(:,4,1))
            grid
            title('modo 4')
            xlabel('Tempo (s)')
            subplot(2,3,5)
            plot(t,y(:,jPu1-ndSig0,1))
            grid
            title('modo u1')
            xlabel('Tempo (s)')
            subplot(2,3,6)
            plot(t,y(:,jPu2-ndSig0,1))
            grid
            title('modo u2')
            xlabel('Tempo (s)')
            %{
            figure
            plot(t,y(:,sPk(1,2)-1,1))
            grid
            figure
            plot(t,y(:,sPk(1,2),1))
            grid
            %}
            %{
            figure
            plot(t,y(:,1,2))
            grid
            figure
            plot(t,y(:,2,2))
            grid
            %}
            %{
            figure
            plot(t,y(:,sPk(1,2)-1,2))
            grid
            figure
            plot(t,y(:,sPk(1,2),2))
            grid
            %}
        end
        
        beta0 = C1;
        beta1 = beta0*Abs;
        beta0Bbs = beta0*Bbs;
        beta2 = beta1*Abs;
        beta1Bbs = beta1*Bbs;
        beta3 = beta2*Abs;
        beta2Bbs = beta2*Bbs;
        beta4 = beta3*Abs;
        beta3Bbs = beta3*Bbs;
        beta5 = beta4*Abs;
        beta4Bbs = beta4*Bbs;
        beta6 = beta5*Abs;
        beta5Bbs = beta5*Bbs;
        beta7 = beta6*Abs;
        beta6Bbs = beta6*Bbs;
        beta8 = beta7*Abs;
        beta7Bbs = beta7*Bbs;
        beta9 = beta8*Abs;
        beta8Bbs = beta8*Bbs;
        beta10 = beta9*Abs;
        beta9Bbs = beta9*Bbs;
        beta11 = beta10*Abs;
        beta10Bbs = beta10*Bbs;
        
        kNM1 = argumentos.kNM1;
        alfa0 = kNM1*C2;
        alfa1 = alfa0*Abs;
        alfa0Bbs = alfa0*Bbs;
        alfa2 = alfa1*Abs;
        alfa1Bbs = alfa1*Bbs;
        alfa3 = alfa2*Abs;
        alfa2Bbs = alfa2*Bbs;
        alfa4 = alfa3*Abs;
        alfa3Bbs = alfa3*Bbs;
        
        ad = [beta0Bbs, beta1Bbs, beta2Bbs, beta3Bbs, beta4Bbs, beta5Bbs, beta6Bbs, beta7Bbs, beta8Bbs, beta9Bbs, beta10Bbs, alfa0Bbs, alfa1Bbs, alfa2Bbs, alfa3Bbs];
        
        %Armazenamento:
        ABvEps(1,1) = {A0}; ABvEps(1,2) = {B0}; 
        ABvEps(1,3) = {A21}; ABvEps(1,4) = {A22}; ABvEps(1,5) = {B2};
        ABvEps(1,6) = {varDif}; ABvEps(1,7) = {epS};
        ABvEps(1,8) = {K2};
        ABvEps(1,9) = {K0}; ABvEps(1,10) = {Klqr};
        
        %Obs.:
        %(1) Convergência não observada no vetor ad
        %(2) A função LQR deu problema para valores de N > 170
        %(3) Isso ocorreu muito prov. pelas diferenças de OG em Krk
        %(4) Linha de ação: separar dinâmicas mais uma vez
    case 5
        %%
        sPk = size(Pk);
        %Aproximação dos modos de corpo rígido:
        eps0 = 10^-8;
        sigTilk(1:ndSig0,1) = eps0*ones(ndSig0,1);
        %Redefinição de Krk, pela aproximação das freq. de CR:
        Krk = diag(sigTilk.');
        Ab0 = [zeros(sPk(1,2)), eye(sPk(1,2)); 
            -Krk, -Crk];
        Bb0 = Bb;
        %Redefinição de Ab, pela aproximação das freq. de CR:
        Ab = Ab0;
        
        %Primeira separação das dinâmicas, segundo os modos de vibração:
        %Obs.:
        %(1) Separação entre modos de CR e modos de vibração
        nSigTilk = length(sigTilk);
        epS = sigTilk(1,1)/sigTilk(ndSig0+1,1);
        A11 = Ab(1:nSigTilk+ndSig0,1:nSigTilk+ndSig0);
        A12 = Ab(1:nSigTilk+ndSig0,nSigTilk+ndSig0+1:2*nSigTilk);
        A21 = epS*Ab(nSigTilk+ndSig0+1:2*nSigTilk,1:nSigTilk+ndSig0);
        A22 = epS*Ab(nSigTilk+ndSig0+1:2*nSigTilk,nSigTilk+ndSig0+1:2*nSigTilk);
        B1 = Bb(1:nSigTilk+ndSig0,:);
        B2 = epS*Bb(nSigTilk+ndSig0+1:2*nSigTilk,:);
        B2 = B2 - B2.*double(abs(B2)<epsZero);
        A0 = A11 - A12*(A22\A21);
        B0 = B1 - A12*(A22\B2);
        varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
        
        %Ganhos parciais:
        Rs0 = diag([10^-2, 10^-2]);
        %(1) Dinâmica lenta:
        sA0 = size(A0);
        Qs = eye(sA0);
        Rs = Rs0;
        [K0,Ss,eS] = lqr(A0,B0,Qs,Rs,Ns);
        info.Ss = Ss;
        info.eS = eS;
        
        ABvEps(1,1) = {A0}; ABvEps(1,2) = {B0}; 
        ABvEps(1,3) = {A21}; ABvEps(1,4) = {A22}; ABvEps(1,5) = {B2};
        ABvEps(1,6) = {varDif}; ABvEps(1,7) = {epS};
        ABvEps(1,8) = {K0};
        
        %%
        %Segunda separação das dinâmicas, segundo os modos de vibração:
        %Obs.:
        %(1) Separação entre modos com OG 1 a 4, e modos com OG de 5 a 6
        Ab = A22;
        Bb = B2;
        nAb = length(Ab);
        idOG = double((-1 < OrdGrand) & (OrdGrand < 5));
        nOG = sum(idOG);
        epS = sigTilk(ndSig0+1,1)/sigTilk(ndSig0+nOG+1,1);
        A11 = Ab(1:nOG,1:nOG);
        A12 = Ab(1:nOG,nOG+1:nAb);
        A21 = epS*Ab(nOG+1:nAb,1:nOG);
        A22 = epS*Ab(nOG+1:nAb,nOG+1:nAb);
        B1 = Bb(1:nOG,:);
        B2 = epS*Bb(nOG+1:nAb,:);
        B2 = B2 - B2.*double(abs(B2)<epsZero);
        A0 = A11 - A12*(A22\A21);
        B0 = B1 - A12*(A22\B2);
        varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
        
        %Ganhos parciais:
        %(1) Dinâmica lenta:
        sA0 = size(A0);
        Qs = eye(sA0);
        Rs = Rs0;
        [K0,Ss,eS] = lqr(A0,B0,Qs,Rs,Ns);
        info.Ss = Ss;
        info.eS = eS;
        
        ABvEps(2,1) = {A0}; ABvEps(2,2) = {B0}; 
        ABvEps(2,3) = {A21}; ABvEps(2,4) = {A22}; ABvEps(2,5) = {B2};
        ABvEps(2,6) = {varDif}; ABvEps(2,7) = {epS};
        ABvEps(2,8) = {K0};
        
        %%
        %(2) Dinâmica rápida:
        sA22 = size(A22);
        Qf = eye(sA22);
        Rf = diag([10^-1, 10^-1]); %Testar diferentes valores
        [K2,Sf,eF] = lqr(A22,B2,Qf,Rf,Nf);
        info.Sf = Sf;
        info.eF = eF;
        
        %Ganho total para projeto do regulador:
        %{
        for j = 2:-1:1
            K0 = ABvEps{j,8};
            A22 = ABvEps{j,4};
            A21 = ABvEps{j,3};
            B2 = ABvEps{j,5};
            K1 = K0 + K2*(A22\(A21 - B2*K0));
            K2 = [K1, K2];
        end
        Klqr = K2;
        %}
        K0 = ABvEps{2,8};
        A22 = ABvEps{2,4};
        A21 = ABvEps{2,3};
        B2 = ABvEps{2,5};
        K1 = K0 + K2*(A22\(A21 - B2*K0));
        K2 = [K1, K2];
        
        K0 = ABvEps{1,8};
        A22 = ABvEps{1,4};
        A21 = ABvEps{1,3};
        B2 = ABvEps{1,5};
        K1 = K0 + K2*(A22\(A21 - B2*K0));
        Klqr0 = [K1, K2];
        
        %Ganho total para manter sistema em repouso:
        Klqr = Klqr0;
        eLQR = eig(Ab0 - Bb0*Klqr);
        %{
        ReaLeLQR = real(eLQR);
        idRealeLQR = double(ReaLeLQR>0);
        nRealeLQR = sum(idRealeLQR);
        maxRealeLQR = max(abs(ReaLeLQR));
        minRealeLQR = min(abs(ReaLeLQR));
        %ok!
        %}
        info.Klqr = Klqr;
        info.eLQR = eLQR;
        
        %%
        %Novo sistema, com aplicação de regulador sobre a dinâmica rápida:
        %Obs.: 
        %(1) Zp = Abs*Z + Bbs*Us
        %(2) Objetivo: definir lei de controle Us para sistema em movimento
        
        Idz = eye(nSigTilk-ndSig0);
        Abs = Ab - Bb*K2*[A22\A21, Idz];
        IdnSig0 = eye(ndSig0);
        Bbs = Bb*(IdnSig0 - K2*(A22\B2));
        
        %Saídas:
        sPk = size(Pk);
        jPtXNM1 = Ng-rEstr;
        C1 = [zeros(1,sPk(1,2)),Pk(jPtXNM1,:)];
        jPuNM1 = Ng-rEstr-1;
        C2 = [Pk(jPuNM1,:),zeros(1,sPk(1,2))];
        %{
        sAb = size(Ab);
        IdAb = eye(sAb);
        Axz = IdAb(:,1:nSigTilk+ndSig0);
        Azz = IdAb(:,nSigTilk+ndSig0+1:2*nSigTilk);
        C1xz = C1*Axz;
        C1zz = C1*Azz;
        C2xz = C2*Axz;
        C2zz = C2*Azz;
        %}
        
        beta0 = C1;
        beta1 = beta0*Abs;
        beta0Bbs = beta0*Bbs;
        beta2 = beta1*Abs;
        beta1Bbs = beta1*Bbs;
        beta3 = beta2*Abs;
        beta2Bbs = beta2*Bbs;
        beta4 = beta3*Abs;
        beta3Bbs = beta3*Bbs;
        
        kNM1 = argumentos.kNM1;
        alfa0 = kNM1*C2;
        alfa1 = alfa0*Abs;
        alfa0Bbs = alfa0*Bbs;
        alfa2 = alfa1*Abs;
        alfa1Bbs = alfa1*Bbs;
        alfa3 = alfa2*Abs;
        alfa2Bbs = alfa2*Bbs;
        alfa4 = alfa3*Abs;
        alfa3Bbs = alfa3*Bbs;
        
        ad = [beta0Bbs, beta1Bbs, beta2Bbs, beta3Bbs, alfa0Bbs, alfa1Bbs, alfa2Bbs, alfa3Bbs];
        
        %Armazenamento:
        ABvEps(1,1) = {A0}; ABvEps(1,2) = {B0}; 
        ABvEps(1,3) = {A21}; ABvEps(1,4) = {A22}; ABvEps(1,5) = {B2};
        ABvEps(1,6) = {varDif}; ABvEps(1,7) = {epS};
        ABvEps(1,8) = {K2};
        ABvEps(1,9) = {K0}; ABvEps(1,10) = {Klqr};
    case 6
        %%
        %PROJETO DO REGULADOR LQR:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jPtX0 = argumentos.jPtX0;
        jPu0 = argumentos.jPu0;
        jPtX1 = argumentos.jPtX1;
        jPu1 = argumentos.jPu1;
        jPtX2 = argumentos.jPtX2;
        jPu2 = argumentos.jPu2;
        jPv = argumentos.jPv;
        jPw = argumentos.jPw;
        %Aproximação dos modos de corpo rígido:%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %eps0 = 10^-6; %Funciona, com "warning": há Re(avl) ~ 10^-4 
        eps0 = 10^-4; %Funciona, sem "warning": cont alguns Re(avl) ~ 10^-4 
        sigTilk(jPtX0,1) = eps0;
        sigTilk(jPu0,1) = eps0;
        %Redefinição de Krk, pela aproximação das freq. de CR:%%%%%%%%%%%%%
        KrkLQR = diag(sigTilk.');
        %Redefinição de Ab, pela aproximação das freq. de CR:%%%%%%%%%%%%%%
        sPk = size(Pk);
        AbLQR = [zeros(sPk(1,2)), eye(sPk(1,2)); 
            -KrkLQR, -Crk];
        BbLQR = Bb;
        %Separação das dinâmicas:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        epS = sigTilk(1,1)/sigTilk(ndSig0+1,1);
        A11 = AbLQR(1:nSigTilk+ndSig0,1:nSigTilk+ndSig0);
        A12 = AbLQR(1:nSigTilk+ndSig0,nSigTilk+ndSig0+1:2*nSigTilk);
        A21 = epS*AbLQR(nSigTilk+ndSig0+1:2*nSigTilk,1:nSigTilk+ndSig0);
        A22 = epS*AbLQR(nSigTilk+ndSig0+1:2*nSigTilk,nSigTilk+ndSig0+1:2*nSigTilk);
        B1 = BbLQR(1:nSigTilk+ndSig0,:);
        B2 = epS*BbLQR(nSigTilk+ndSig0+1:2*nSigTilk,:);
        %B2 = B2 - B2.*double(abs(B2)<epsZero);
        A0 = A11 - A12*(A22\A21);
        B0 = B1 - A12*(A22\B2);
        %varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));
        %Ganhos parciais:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(1) Dinâmica lenta:
        sA0 = size(A0);
        Qs = eye(sA0);
        Rs = diag([10^-3, 10^-3]); %Testar diferentes valores
        [K0,Ss,eS] = lqr(A0,B0,Qs,Rs,Ns);
        info.Ss = Ss;
        info.eS = eS;
        %(2) Dinâmica rápida:
        sA22 = size(A22);
        Qf = eye(sA22);
        %Pesos em uTilrk0:
        Qf(jPtX1-ndSig0,jPtX1-ndSig0) = 10^3;
        Qf(jPtX2-ndSig0,jPtX2-ndSig0) = 10^3;
        Qf(jPv-ndSig0,jPv-ndSig0) = 10^3;
        Qf(jPw-ndSig0,jPw-ndSig0) = 10^3;
        Qf(jPu1-ndSig0,jPu1-ndSig0) = 10^3;
        Qf(jPu2-ndSig0,jPu2-ndSig0) = 10^3;
        Rf = diag([10^-3, 10^-3]); %Testar diferentes valores
        [K2,Sf,eF] = lqr(A22,B2,Qf,Rf,Nf);
        info.Sf = Sf;
        info.eF = eF;
        %Ganho total:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K1 = K0 + K2*(A22\(A21 - B2*K0));
        Klqr = [K1, K2];
        %Autalores do sistema:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(1) Planta (metamodelo) e Controlador (metamodelo):
        eLQR = eig(AbLQR - BbLQR*Klqr);
        %{
        ReaLeLQR = real(eLQR);
        idRealeLQR = double(ReaLeLQR>0);
        nRealeLQR = sum(idRealeLQR);
        maxRealeLQR = max(abs(ReaLeLQR));
        minRealeLQR = min(abs(ReaLeLQR));
        log10dOG = log10(maxRealeLQR) - log10(minRealeLQR);
        %ok!
        %}
        info.eLQR = eLQR;
        %(2) Planta (modelo sem aprox de CR) e Controlador (metamodelo):
        ePLQR = eig(Ab - Bb*Klqr);
        %{
        ReaLePLQR = real(ePLQR);
        idRealePLQR = double(ReaLePLQR>0);
        nRealePLQR = sum(idRealePLQR);
        maxRealePLQR = max(abs(ReaLePLQR));
        minRealePLQR = min(abs(ReaLePLQR));
        log10dOG = log10(maxRealePLQR) - log10(minRealePLQR);
        %ok!
        %}
        info.ePLQR = ePLQR;
        %Análise linear do sistema "Planta + Regulador":%%%%%%%%%%%%%%%%%%%
        testeSS = 1;
        %(1) Planta e Controlador: metamodelo
        if testeSS>0
            nAbLQR = length(AbLQR);
            sysK = ss([],[],[],Klqr);
            Cb = eye(nAbLQR);
            sysPb = ss(AbLQR,BbLQR,Cb,[]);
            sys = feedback(sysPb,sysK);
            %{
            sBb = size(Bb);
            feedin = 1:sBb(1,2);
            feedout = 1:nAb;
            sys = feedback(sysP,sysK,feedin,feedout,-1);
            %}
            t = 0:0.0001:50;
            [y,t] = impulse(sys,t);
            %Reposta a função impulso aplicada a Tteta:
            figure
            subplot(2,4,1)
            plot(t,y(:,jPtX1,1))
            grid
            title('modo 1')
            xlabel('Tempo (s)')
            subplot(2,4,2)
            plot(t,y(:,jPv,1))
            grid
            title('modo 2')
            xlabel('Tempo (s)')
            subplot(2,4,3)
            plot(t,y(:,jPw,1))
            grid
            title('modo 3')
            xlabel('Tempo (s)')
            subplot(2,4,4)
            plot(t,y(:,jPtX2,1))
            grid
            title('modo 4')
            xlabel('Tempo (s)')
            subplot(2,4,5)
            plot(t,y(:,jPu1,1))
            grid
            title('modo u1')
            xlabel('Tempo (s)')
            subplot(2,4,6)
            plot(t,y(:,jPu2,1))
            grid
            title('modo u2')
            xlabel('Tempo (s)')
            subplot(2,4,7)
            plot(t,y(:,jPu0,1))
            grid
            title('modo u0')
            xlabel('Tempo (s)')
            subplot(2,4,8)
            plot(t,y(:,jPtX0,1))
            grid
            title('modo tetaXp0')
            xlabel('Tempo (s)')
            %Reposta a função impulso aplicada a Tp:
            figure
            subplot(2,4,1)
            plot(t,y(:,jPtX1,2))
            grid
            title('modo 1')
            xlabel('Tempo (s)')
            subplot(2,4,2)
            plot(t,y(:,jPv,2))
            grid
            title('modo 2')
            xlabel('Tempo (s)')
            subplot(2,4,3)
            plot(t,y(:,jPw,2))
            grid
            title('modo 3')
            xlabel('Tempo (s)')
            subplot(2,4,4)
            plot(t,y(:,jPtX2,2))
            grid
            title('modo 4')
            xlabel('Tempo (s)')
            subplot(2,4,5)
            plot(t,y(:,jPu1,2))
            grid
            title('modo u1')
            xlabel('Tempo (s)')
            subplot(2,4,6)
            plot(t,y(:,jPu2,2))
            grid
            title('modo u2')
            xlabel('Tempo (s)')
            subplot(2,4,7)
            plot(t,y(:,jPu0,2))
            grid
            title('modo u0')
            xlabel('Tempo (s)')
            subplot(2,4,8)
            plot(t,y(:,jPtX0,2))
            grid
            title('modo tetaXp0')
            xlabel('Tempo (s)')
        end
        %(2) Planta (modelo sem aprox de CR) e Controlador (metamodelo):
        if testeSS>0
            nAb = length(Ab);
            Cb = eye(nAb);
            sysP = ss(Ab,Bb,Cb,[]);
            sysK = ss([],[],[],Klqr);
            sys = feedback(sysP,sysK);
            %{
            sBb = size(Bb);
            feedin = 1:sBb(1,2);
            feedout = 1:nAb;
            sys = feedback(sysP,sysK,feedin,feedout,-1);
            %}
            t = 0:0.0001:50;
            [y,t] = impulse(sys,t);
            %Reposta a função impulso aplicada a Tteta:
            figure
            subplot(2,4,1)
            plot(t,y(:,jPtX1,1))
            grid
            title('modo 1')
            xlabel('Tempo (s)')
            subplot(2,4,2)
            plot(t,y(:,jPv,1))
            grid
            title('modo 2')
            xlabel('Tempo (s)')
            subplot(2,4,3)
            plot(t,y(:,jPw,1))
            grid
            title('modo 3')
            xlabel('Tempo (s)')
            subplot(2,4,4)
            plot(t,y(:,jPtX2,1))
            grid
            title('modo 4')
            xlabel('Tempo (s)')
            subplot(2,4,5)
            plot(t,y(:,jPu1,1))
            grid
            title('modo u1')
            xlabel('Tempo (s)')
            subplot(2,4,6)
            plot(t,y(:,jPu2,1))
            grid
            title('modo u2')
            xlabel('Tempo (s)')
            subplot(2,4,7)
            plot(t,y(:,jPu0,1))
            grid
            title('modo u0')
            xlabel('Tempo (s)')
            subplot(2,4,8)
            plot(t,y(:,jPtX0,1))
            grid
            title('modo tetaXp0')
            xlabel('Tempo (s)')
            %Reposta a função impulso aplicada a Tp:
            figure
            subplot(2,4,1)
            plot(t,y(:,jPtX1,2))
            grid
            title('modo 1')
            xlabel('Tempo (s)')
            subplot(2,4,2)
            plot(t,y(:,jPv,2))
            grid
            title('modo 2')
            xlabel('Tempo (s)')
            subplot(2,4,3)
            plot(t,y(:,jPw,2))
            grid
            title('modo 3')
            xlabel('Tempo (s)')
            subplot(2,4,4)
            plot(t,y(:,jPtX2,2))
            grid
            title('modo 4')
            xlabel('Tempo (s)')
            subplot(2,4,5)
            plot(t,y(:,jPu1,2))
            grid
            title('modo u1')
            xlabel('Tempo (s)')
            subplot(2,4,6)
            plot(t,y(:,jPu2,2))
            grid
            title('modo u2')
            xlabel('Tempo (s)')
            subplot(2,4,7)
            plot(t,y(:,jPu0,2))
            grid
            title('modo u0')
            xlabel('Tempo (s)')
            subplot(2,4,8)
            plot(t,y(:,jPtX0,2))
            grid
            title('modo tetaXp0')
            xlabel('Tempo (s)')
        end
        
        %%
        %PROJETO DO CONTROLE POR SEGUIMENTO DE TRAJETORIA
        %Obs:
        %(1) Absorção de vibrações via REGULADOR LQR
        %(2) Seguimento de trajetória via MODOS DESLIZANTES
        
        %Matriz permuta de linhas Zb = Perl*Z e Z = Perl.'*Zb:
        sAb = size(Ab);
        IdAb = eye(sAb);
        Perl = [IdAb(1:ndSig0,:); 
            IdAb(sPk(1,2)+1:sPk(1,2)+ndSig0,:); 
            IdAb(ndSig0+1:sPk(1,2),:); 
            IdAb(sPk(1,2)+ndSig0+1:2*sPk(1,2),:)];
        info.Perl = Perl;
        %Aproximação dos modos de corpo rígido:%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        eps0 = 10^-8;
        sigTilk(jPtX0,1) = eps0;
        sigTilk(jPu0,1) = eps0;
        %Redefinição de Krk, pela aproximação das freq. de CR:%%%%%%%%%%%%%
        KrkLQR = diag(sigTilk.');
        %Redefinição de Ab, pela aproximação das freq. de CR:%%%%%%%%%%%%%%
        sPk = size(Pk);
        AbLQR = [zeros(sPk(1,2)), eye(sPk(1,2)); 
            -KrkLQR, -Crk];
        %Equação de dinâmica com rearranjo: Zbp = Abb*Zb + Bbb*Uc
        AbbLQR = Perl*AbLQR*Perl.';
        BbbLQR = Perl*BbLQR;
        %Dinâmica rápida (modos de vibração):
        A21 = epS*AbbLQR(2*ndSig0+1:2*sPk(1,2),1:2*ndSig0);
        A22 = epS*AbbLQR(2*ndSig0+1:2*sPk(1,2),2*ndSig0+1:2*sPk(1,2));
        B2 = epS*BbbLQR(2*ndSig0+1:2*sPk(1,2),:);
        %Ganhos parciais:
        %(1) Dinâmica rápida:
        sA22 = size(A22);
        Qbf = eye(sA22);
        %Pesos em uTilrk0:
        Qbf(jPtX1-ndSig0,jPtX1-ndSig0) = 10^3;
        Qbf(jPtX2-ndSig0,jPtX2-ndSig0) = 10^3;
        Qbf(jPv-ndSig0,jPv-ndSig0) = 10^3;
        Qbf(jPw-ndSig0,jPw-ndSig0) = 10^3;
        Qbf(jPu1-ndSig0,jPu1-ndSig0) = 10^3;
        Qbf(jPu2-ndSig0,jPu2-ndSig0) = 10^3;
        %Pesos em vTilrk0:
        nrkv = sPk(1,2) - ndSig0;
        Qbf(nrkv+jPtX1-ndSig0,nrkv+jPtX1-ndSig0) = 10^3;
        Qbf(nrkv+jPtX2-ndSig0,nrkv+jPtX2-ndSig0) = 10^3;
        Qbf(nrkv+jPv-ndSig0,nrkv+jPv-ndSig0) = 10^3;
        Qbf(nrkv+jPw-ndSig0,nrkv+jPw-ndSig0) = 10^3;
        Qbf(nrkv+jPu1-ndSig0,nrkv+jPu1-ndSig0) = 10^3;
        Qbf(nrkv+jPu2-ndSig0,nrkv+jPu2-ndSig0) = 10^3;
        Rbf = diag([10^-3, 10^-3]); %Testar diferentes valores
        [K2,Sf,eF] = lqr(A22,B2,Qbf,Rbf,Nf);
        info.Sf = Sf;
        info.eF = eF;
        %Novo sistema, com aplicação de regulador sobre a dinâmica rápida:
        %Obs.: Zp = Abs*Z + Bbs*Us
        Idz = eye(2*(nSigTilk-ndSig0));
        Abs = AbbLQR - BbbLQR*K2*[A22\A21, Idz];
        IdnSig0 = eye(ndSig0);
        Bbs = BbbLQR*(IdnSig0 - K2*(A22\B2));
        %Teste:
        eAbs = eig(Abs);
        info.eAbs = eAbs;
        %{
        ReaLeAbs = real(eAbs);
        idRealeAbs = double(ReaLeAbs>0);
        nRealeAbs = sum(idRealeAbs);
        maxRealeAbs = max(abs(ReaLeAbs));
        minRealeAbs = min(abs(ReaLeAbs));
        log10dOG = log10(maxRealeAbs) - log10(minRealeAbs);
        ad = 1;
        %ok!
        %}
        
        %Saídas:
        C1 = [zeros(1,sPk(1,2)),Pk(2,:)]*Perl.';
        C2 = [Pk(1,:),zeros(1,sPk(1,2))]*Perl.';
        
        beta0 = C1;
        beta1 = beta0*Abs;
        beta0Bbs = beta0*Bbs;
        %{
        beta2 = beta1*Abs;
        beta1Bbs = beta1*Bbs;
        beta3 = beta2*Abs;
        beta2Bbs = beta2*Bbs;
        beta4 = beta3*Abs;
        beta3Bbs = beta3*Bbs;
        beta5 = beta4*Abs;
        beta4Bbs = beta4*Bbs;
        beta6 = beta5*Abs;
        beta5Bbs = beta5*Bbs;
        beta7 = beta6*Abs;
        beta6Bbs = beta6*Bbs;
        beta8 = beta7*Abs;
        beta7Bbs = beta7*Bbs;
        beta9 = beta8*Abs;
        beta8Bbs = beta8*Bbs;
        beta10 = beta9*Abs;
        beta9Bbs = beta9*Bbs;
        beta11 = beta10*Abs;
        beta10Bbs = beta10*Bbs;
        %}
        alfa0 = C2;
        alfa1 = alfa0*Abs;
        alfa2 = alfa1*Abs;
        alfa1Bbs = alfa1*Bbs;
        %{
        alfa0Bbs = alfa0*Bbs;
        alfa3 = alfa2*Abs;
        alfa2Bbs = alfa2*Bbs;
        alfa4 = alfa3*Abs;
        alfa3Bbs = alfa3*Bbs;
        ad = [beta0Bbs, beta1Bbs, beta2Bbs, beta3Bbs, beta4Bbs, beta5Bbs, beta6Bbs, beta7Bbs, beta8Bbs, beta9Bbs, beta10Bbs, alfa0Bbs, alfa1Bbs, alfa2Bbs, alfa3Bbs];
        %}
        
        %Armazenamento:
        %(1) Ganho do regulador:
        ABvEps(1,1) = {Klqr};
        %(2) Matrizes da dinâmica, com absorvedor de vibração:
        ABvEps(1,2) = {Abs}; ABvEps(1,3) = {Bbs};
        %(3) Matrizes da dinâmica externa:
        ABvEps(1,4) = {beta0}; ABvEps(1,5) = {beta1}; 
        ABvEps(1,6) = {beta0Bbs};
        ABvEps(1,7) = {alfa0}; ABvEps(1,8) = {alfa1}; ABvEps(1,9) = {alfa2}; 
        ABvEps(1,10) = {alfa1Bbs};
    otherwise
    disp('other value')
end
%}
%%
%PROJETO DO REGULADOR LQR:

%Aproximação dos modos de corpo rígido:
epsTx0 = 10^-4;
epsU0 = 10^-4;
sigTilk(jPtX0,1) = epsTx0;
sigTilk(jPu0,1) = epsU0;
%Redefinição de Krk, pela aproximação das freq. de CR:
KrkLQR = diag(sigTilk.');
%Redefinição de Ab, pela aproximação das freq. de CR:
sPk = size(Pk);
AbLQR = [zeros(sPk(1,2)), eye(sPk(1,2)); 
    -KrkLQR, -Crk];
BbLQR = Bb;

%%
%Separação das dinâmicas:
%Obs.:
%(1) Primeira sepação entre as dinâmicas de corpo rígido e de vibração
%(2) Matriz inclui bloco referentes a urk, e linhas e colunas ref. a vrk0
epS = sigTilk(1,1)/sigTilk(ndSig0+1,1);
A11 = AbLQR(1:nSigTilk+ndSig0,1:nSigTilk+ndSig0);
A12 = AbLQR(1:nSigTilk+ndSig0,nSigTilk+ndSig0+1:2*nSigTilk);
A21 = epS*AbLQR(nSigTilk+ndSig0+1:2*nSigTilk,1:nSigTilk+ndSig0);
A22 = epS*AbLQR(nSigTilk+ndSig0+1:2*nSigTilk,nSigTilk+ndSig0+1:2*nSigTilk);
B1 = BbLQR(1:nSigTilk+ndSig0,:);
B2 = epS*BbLQR(nSigTilk+ndSig0+1:2*nSigTilk,:);
%B2 = B2 - B2.*double(abs(B2)<epsZero);
A0 = A11 - A12*(A22\A21);
B0 = B1 - A12*(A22\B2);
%varDif = @(Xs,Us)(-A22\(A21*Xs + B2*Us));

%%
%Ganhos parciais:
FtTx1 = 1;
FtTx2 = 1;
FtV = 1;
FtW = 1;
FtU1 = 10^1;
FtU2 = 10^1;

%(1) Dinâmica lenta:
sA0 = size(A0);
Qs = eye(sA0);
%
Qs(jPtX1,jPtX1) = FtTx1;
Qs(jPtX2,jPtX2) = FtTx2;
Qs(jPv,jPv) = FtV;
Qs(jPw,jPw) = FtW;
Qs(jPu1,jPu1) = FtU1;
Qs(jPu2,jPu2) = FtU2;
%Obs.:
%(1) O acréscimo de peso nos modos principais aumentou o módulo da parte 
%real máxima dos autovalores global (planta+regulador)
%(2) Não houve alteração no módulo da parte real mínima
%}

Rs = diag([10^-3, 10^-3]);
[K0,Ss,eS] = lqr(A0,B0,Qs,Rs,Ns);
info.Ss = Ss;
info.eS = eS;

%(2) Dinâmica rápida:
sA22 = size(A22);
Qf = eye(sA22);

%Pesos em uTilrkv:
Qf(jPtX1-ndSig0,jPtX1-ndSig0) = FtTx1;
Qf(jPtX2-ndSig0,jPtX2-ndSig0) = FtTx2;
Qf(jPv-ndSig0,jPv-ndSig0) = FtV;
Qf(jPw-ndSig0,jPw-ndSig0) = FtW;
Qf(jPu1-ndSig0,jPu1-ndSig0) = FtU1;
Qf(jPu2-ndSig0,jPu2-ndSig0) = FtU2;

%Pesos em Uf:
Rf = diag([10^-3, 10^-3]);

%Ganho:
[K2,Sf,eF] = lqr(A22,B2,Qf,Rf,Nf);
info.Sf = Sf;
info.eF = eF;

%Ganho total:
K1 = K0 + K2*(A22\(A21 - B2*K0));
Klqr = [K1, K2];

%Autalores do sistema:
%(1) Planta (metamodelo) e Controlador (metamodelo):
eLQR = eig(AbLQR - BbLQR*Klqr);
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
info.eLQR = eLQR;
%(2) Planta (modelo sem aprox de CR) e Controlador (metamodelo):
ePLQR = eig(Ab - Bb*Klqr);
%{
ReaLePLQR = real(ePLQR);
idRealePLQR = double(ReaLePLQR>0);
nRealePLQR = sum(idRealePLQR);
maxRealePLQR = max(abs(ReaLePLQR));
minRealePLQR = min(abs(ReaLePLQR));
log10dOG = log10(maxRealePLQR) - log10(minRealePLQR);
%ok!
%Obs.:
%(1) Todos os autovalores com parte real negativa
%(2) Módulo da menor parte real: ~1.37*10^-4
%09/06/2022
%}
info.ePLQR = ePLQR;

%%
%Análise linear do sistema "Planta + Regulador" (gráficos):
%{
testeSS = 0;
%(1) Planta e Controlador: metamodelo
if testeSS>0
    nAbLQR = length(AbLQR);
    sysK = ss([],[],[],Klqr);
    Cb = eye(nAbLQR);
    sysPb = ss(AbLQR,BbLQR,Cb,[]);
    sys = feedback(sysPb,sysK);
    %{
    sBb = size(Bb);
    feedin = 1:sBb(1,2);
    feedout = 1:nAb;
    sys = feedback(sysP,sysK,feedin,feedout,-1);
    %}
    t = 0:0.0001:50;
    [y,t] = impulse(sys,t);
    %Reposta a função impulso aplicada a Tteta:
    figure
    subplot(2,4,1)
    plot(t,y(:,jPtX1,1))
    grid
    title('modo 1')
    xlabel('Tempo (s)')
    subplot(2,4,2)
    plot(t,y(:,jPv,1))
    grid
    title('modo 2')
    xlabel('Tempo (s)')
    subplot(2,4,3)
    plot(t,y(:,jPw,1))
    grid
    title('modo 3')
    xlabel('Tempo (s)')
    subplot(2,4,4)
    plot(t,y(:,jPtX2,1))
    grid
    title('modo 4')
    xlabel('Tempo (s)')
    subplot(2,4,5)
    plot(t,y(:,jPu1,1))
    grid
    title('modo u1')
    xlabel('Tempo (s)')
    subplot(2,4,6)
    plot(t,y(:,jPu2,1))
    grid
    title('modo u2')
    xlabel('Tempo (s)')
    subplot(2,4,7)
    plot(t,y(:,jPu0,1))
    grid
    title('modo u0')
    xlabel('Tempo (s)')
    subplot(2,4,8)
    plot(t,y(:,jPtX0,1))
    grid
    title('modo tetaXp0')
    xlabel('Tempo (s)')
    %Reposta a função impulso aplicada a Tp:
    figure
    subplot(2,4,1)
    plot(t,y(:,jPtX1,2))
    grid
    title('modo 1')
    xlabel('Tempo (s)')
    subplot(2,4,2)
    plot(t,y(:,jPv,2))
    grid
    title('modo 2')
    xlabel('Tempo (s)')
    subplot(2,4,3)
    plot(t,y(:,jPw,2))
    grid
    title('modo 3')
    xlabel('Tempo (s)')
    subplot(2,4,4)
    plot(t,y(:,jPtX2,2))
    grid
    title('modo 4')
    xlabel('Tempo (s)')
    subplot(2,4,5)
    plot(t,y(:,jPu1,2))
    grid
    title('modo u1')
    xlabel('Tempo (s)')
    subplot(2,4,6)
    plot(t,y(:,jPu2,2))
    grid
    title('modo u2')
    xlabel('Tempo (s)')
    subplot(2,4,7)
    plot(t,y(:,jPu0,2))
    grid
    title('modo u0')
    xlabel('Tempo (s)')
    subplot(2,4,8)
    plot(t,y(:,jPtX0,2))
    grid
    title('modo tetaXp0')
    xlabel('Tempo (s)')
end

%(2) Planta (modelo sem aprox de CR) e Controlador (metamodelo):
nAb = length(Ab);
if testeSS>0
    Cb = eye(nAb);
    sysP = ss(Ab,Bb,Cb,[]);
    sysK = ss([],[],[],Klqr);
    sys = feedback(sysP,sysK);
    %{
    sBb = size(Bb);
    feedin = 1:sBb(1,2);
    feedout = 1:nAb;
    sys = feedback(sysP,sysK,feedin,feedout,-1);
    %}
    t = 0:0.0001:50;
    [y,t] = impulse(sys,t);
    %Reposta a função impulso aplicada a Tteta:
    figure
    subplot(2,4,1)
    plot(t,y(:,jPtX1,1))
    grid
    title('modo 1')
    xlabel('Tempo (s)')
    subplot(2,4,2)
    plot(t,y(:,jPv,1))
    grid
    title('modo 2')
    xlabel('Tempo (s)')
    subplot(2,4,3)
    plot(t,y(:,jPw,1))
    grid
    title('modo 3')
    xlabel('Tempo (s)')
    subplot(2,4,4)
    plot(t,y(:,jPtX2,1))
    grid
    title('modo 4')
    xlabel('Tempo (s)')
    subplot(2,4,5)
    plot(t,y(:,jPu1,1))
    grid
    title('modo u1')
    xlabel('Tempo (s)')
    subplot(2,4,6)
    plot(t,y(:,jPu2,1))
    grid
    title('modo u2')
    xlabel('Tempo (s)')
    subplot(2,4,7)
    plot(t,y(:,jPu0,1))
    grid
    title('modo u0')
    xlabel('Tempo (s)')
    subplot(2,4,8)
    plot(t,y(:,jPtX0,1))
    grid
    title('modo tetaXp0')
    xlabel('Tempo (s)')
    %Reposta a função impulso aplicada a Tp:
    figure
    subplot(2,4,1)
    plot(t,y(:,jPtX1,2))
    grid
    title('modo 1')
    xlabel('Tempo (s)')
    subplot(2,4,2)
    plot(t,y(:,jPv,2))
    grid
    title('modo 2')
    xlabel('Tempo (s)')
    subplot(2,4,3)
    plot(t,y(:,jPw,2))
    grid
    title('modo 3')
    xlabel('Tempo (s)')
    subplot(2,4,4)
    plot(t,y(:,jPtX2,2))
    grid
    title('modo 4')
    xlabel('Tempo (s)')
    subplot(2,4,5)
    plot(t,y(:,jPu1,2))
    grid
    title('modo u1')
    xlabel('Tempo (s)')
    subplot(2,4,6)
    plot(t,y(:,jPu2,2))
    grid
    title('modo u2')
    xlabel('Tempo (s)')
    subplot(2,4,7)
    plot(t,y(:,jPu0,2))
    grid
    title('modo u0')
    xlabel('Tempo (s)')
    subplot(2,4,8)
    plot(t,y(:,jPtX0,2))
    grid
    title('modo tetaXp0')
    xlabel('Tempo (s)')
end
%}

%%
%PROJETO DO CONTROLE POR SEGUIMENTO DE TRAJETORIA
%Obs:
%(1) Absorção de vibrações via REGULADOR LQR
%(2) Seguimento de trajetória via MODOS DESLIZANTES

%Matriz permuta de linhas Zb = Perl*Z e Z = Perl.'*Zb:
sAb = size(Ab);
IdAb = eye(sAb);
Perl = [IdAb(1:ndSig0,:); 
    IdAb(sPk(1,2)+1:sPk(1,2)+ndSig0,:); 
    IdAb(ndSig0+1:sPk(1,2),:); 
    IdAb(sPk(1,2)+ndSig0+1:2*sPk(1,2),:)];
info.Perl = Perl;

%%
%Aproximação dos modos de corpo rígido:
epsTx0 = 10^-3;
epsU0 = 10^-3;
sigTilk(jPtX0,1) = epsTx0;
sigTilk(jPu0,1) = epsU0;

%%
%Redefinição de Krk, pela aproximação das freq. de CR:
KrkLQR = diag(sigTilk.');

%Redefinição de Ab, pela aproximação das freq. de CR:
sPk = size(Pk);
AbLQR = [zeros(sPk(1,2)), eye(sPk(1,2)); 
    -KrkLQR, -Crk];

%Equação de dinâmica com rearranjo: Zbp = Abb*Zb + Bbb*Uc
AbbLQR = Perl*AbLQR*Perl.';
BbbLQR = Perl*BbLQR;

%Dinâmica rápida (modos de vibração):
%epS = sigTilk(1,1)/sigTilk(ndSig0+1,1);
epS = sigTilk(1,1)/sigTilk(jPu1,1);
A21 = epS*AbbLQR(2*ndSig0+1:2*sPk(1,2),1:2*ndSig0);
A22 = epS*AbbLQR(2*ndSig0+1:2*sPk(1,2),2*ndSig0+1:2*sPk(1,2));
B2 = epS*BbbLQR(2*ndSig0+1:2*sPk(1,2),:);

%Ganhos parciais:
%(1) Dinâmica rápida:
sA22 = size(A22);
Qbf = eye(sA22);
%(1.a) Pesos em uTilrkv:
Qbf(jPtX1-ndSig0,jPtX1-ndSig0) = FtTx1;
Qbf(jPtX2-ndSig0,jPtX2-ndSig0) = FtTx2;
Qbf(jPv-ndSig0,jPv-ndSig0) = FtV;
Qbf(jPw-ndSig0,jPw-ndSig0) = FtW;
Qbf(jPu1-ndSig0,jPu1-ndSig0) = FtU1;
Qbf(jPu2-ndSig0,jPu2-ndSig0) = FtU2;
%(1.b) Pesos em vTilrkv:
nrkv = sPk(1,2) - ndSig0;
Qbf(nrkv+jPtX1-ndSig0,nrkv+jPtX1-ndSig0) = FtTx1;
Qbf(nrkv+jPtX2-ndSig0,nrkv+jPtX2-ndSig0) = FtTx2;
Qbf(nrkv+jPv-ndSig0,nrkv+jPv-ndSig0) = FtV;
Qbf(nrkv+jPw-ndSig0,nrkv+jPw-ndSig0) = FtW;
Qbf(nrkv+jPu1-ndSig0,nrkv+jPu1-ndSig0) = FtU1;
Qbf(nrkv+jPu2-ndSig0,nrkv+jPu2-ndSig0) = FtU2;
%(1.c) Pesos em Uf:
Rbf = diag([10^-3, 10^-3]);

%Ganho:
[K2,Sf,eF] = lqr(A22,B2,Qbf,Rbf,Nf);
info.Sf = Sf;
info.eF = eF;

%Novo sistema, com aplicação de regulador sobre a dinâmica rápida:
%Obs.: Zp = Abs*Z + Bbs*Us
Idz = eye(2*(nSigTilk-ndSig0));
Abs = AbbLQR - BbbLQR*K2*[A22\A21, Idz];
IdnSig0 = eye(ndSig0);
Bbs = BbbLQR*(IdnSig0 - K2*(A22\B2));

%Teste:
eAbs = eig(Abs);
info.eAbs = eAbs;
%{
ReaLeAbs = real(eAbs);
idRealeAbs = double(ReaLeAbs>0);
nRealeAbs = sum(idRealeAbs);
maxRealeAbs = max(abs(ReaLeAbs));
minRealeAbs = min(abs(ReaLeAbs));
log10dOG = log10(maxRealeAbs) - log10(minRealeAbs);
ad = 1;
%ok!
%Obs.:
%(1) Todos os autovalores apresentaram parte real negativa
%(2) Alguns resultados interessantes:
%    (a) eps0 = 10^-1: minRealeAbs = 1.37*10^-4
%    (b) eps0 = 10^-2: minRealeAbs = 1.37*10^-4
%    (c) eps0 = 10^-3: minRealeAbs = 1.37*10^-4
%    (d) eps0 = 10^-4: minRealeAbs = 1.53*10^-5
%09/06/2022
%}

%%
%PROJETO DO CONTROLE POR SEGUIMENTO DE TRAJETORIA
%Obs:
%(1) Absorção de vibrações via REGULADOR LQR
%(2) Seguimento de trajetória via MODOS DESLIZANTES
%(3) Quebra das dinâmicas de vibração

%Matriz permuta de linhas Zbb = Perlqv*Z e Z = Perlqv.'*Zbb:
sAb = size(Ab);
IdAb = eye(sAb);
Perlqv = [IdAb(1:ndSig0,:);
    IdAb(sPk(1,2)+1:sPk(1,2)+ndSig0,:);
    IdAb(jPtX1:jPtX2,:);
    IdAb(sPk(1,2)+jPtX1:sPk(1,2)+jPtX2,:);
    IdAb(jPtX2+1:sPk(1,2),:);
    IdAb(sPk(1,2)+jPtX2+1:2*sPk(1,2),:)];
info.Perlqv = Perlqv;

%%
%Aproximação dos modos de corpo rígido:
epsTx0 = 10^-3;
epsU0 = 10^-3;
sigTilk(jPtX0,1) = epsTx0;
sigTilk(jPu0,1) = epsU0;

%%
%Redefinição de Krk, pela aproximação das freq. de CR:
KrkLQR = diag(sigTilk.');

%Redefinição de Ab, pela aproximação das freq. de CR:
sPk = size(Pk);
AbLQR = [zeros(sPk(1,2)), eye(sPk(1,2)); 
    -KrkLQR, -Crk];

%Equação de dinâmica com rearranjo: Zbp = Abb*Zb + Bbb*Uc
AbbLQR = Perlqv*AbLQR*Perlqv.';
BbbLQR = Perlqv*BbLQR;

%%
%Separação das dinâmicas:
%Obs.:
%(1) Primeira separação entre as dinâmicas: corpo rígido e vibração
epS1 = sigTilk(1,1)/sigTilk(jPtX1,1);
Abb111 = AbbLQR(1:2*ndSig0,1:2*ndSig0);
Abb121 = AbbLQR(1:2*ndSig0,2*ndSig0+1:2*nSigTilk);
Abb211 = AbbLQR(2*ndSig0+1:2*nSigTilk,1:2*ndSig0);
Abb221 = AbbLQR(2*ndSig0+1:2*nSigTilk,2*ndSig0+1:2*nSigTilk);
Bbb11 = BbbLQR(1:2*ndSig0,:);
Bbb21 = BbbLQR(2*ndSig0+1:2*nSigTilk,:);
%(1a) xpS1 = A01*xS1 + B01*Us1 e epS1*zpF1 = A221*zF1 + B21*Uf1:
A111 = Abb111;
A121 = Abb121;
A211 = epS1*Abb211;
A221 = epS1*Abb221;
B11 = Bbb11;
B21 = epS1*Bbb21;
A01 = A111 - A121*(A221\A211);
B01 = B11 - A121*(A221\B21);
fzS1 = @(xS1,Us1)(-A221\(A211*xS1 + B21*Us1));
%(2) Segunda separação entre as dinâmicas: torção + reversos e demais
sA221 = size(A221);
epS2 = sigTilk(jPtX1,1)/sigTilk(jPu1,1);
Abb112 = A221(1:2*(jPtX2-ndSig0),1:2*(jPtX2-ndSig0));
Abb122 = A221(1:2*(jPtX2-ndSig0),2*(jPtX2-ndSig0)+1:sA221(1,1));
Abb212 = A221(2*(jPtX2-ndSig0)+1:sA221(1,1),1:2*(jPtX2-ndSig0));
Abb222 = A221(2*(jPtX2-ndSig0)+1:sA221(1,1),2*(jPtX2-ndSig0)+1:sA221(1,1));
Bbb12 = B21(1:2*(jPtX2-ndSig0),:);
Bbb22 = B21(2*(jPtX2-ndSig0)+1:sA221(1,1),:);
%(2a) xpS2 = A02*xS2 + B02*Us2 e epS2*zpF2 = A222*zF2 + B22*Uf2:
A112 = Abb112;
A122 = Abb122;
A212 = epS2*Abb212;
A222 = epS2*Abb222;
B12 = Bbb12;
B22 = epS2*Bbb22;
A02 = A112 - A122*(A222\A212);
B02 = B12 - A122*(A222\B22);
fzS2 = @(xS2,Us2)(-A222\(A212*xS2 + B22*Us2));

%%
%Ganhos parciais:
%(1) Dinâmica rápida da segunda separação:
sA222 = size(A222);
Qbf22 = eye(sA222);
%(1.a) Pesos em uTilrkv2:
Qbf22(jPu1-jPtX2,jPu1-jPtX2) = FtU1;
Qbf22(jPu2-jPtX2,jPu2-jPtX2) = FtU2;
%(1.b) Pesos em vTilrkv2:
nrkv = jPu2 - jPtX2;
Qbf22(nrkv+jPu1-jPtX2,nrkv+jPu1-jPtX2) = FtU1;
Qbf22(nrkv+jPu2-jPtX2,nrkv+jPu2-jPtX2) = FtU2;
%(1.c) Pesos em Uf:
Rbf22 = diag([10^-3, 10^-3]);
%(1.d) Ganho: Uf2 = -K22*zF2
[K22,Sf22,eF22] = lqr(A222,B22,Qbf22,Rbf22,Nf);
info.Sf22 = Sf22;
info.eF22 = eF22;
%(2) Dinâmica lenta da segunda separação:
sA02 = size(A02);
Qbf02 = eye(sA02);
%(2.a) Pesos em uTilrkv1:
Qbf02(jPtX1-ndSig0,jPtX1-ndSig0) = FtTx1;
Qbf02(jPtX2-ndSig0,jPtX2-ndSig0) = FtTx2;
Qbf02(jPv-ndSig0,jPv-ndSig0) = FtV;
Qbf02(jPw-ndSig0,jPw-ndSig0) = FtW;
%(2.b) Pesos em vTilrkv1:
nrkv = jPtX2 - ndSig0;
Qbf02(nrkv+jPtX1-ndSig0,nrkv+jPtX1-ndSig0) = FtTx1;
Qbf02(nrkv+jPtX2-ndSig0,nrkv+jPtX2-ndSig0) = FtTx2;
Qbf02(nrkv+jPv-ndSig0,nrkv+jPv-ndSig0) = FtV;
Qbf02(nrkv+jPw-ndSig0,nrkv+jPw-ndSig0) = FtW;
%(2.c) Pesos em Uf:
Rbf02 = diag([10^-3, 10^-3]);
%(2.d) Ganho: Us2 = -K02*xS2
[K02,Sf02,eF02] = lqr(A02,B02,Qbf02,Rbf02,Nf);
info.Sf02 = Sf02;
info.eF02 = eF02;

%%
%Ganho total para regulador dos modos de vibração: Uf1 = - K21*zF1
K21 = [K02 + K22*(A222\(A212 - B22*K02)), K22];

%Novo sistema, com aplicação de regulador sobre a dinâmica rápida:
%Obs.: Zbbp = Abbs*Zbb + Bbbs*Us
Iuc = eye(2);
Kz = [K21*(A221\A211), K21];
Kuc = Iuc - K21*(A221\B21);
Abbs = AbbLQR - BbbLQR*Kz;
Bbbs = BbbLQR*Kuc;

%Teste:
eAbbs = eig(Abbs);
info.eAbbs = eAbbs;
%{
ReaLeAbs = real(eAbbs);
idRealeAbs = double(ReaLeAbs>0);
nRealeAbs = sum(idRealeAbs);
maxRealeAbs = max(abs(ReaLeAbs));
minRealeAbs = min(abs(ReaLeAbs));
log10dOG = log10(maxRealeAbs) - log10(minRealeAbs);
ad = 1;
%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FORMA NORMAL COM LQR
%Obs.: Zb é o vetor original de estado

%Saídas:
C1s = [zeros(1,sPk(1,2)),Pk(2,:)]*Perl.';
C2s = [Pk(1,:),zeros(1,sPk(1,2))]*Perl.';

%Dinâmica externa (sistema com LQR):
%Obs.: yTil = aTils(X) + bTilb*Us
beta0s = C1s;
beta1s = beta0s*Abs;
beta0Bbs = beta0s*Bbs;
%{
beta2 = beta1*Abs;
beta1Bbs = beta1*Bbs;
beta3 = beta2*Abs;
beta2Bbs = beta2*Bbs;
beta4 = beta3*Abs;
beta3Bbs = beta3*Bbs;
beta5 = beta4*Abs;
beta4Bbs = beta4*Bbs;
beta6 = beta5*Abs;
beta5Bbs = beta5*Bbs;
beta7 = beta6*Abs;
beta6Bbs = beta6*Bbs;
beta8 = beta7*Abs;
beta7Bbs = beta7*Bbs;
beta9 = beta8*Abs;
beta8Bbs = beta8*Bbs;
beta10 = beta9*Abs;
beta9Bbs = beta9*Bbs;
beta11 = beta10*Abs;
beta10Bbs = beta10*Bbs;
%}
alfa0s = C2s;
alfa1s = alfa0s*Abs;
alfa2s = alfa1s*Abs;
alfa1Bbs = alfa1s*Bbs;
%{
alfa0Bbs = alfa0*Bbs;
alfa3 = alfa2*Abs;
alfa2Bbs = alfa2*Bbs;
alfa4 = alfa3*Abs;
alfa3Bbs = alfa3*Bbs;
ad1 = [beta0Bbs, beta1Bbs, beta2Bbs, beta3Bbs, beta4Bbs, beta5Bbs, beta6Bbs, beta7Bbs, beta8Bbs, beta9Bbs, beta10Bbs, alfa0Bbs, alfa1Bbs, alfa2Bbs, alfa3Bbs];
ad2 = [beta0Bbs, alfa0Bbs, alfa1Bbs];
%}

%Dinâmica interna (sistema com LQR):
%Obs.:
%(1) mi = Ami*Zb
%(2) mip = GamaMi*mi + SigmaMi*yTil + Ami*Bbs*Us
[maxBb1,jM1] = max(abs(Bbs(:,1)));
[maxBb2,jM2] = max(abs(Bbs(:,2)));
info.maxBb1s = maxBb1;
info.maxBb2s = maxBb2;
bn11 = Bbs(jM1,1);
bn22 = Bbs(jM2,2);
bn21 = Bbs(jM2,1);
bn12 = Bbs(jM1,2);
Bn = [bn11, bn21; bn12, bn22];
invBn = Bn\eye(2);
ibn11 = invBn(1,1);
ibn22 = invBn(2,2);
ibn21 = invBn(2,1);
ibn12 = invBn(1,2);
AmiTotals = zeros(nAb-2,nAb);
j = 1;
for i=1:nAb
    if i ~= jM1 && i ~= jM2
        AmiTotals(j,:) = IdAb(i,:) -...
            IdAb(jM1,:)*(Bbs(i,1)*ibn11 + Bbs(i,2)*ibn21) -...
            IdAb(jM2,:)*(Bbs(i,1)*ibn12 + Bbs(i,2)*ibn22);
        j = j + 1;
    end
end
argumentos.AmiTotals = AmiTotals;
%{
%Teste:
ad = AmiTotals*Bbs;
%ok!
%}

%Difeomorfismo:
%Obs.: X = sigDif*Zb
sigDif0s = [beta0s; alfa0s; alfa1s];
[ AmiS ] = fAG( sigDif0s, AmiTotals, Abs );
sigDifs = [sigDif0s; AmiS];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FORMA NORMAL COM LQR (INCLUINDO SEPARAÇÃO ENTRE DINÂMICAS DE VIBRAÇÃO)
%Obs.: Zbb é o vetor original de estado

%Saídas:
C1s2 = [zeros(1,sPk(1,2)),Pk(2,:)]*Perlqv.';
C2s2 = [Pk(1,:),zeros(1,sPk(1,2))]*Perlqv.';

%Dinâmica externa (sistema com LQR):
%Obs.: yTil = aTils(X) + bTilb*Us
beta0s2 = C1s2;
beta1s2 = beta0s2*Abbs;
beta0Bbs2 = beta0s2*Bbbs;
%{
beta2 = beta1*Abs;
beta1Bbs = beta1*Bbs;
beta3 = beta2*Abs;
beta2Bbs = beta2*Bbs;
beta4 = beta3*Abs;
beta3Bbs = beta3*Bbs;
beta5 = beta4*Abs;
beta4Bbs = beta4*Bbs;
beta6 = beta5*Abs;
beta5Bbs = beta5*Bbs;
beta7 = beta6*Abs;
beta6Bbs = beta6*Bbs;
beta8 = beta7*Abs;
beta7Bbs = beta7*Bbs;
beta9 = beta8*Abs;
beta8Bbs = beta8*Bbs;
beta10 = beta9*Abs;
beta9Bbs = beta9*Bbs;
beta11 = beta10*Abs;
beta10Bbs = beta10*Bbs;
%}
alfa0s2 = C2s2;
alfa1s2 = alfa0s2*Abbs;
alfa2s2 = alfa1s2*Abbs;
alfa1Bbs2 = alfa1s2*Bbbs;
%{
alfa0Bbs = alfa0*Bbs;
alfa3 = alfa2*Abs;
alfa2Bbs = alfa2*Bbs;
alfa4 = alfa3*Abs;
alfa3Bbs = alfa3*Bbs;
ad1 = [beta0Bbs, beta1Bbs, beta2Bbs, beta3Bbs, beta4Bbs, beta5Bbs, beta6Bbs, beta7Bbs, beta8Bbs, beta9Bbs, beta10Bbs, alfa0Bbs, alfa1Bbs, alfa2Bbs, alfa3Bbs];
ad2 = [beta0Bbs, alfa0Bbs, alfa1Bbs];
%}

%Dinâmica interna (sistema com LQR):
%Obs.:
%(1) mi = Ami*Zb
%(2) mip = GamaMi*mi + SigmaMi*yTil + Ami*Bbs*Us
[maxBbb1,jM1] = max(abs(Bbbs(:,1)));
[maxBbb2,jM2] = max(abs(Bbbs(:,2)));
info.maxBbb1s = maxBbb1;
info.maxBbb2s = maxBbb2;
bn11 = Bbbs(jM1,1);
bn22 = Bbbs(jM2,2);
bn21 = Bbbs(jM2,1);
bn12 = Bbbs(jM1,2);
Bn = [bn11, bn21; bn12, bn22];
invBn = Bn\eye(2);
ibn11 = invBn(1,1);
ibn22 = invBn(2,2);
ibn21 = invBn(2,1);
ibn12 = invBn(1,2);
AmiTotals2 = zeros(nAb-2,nAb);
j = 1;
for i=1:nAb
    if i ~= jM1 && i ~= jM2
        AmiTotals2(j,:) = IdAb(i,:) -...
            IdAb(jM1,:)*(Bbbs(i,1)*ibn11 + Bbbs(i,2)*ibn21) -...
            IdAb(jM2,:)*(Bbbs(i,1)*ibn12 + Bbbs(i,2)*ibn22);
        j = j + 1;
    end
end
argumentos.AmiTotals2 = AmiTotals2;
%{
%Teste:
ad = AmiTotals2*Bbbs;
%ok!
%}

%Difeomorfismo:
%Obs.: X = sigDif*Zb
sigDif0s2 = [beta0s2; alfa0s2; alfa1s2];
[ AmiS2 ] = fAG( sigDif0s2, AmiTotals2, Abbs );
sigDifs2 = [sigDif0s2; AmiS2];

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
%beta5NM1 = beta4NM1*Ab;
%beta4NM1Bb = beta4NM1*Bb;
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
alfa1NM1 = alfa0NM1*Ab;
%alfa2NM1 = alfa1NM1*Ab;
%alfa1NM1Bb = alfa1NM1*Bb;
%{
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
%sigDif0NM1 = [beta0NM1; beta1NM1; beta2NM1; beta3NM1; beta4NM1; alfa0NM1; alfa1NM1];
%{
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
%
lmiNM1 = 13; %Linha retirada de AmiTotal, com "melhor" resultado geral
AmiNM1 = [AmiTotal(1:lmiNM1-1,:); AmiTotal(lmiNM1+1:nAb-2,:)];
%}
%sigDifNM1 = [sigDif0NM1; AmiNM1];
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