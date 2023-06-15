function [ argumentos, repositorio ] = zzintegracaonmk( argumentos, repositorio )
simulacao = argumentos.simulacao;
t0 = argumentos.t0;
dt = argumentos.dt_newmark;
T = argumentos.T;
nT = argumentos.nT;
fupnM1 = repositorio.fupnM1;
fuppnM1 = repositorio.fuppnM1;
fi1Nmk = repositorio.fi1Nmk; %fi1Nmk = @(dtnM1,un,upn,uppn)
fi2Nmk = repositorio.fi2Nmk; %fi2Nmk = @(dtnM1,un,upn,uppn)
fdupduTnmk = repositorio.fdupduTnmk; %fdupduTnmk = @(dtnM1)
%fduppduTnmk = repositorio.fduppduTnmk;
ind_LVC = argumentos.ind_LVC;

switch simulacao
    case 1
        M = argumentos.M;
        C = argumentos.C;
        K = argumentos.K;
        B = argumentos.B;
        W0 = argumentos.W0;
        fTilc = repositorio.fTilc; %fTilc = @(ubr,ubpr)
        %{
        fUc = repositorio.fUc; %fUc = @(t,u1bp)
        %}
        fUc = repositorio.fUc; %fUc = @(t)
        %Integração:
        fTilcEta = repositorio.fTilcEta;
        DfTilcEtaDEtaT = repositorio.DfTilcEtaDEtaT;
        Knmk = repositorio.Knmk;
        %Outros parâmetros e funções:
        Ng = argumentos.Ng;
        r = argumentos.rEstr;
        %{
        Idr = eye(Ng-r);
        %}
    case 2
        M = argumentos.M;
        C = argumentos.C;
        K = argumentos.K;
        B = argumentos.B;
        W0 = argumentos.W0;
        fTilc = repositorio.fTilc; %fTilc = @(ubrk,ubprk)
        %{
        fUc = repositorio.fUc; %fUc = @(t,u1bp)
        %}
        fUc = repositorio.fUc; %fUc = @(t)
        %Integração:
        fTilcEta = repositorio.fTilcEta;
        DfTilcEtaDEtaT = repositorio.DfTilcEtaDEtaT;
        Knmk = repositorio.Knmk;
        %Outros parâmetros e funções:
        Ng = argumentos.Ng;
        r = argumentos.rEstr;
        Pk = argumentos.Pk;
        sPk = size(Pk);
        npk = sPk(1,1);
        mpk = sPk(1,2);
        %{
        Idr = eye(Ng-r);
        %}
    case 3
        M = argumentos.M;
        C = argumentos.C;
        K = argumentos.K;
        B = argumentos.B;
        W0 = argumentos.W0;
        fTilc = repositorio.fTilc; %fTilc = @(ubr,ubpr,ubppr)
        fUc = repositorio.fUc;
        %Integração:
        fTilcEta = repositorio.fTilcEta;
        DfTilcEtaDEtaT = repositorio.DfTilcEtaDetaT;
        Knmk = repositorio.Knmk;
        %Outros parâmetros e funções:
        Ng = argumentos.Ng;
        r = argumentos.rEstr;
        %{
        Idr = eye(Ng-r);
        %}
    case 4
        M = argumentos.M;
        C = argumentos.C;
        K = argumentos.K;
        B = argumentos.B;
        W0 = argumentos.W0;
        fTilc = repositorio.fTilc; %fTilc = @(ubr,ubpr,ubppr)
        fUc = repositorio.fUc;
        %Integração:
        fTilcEta = repositorio.fTilcEta;
        DfTilcEtaDEtaT = repositorio.DfTilcEtaDEtaT;
        Bz = argumentos.Bz;
        Kz = repositorio.Kz;
        W0z = argumentos.W0z;
        %Outros parâmetros e funções:
        Ng = argumentos.Ng;
        r = argumentos.rEstr;
        %{
        IdNg = eye(Ng);
        %}
        Auz = argumentos.Auz;
        Almbz = argumentos.Almbz;
        fiNmkz = @(dtnM1,uTilbn,uTilbpn,uTilbppn)(...
            [M*fi2Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn) + C*fi1Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn); 
             zeros(r,1)]);
         kNmk = argumentos.k_newmark;
    case 5
        Afix = argumentos.Afix;
        Amix = argumentos.Amix;
        GamaFiBdxnM1 = repositorio.GamaFiBdxnM1; %GamaFiBdxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)
        KlcdxnM1 = repositorio.KlcdxnM1; %KlcdxnM1 = @(tnM1,xTilnM1)
        Kfi = argumentos.Kfi;
        GamaMis = argumentos.GamaMis;
        SigmaMis = argumentos.SigmaMis;
        fmidp = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(Amix.'*fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn));
        dGamaFidfipBdTxnM1 = repositorio.dGamaFidfipBdTxnM1; %dGamaFidfipBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)
        dGamaEtadEtaT = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            dGamaFidfipBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)*Afix.'*fdupduTnmk(dtnM1));
        dKlcddmidTxnM1 = repositorio.dKlcddmidTxnM1; %dKlcddmidTxnM1 = @(tnM1,xTilnM1)
        Kfix = Kfi*Afix.';
        dKlcdEtadEtaT = @(tnM1,eta)(dKlcddmidTxnM1(tnM1,eta)*Amix.');
        fdmidpdEtaT = @(dtnM1)(Amix.'*fdupduTnmk(dtnM1));
        %Outros parâmetros e expressões:
        Klcd = repositorio.Klcds; %Klcd = @(t,mid)
        %dKlcddmidT = repositorio.dKlcddmidT; %dKlcddmidT = @(t,mid)
        %DfiBd = repositorio.DfiBd; %DfiBd = @(fipBdTil)
        ydTil = repositorio.ydTils;
        %ypdTil = repositorio.ypdTil;
        Perl = argumentos.Perl;
        jPu2 = argumentos.jPu2;
        AmiNM1 = argumentos.AmiNM1;
        fZd = repositorio.fZd; %fZbd = @(t,mid)
        fZdx = @(t,xTil)(fZd(t,Amix.'*xTil));
        ffippBd = repositorio.ffippBd; %ffippBd = @(t,mid,mipd,fipBd)
        fmipd = repositorio.fmipd; %fmidp = @(t,mid)
        fmippd = repositorio.fmippd; %fmippd = @(t,mipd)
    case 6
        Afix = argumentos.Afix;
        Amix = argumentos.Amix;
        GamaFiBdxnM1 = repositorio.GamaFiBdxnM1; %GamaFiBdxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)
        KlcdxnM1 = repositorio.KlcdxnM1; %KlcdxnM1 = @(tnM1,xTilnM1)
        Kfi = argumentos.Kfi;
        GamaMi = argumentos.GamaMi;
        SigmaMi = argumentos.SigmaMi;
        fmidp = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(Amix.'*fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn));
        dGamaFidfipBdTxnM1 = repositorio.dGamaFidfipBdTxnM1; %dGamaFidfipBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)
        dGamaEtadEtaT = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            dGamaFidfipBdTxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)*Afix.'*fdupduTnmk(dtnM1));
        dKlcddmidTxnM1 = repositorio.dKlcddmidTxnM1; %dKlcddmidTxnM1 = @(tnM1,xTilnM1)
        Kfix = Kfi*Afix.';
        dKlcdEtadEtaT = @(tnM1,eta)(dKlcddmidTxnM1(tnM1,eta)*Amix.');
        fdmidpdEtaT = @(dtnM1)(Amix.'*fdupduTnmk(dtnM1));
        %Outros parâmetros e expressões:
        Klcd = repositorio.Klcd; %Klcd = @(t,mid)
        %dKlcddmidT = repositorio.dKlcddmidT; %dKlcddmidT = @(t,mid)
        %DfiBd = repositorio.DfiBd; %DfiBd = @(fipBdTil)
        ydTil = repositorio.ydTil;
        %ypdTil = repositorio.ypdTil;
        Perl = argumentos.Perl;
        jPu2 = argumentos.jPu2;
        AmiS = argumentos.AmiS;
        fZbd = repositorio.fZbd; %fZbd = @(t,mid)
        fZbdx = @(t,xTil)(fZbd(t,Amix.'*xTil));
        ffippBd = repositorio.ffippBdFN; %ffippBd = @(t,mid,mipd,fipBd)
        fmipd = repositorio.fmipd; %fmidp = @(t,mid)
        fmippd = repositorio.fmippd; %fmippd = @(t,mipd)
    case 7
        
    case 8
        
    case 9
        
    case 10
        
    case 11
        
    case 12
        
    otherwise
        disp('other value')
end

tn = t0;
switch simulacao
    case 1
        uTilbr0 = zeros(Ng-r,1);
        uTilbpr0 = zeros(Ng-r,1);
        %{
        u1bp0 = Idr(1,:)*uTilbpr0;
        Uc0 = fUc(tn,u1bp0);
        %}
        Uc0 = fUc(tn);
        uTilbppr0 = M\(B*Uc0 - C*uTilbpr0 - K*uTilbr0 - fTilc(uTilbr0,uTilbpr0));
        ZZ = zeros(nT,3*(Ng-r) + 2);
        ZZ(1,:) = [uTilbr0.', uTilbpr0.', uTilbppr0.', Uc0.'];
        uTilbn = uTilbr0;
        uTilbpn = uTilbpr0;
        uTilbppn = uTilbppr0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            uTilbnM10 = uTilbn + dtnM1*uTilbpn + 0.5*dtnM1^2*uTilbppn;
            %{
            %Velocidade longitudinal do nó 1:
            fu1bpnM1 = @(eta)(Idr(1,:)*fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn));
            %}
            %Integração:
            Eta0 = uTilbnM10;
            %{
            Psi = @(eta)(Knmk(dtnM1)*eta + fTilcEta(dtnM1,eta,uTilbn,uTilbpn,uTilbppn) - B*fUc(tnM1,fu1bpnM1(eta)) +...
                M*fi2Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn) + C*fi1Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn));
            %}
            Psi = @(eta)(Knmk(dtnM1)*eta + fTilcEta(dtnM1,eta,uTilbn,uTilbpn,uTilbppn) - B*fUc(tnM1) - W0 +...
                M*fi2Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn) + C*fi1Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn));
            dPsidEtaT = @(eta)(Knmk(dtnM1) + DfTilcEtaDEtaT(dtnM1,eta,uTilbn,uTilbpn,uTilbppn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            uTilbnM1 = Eta_conv;
            uTilbpnM1 = fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn);
            uTilbppnM1 = fuppnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn);
            %{
            UcnM1 = fUc(tnM1,fu1bpnM1(uTilbnM1));
            %}
            UcnM1 = fUc(tnM1);
            ZZ(i,:) = [uTilbnM1.', uTilbpnM1.', uTilbppnM1.', UcnM1.'];
            uTilbn = uTilbnM1;
            uTilbpn = uTilbpnM1;
            uTilbppn = uTilbppnM1;
        end
    case 2
        uTilbr0 = Pk.'*zeros(Ng-r,1);
        uTilbpr0 = Pk.'*zeros(Ng-r,1);
        %{
        u1bp0 = Idr(1,:)*(Pk*uTilbpr0);
        Uc0 = fUc(tn,u1bp0);
        %}
        Uc0 = fUc(tn);
        uTilbppr0 = M\(B*Uc0 + W0 - C*uTilbpr0 - K*uTilbr0 - fTilc(uTilbr0,uTilbpr0));
        ZZ = zeros(nT,3*(Ng-r) + 2)*...
            [Pk, zeros(npk,mpk), zeros(npk,mpk), zeros(npk,2); 
             zeros(npk,mpk), Pk, zeros(npk,mpk), zeros(npk,2); 
             zeros(npk,mpk), zeros(npk,mpk), Pk, zeros(npk,2); 
             zeros(2,mpk), zeros(2,mpk), zeros(2,mpk), eye(2)];
        ZZ(1,:) = [uTilbr0.', uTilbpr0.', uTilbppr0.', Uc0.'];
        uTilbn = uTilbr0;
        uTilbpn = uTilbpr0;
        uTilbppn = uTilbppr0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            uTilbnM10 = uTilbn + dtnM1*uTilbpn + 0.5*dtnM1^2*uTilbppn;
            %{
            %Velocidade longitudinal do nó 1:
            fu1bpnM1 = @(eta)(Idr(1,:)*(Pk*fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)));
            %}
            %Integração:
            Eta0 = uTilbnM10;
            %{
            Psi = @(eta)(Knmk(dtnM1)*eta + fTilcEta(dtnM1,eta,uTilbn,uTilbpn,uTilbppn) - B*fUc(tnM1,fu1bpnM1(eta)) +...
                M*fi2Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn) + C*fi1Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn));
            %}
            Psi = @(eta)(Knmk(dtnM1)*eta + fTilcEta(dtnM1,eta,uTilbn,uTilbpn,uTilbppn) - B*fUc(tnM1) - W0 +...
                M*fi2Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn) + C*fi1Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn));
            dPsidEtaT = @(eta)(Knmk(dtnM1) + DfTilcEtaDEtaT(dtnM1,eta,uTilbn,uTilbpn,uTilbppn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            uTilbnM1 = Eta_conv;
            uTilbpnM1 = fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn);
            uTilbppnM1 = fuppnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn);
            %{
            UcnM1 = fUc(tnM1,fu1bpnM1(uTilbnM1));
            %}
            UcnM1 = fUc(tnM1);
            ZZ(i,:) = [uTilbnM1.', uTilbpnM1.', uTilbppnM1.', UcnM1.'];
            uTilbn = uTilbnM1;
            uTilbpn = uTilbpnM1;
            uTilbppn = uTilbppnM1;
        end
    case 3
        uTilbr0 = zeros(Ng-r,1);
        uTilbpr0 = zeros(Ng-r,1);
        uTilbppr0 = zeros(Ng-r,1);
        %{
        u1bp0 = Idr(1,:)*uTilbpr0;
        Uc0 = fUc(tn,u1bp0);
        %}
        Uc0 = fUc(tn);
        uTilbppr0 = M\(B*Uc0 + W0 - C*uTilbpr0 - K*uTilbr0 - fTilc(uTilbr0,uTilbpr0,uTilbppr0));
        ZZ = zeros(nT,3*(Ng-r) + 2);
        ZZ(1,:) = [uTilbr0.', uTilbpr0.', uTilbppr0.', Uc0.'];
        uTilbn = uTilbr0;
        uTilbpn = uTilbpr0;
        uTilbppn = uTilbppr0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            uTilbnM10 = uTilbn + dtnM1*uTilbpn + 0.5*dtnM1^2*uTilbppn;
            %{
            %Velocidade longitudinal do nó 1:
            fu1bpnM1 = @(eta)(Idr(1,:)*fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn));
            %}
            %Integração:
            Eta0 = uTilbnM10;
            %{
            Psi = @(eta)(Knmk(dtnM1)*eta + fTilcEta(dtnM1,eta,uTilbn,uTilbpn,uTilbppn) - B*fUc(tnM1,fu1bpnM1(eta)) +...
                M*fi2Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn) + C*fi1Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn));
            %}
            Psi = @(eta)(Knmk(dtnM1)*eta + fTilcEta(dtnM1,eta,uTilbn,uTilbpn,uTilbppn) - B*fUc(tnM1) - W0 +...
                M*fi2Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn) + C*fi1Nmk(dtnM1,uTilbn,uTilbpn,uTilbppn));
            dPsidEtaT = @(eta)(Knmk(dtnM1) + DfTilcEtaDEtaT(dtnM1,eta,uTilbn,uTilbpn,uTilbppn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            uTilbnM1 = Eta_conv;
            uTilbpnM1 = fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn);
            uTilbppnM1 = fuppnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn);
            %{
            UcnM1 = fUc(tnM1,fu1bpnM1(uTilbnM1));
            %}
            UcnM1 = fUc(tnM1);
            ZZ(i,:) = [uTilbnM1.', uTilbpnM1.', uTilbppnM1.', UcnM1.'];
            uTilbn = uTilbnM1;
            uTilbpn = uTilbpnM1;
            uTilbppn = uTilbppnM1;
        end
    case 4
        uTilb0 = zeros(Ng,1);
        uTilbp0 = zeros(Ng,1);
        uTilbpp0 = zeros(Ng,1);
        %{
        u1bp0 = IdNg(1,:)*uTilbp0;
        Uc0 = fUc(tn,u1bp0);
        %}
        Uc0 = fUc(tn);
        uTilbpp0 = M\(B*Uc0 + W0 - C*uTilbp0 - K*uTilb0 - fTilc(uTilb0,uTilbp0,uTilbpp0));
        lambdaR0 = zeros(r,1);
        ZZ = zeros(nT,3*Ng + 2 + r);
        ZZ(1,:) = [uTilb0.', uTilbp0.', uTilbpp0.', Uc0.', lambdaR0.'];
        uTilbn = uTilb0;
        uTilbpn = uTilbp0;
        uTilbppn = uTilbpp0;
        lambdaRn = lambdaR0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            %Chute inicial:
            uTilbnM10 = uTilbn + dtnM1*uTilbpn + 0.5*dtnM1^2*uTilbppn;
            lambdaRnM10 = lambdaRn;
            %{
            %Velocidade longitudinal do nó 1:
            fu1bpnM1 = @(eta)(IdNg(1,:)*fupnM1(dtnM1,Auz.'*eta,uTilbn,uTilbpn,uTilbppn));
            %}
            %Integração:
            Eta0 = [uTilbnM10; lambdaRnM10];
            %{
            Psi = @(eta)(Kz(dtnM1)*eta + fTilcEta(dtnM1,eta,uTilbn,uTilbpn,uTilbppn) - Bz*fUc(tnM1,fu1bpnM1(eta)) +...
                fiNmkz(dtnM1,uTilbn,uTilbpn,uTilbppn));
            %}
            Psi = @(eta)(Kz(dtnM1)*eta + fTilcEta(dtnM1,eta,uTilbn,uTilbpn,uTilbppn) - Bz*fUc(tnM1) - W0z +...
                fiNmkz(dtnM1,uTilbn,uTilbpn,uTilbppn));
            dPsidEtaT = @(eta)(Kz(dtnM1) + DfTilcEtaDEtaT(dtnM1,eta,uTilbn,uTilbpn,uTilbppn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            %Distinção de saídas:
            uTilbnM1 = Auz.'*znM1;
            uTilbpnM1 = fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn);
            uTilbppnM1 = fuppnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn);
            %{
            UcnM1 = fUc(tnM1,fu1bpnM1(znM1));
            %}
            UcnM1 = fUc(tnM1);
            lambdaRnM1 = Almbz.'*znM1;
            %Armazenamento:
            ZZ(i,:) = [uTilbnM1.', uTilbpnM1.', uTilbppnM1.', UcnM1.', kNmk*lambdaRnM1.'];
            %Atualização:
            uTilbn = uTilbnM1;
            uTilbpn = uTilbpnM1;
            uTilbppn = uTilbppnM1;
            lambdaRn = lambdaRnM1;
        end
    case 5
        t0 = 0;
        ydTil0 = ydTil(t0);
        uTilb0 = zeros(jPu2,1);
        uTilbp0 = zeros(jPu2,1);
        Z0 = [uTilb0; uTilbp0];
        Zb0 = Perl*Z0;
        mid0 = AmiS*Zb0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBd0 = Kfi\Klcd(t0,mid0);
        fipBd0 = [0; 0];
        fippBd0 = ffippBd(t0,mid0,mipd0,fipBd0);
        xTil0 = [fiBd0; mid0];
        xTilp0 = [fipBd0; mipd0];
        xTilpp0 = [fippBd0; mippd0];
        Zbd0 = fZbdx(t0,xTil0);
        Zd0 = Perl*Zbd0;
        ZZ = zeros(nT,3*length(xTil0) + length(ydTil0) + length(Zd0));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', ydTil0.', Zd0.'];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            Eta0 = xTilnM10;
            Psi = @(eta)([GamaFiBdxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + Kfix*eta - KlcdxnM1(tnM1,eta);
                fmidp(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMis*Amix.'*eta - SigmaMis*ydTil(tnM1)]);
            dPsidEtaT = @(eta)([dGamaEtadEtaT(dtnM1,eta,xTiln,xTilpn,xTilppn) + Kfix - dKlcdEtadEtaT(tnM1,eta);
                fdmidpdEtaT(dtnM1) - GamaMis*Amix.']);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            fipBdnM1 = Afix.'*xTilpnM1;
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fippBdnM1 = ffippBd(tnM1,midnM1,mipdnM1,fipBdnM1);
            xTilppnM1 = [fippBdnM1; mippdnM1];
            %}
            %Trajetória desejada:
            ydTilnM1 = ydTil(tnM1);
            %Vetor de estados desejado:
            ZbdnM1 = fZbdx(tnM1,xTilnM1);
            ZdnM1 = Perl.'*ZbdnM1;
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', ydTilnM1.', ZdnM1.'];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 6
        t0 = 0;
        ydTil0 = ydTil(t0);
        uTilb0 = zeros(jPu2,1);
        uTilbp0 = zeros(jPu2,1);
        Z0 = [uTilb0; uTilbp0];
        mid0 = AmiNM1*Z0;
        mipd0 = fmipd(t0,mid0);
        mippd0 = fmippd(t0,mipd0);
        fiBd0 = Kfi\Klcd(t0,mid0);
        fipBd0 = [0; 0];
        fippBd0 = ffippBd(t0,mid0,mipd0,fipBd0);
        xTil0 = [fiBd0; mid0];
        xTilp0 = [fipBd0; mipd0];
        xTilpp0 = [fippBd0; mippd0];
        Zd0 = fZdx(t0,xTil0);
        ZZ = zeros(nT,3*length(xTil0) + length(ydTil0) + length(Zd0));
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', ydTil0.', Zd0.'];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            Eta0 = xTilnM10;
            Psi = @(eta)([GamaFiBdxnM1(dtnM1,eta,xTiln,xTilpn,xTilppn) + Kfix*eta - KlcdxnM1(tnM1,eta);
                fmidp(dtnM1,eta,xTiln,xTilpn,xTilppn) - GamaMi*Amix.'*eta - SigmaMi*ydTil(tnM1)]);
            dPsidEtaT = @(eta)([dGamaEtadEtaT(dtnM1,eta,xTiln,xTilpn,xTilppn) + Kfix - dKlcdEtadEtaT(tnM1,eta);
                fdmidpdEtaT(dtnM1) - GamaMi*Amix.']);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            %
            %Maior precisão na estimação de fippBdnM1:
            fipBdnM1 = Afix.'*xTilpnM1;
            midnM1 = Amix.'*xTilnM1;
            mipdnM1 = fmipd(tnM1,midnM1);
            mippdnM1 = fmippd(tnM1,mipdnM1);
            fippBdnM1 = ffippBd(tnM1,midnM1,mipdnM1,fipBdnM1);
            xTilppnM1 = [fippBdnM1; mippdnM1];
            %}
            %Trajetória desejada:
            ydTilnM1 = ydTil(tnM1);
            %Vetor de estados desejado:
            ZdnM1 = fZdx(tnM1,xTilnM1);
            %Armazenamento:
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', ydTilnM1.', ZdnM1.'];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
        end
    case 7
        uTilb0 = zeros(6,1);
        uTilc0 = fuTilc(uTilb0);
        uTilbp0 = zeros(6,1);
        uTilcp0 = fuTilc(uTilbp0);
        Z0 = [uTilc0; uTilcp0];
        mid0 = Ami*Z0;
        Zd0 = fcZdHat(t0,mid0,Z0);
        fiBd0 = Kfi^-1*Kcd(t0,mid0,Zd0);
        mipd0 = fcMidpHat(t0,mid0,Zd0);
        fipBd0 = 0;
        mippd0 = Ami*AbHat*fcZpdHat(t0,mid0,Zd0);
        fippBd0 = -Kfi*fipBd0 + DKcDXdT(t0,mid0,Zd0)*fcXpdHat(t0,mid0,Zd0);
        Uc0 = Uccu(t0,uTilc0,uTilcp0,fiBd0,mid0,Zd0);
        tetaP0 = tetaPuc(uTilc0);
        uTilbpp0 = zeros(6,1);
        uTilbpp0 = M1\(B6(uTilb0)*[Uc0; tetaP0] - C*uTilbp0 - K*uTilb0 - fTilc(uTilb0,uTilbp0,uTilbpp0));
        xTil0 = [uTilb0; fiBd0; mid0];
        xTilp0 = [uTilbp0; fipBd0; mipd0];
        xTilpp0 = [uTilbpp0; fippBd0; mippd0];
        X0 = fcXhatx(xTil0,xTilp0);
        Upc0 = Upccx(t0,xTil0,xTilp0,Zd0);
        ZZ = zeros(nT,3*(6+1+1)+4+1+1+1+4);
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', Zd0.', Uc0, Upc0, tetaP0, X0.'];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        Zdn = Zd0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            Eta0 = xTilnM10;
            Psi = @(eta)([KnmkEta(dtnM1)*eta + fTilcEta(dtnM1,eta,xTiln,xTilpn,xTilppn) - B6Eta(eta)*[UccEta(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,Zdn); tetaPEta(dtnM1,eta,xTiln,xTilpn,xTilppn)] + M1x*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                GamaFiBdEta(dtnM1,eta,xTiln,xTilpn,xTilppn) + Kfix*eta - KcdEta(tnM1,eta,Zdn);
                fmidpEta(dtnM1,eta,xTiln,xTilpn,xTilppn) - Ami*AbHat*fcZdHatEta(tnM1,eta,Zdn) - (Ami*Bhat)*UccMidEta(tnM1,eta,Zdn)]);
            dPsidEtaT = @(eta)([KnmkEta(dtnM1) + DfTilcEtaDetaT(dtnM1,eta,xTiln,xTilpn,xTilppn) - dB6UccTetaPEtadEtaT(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zdn);
                dGamaEtadEtaT(dtnM1,eta,xTiln,xTilpn,xTilppn) + Kfix - dKcEtadEtaT(tnM1,eta,Zdn);
                fdmidpdEtaT(dtnM1) - Ami*AbHat*dfcZdHatEtadEtaT(tnM1,eta,Zdn) - (Ami*Bhat)*dUccMidEtadEtaT(tnM1,eta,Zdn)]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            ZdnM1 = fcZdHatx(tnM1,xTilnM1,Zdn);
            XnM1 = fcXhatx(xTilnM1,xTilpnM1);
            UcnM1 = Uccx(tnM1,xTilnM1,xTilpnM1,Zdn);
            UpcnM1 = Upccx(tnM1,xTilnM1,xTilpnM1,Zdn);
            tetaPnM1 = tetaPx(xTilpnM1);
            if ind_LVC
                %Limitação das variáveis cíclicas teta2 e teta3:
                xTilnM1 = xTilnM1 - 2*pi*[1; zeros(4,1); 1; 0; q]*double(UcnM1>=2*pi);
                ZdnM1 = ZdnM1 - 2*pi*[ones(2,1); zeros(2,1)]*double(UcnM1>=2*pi);
                XnM1 = XnM1 - 2*pi*[zeros(3,1); q]*double(UcnM1>=2*pi);
                UcnM1 = UcnM1 - 2*pi*double(UcnM1>=2*pi);
            end
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', ZdnM1.', UcnM1, UpcnM1, tetaPnM1, XnM1.'];
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            Zdn = ZdnM1;
        end
    case 8
        uTilbr0 = zeros(Ng-2,1);
        Uc0 = fUc(tn);
        uTilb0 = fubVc(uTilbr0,Uc0);
        uTilbpr0 = zeros(Ng-2,1);
        Upc0 = fUpc(tn);
        uTilbp0 = fubVc(uTilbpr0,Upc0);
        uTilbppr0 = zeros(Ng-2,1);
        Uppc0 = fUppc(tn);
        uTilbpp0 = fubVc(uTilbppr0,Uppc0); %Válido apenas para a próxima expressão:
        lambdaTil0 = zeros(rEstr,1);
        uTilbpp0 = M1gEF\(-CgEF*uTilbp0 - KgEF*uTilb0 - fcEFb(uTilb0,uTilbp0,uTilbpp0) -...
            kNmk*GamaR.'*lambdaTil0);
        ZZ = zeros(nT,3*Ng + 8 + 3*2);
        ZZ(1,:) = [uTilb0.', uTilbp0.', uTilbpp0.', lambdaTil0.',Uc0.',Upc0.',Uppc0.'];
        uTilbrn = Aur.'*uTilb0;
        uTilbprn = Aur.'*uTilbp0;
        uTilbpprn = Aur.'*uTilbpp0;
        lambdaTiln = lambdaTil0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            if ind_LVC
                %Limitação das variáveis cíclicas:
                Ucn = fUc(tnM1-dtnM1) - 2*pi*[1;0]*floor([1,0]*fUc(tnM1-dtnM1)/(2*pi));
                UcnM1 = fUc(tnM1) - 2*pi*[1;0]*floor([1,0]*fUc(tnM1)/(2*pi));
                uTilbrn = uTilbrn - 2*pi*IdLVC*double(([1,0]*UcnM1)<([1,0]*Ucn));
            else
                %Não limitação das variáveis cíclicas:
                Ucn = fUc(tnM1-dtnM1);
                UcnM1 = fUc(tnM1);
            end
            Upcn = fUpc(tnM1-dtnM1);
            Uppcn = fUppc(tnM1-dtnM1);
            uTilbrnM10 = uTilbrn + dtnM1*uTilbprn + 0.5*dtnM1^2*uTilbpprn;
            lambdaTilnM10 = lambdaTiln;
            znM10 = [uTilbrnM10;lambdaTilnM10];
            Eta0 = znM10;
            Psi = @(eta)(Psic(dtnM1,eta,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn));
            dPsidEtaT = @(eta)(dPsicdznM1T(dtnM1,eta,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            uTilbrnM1 = Azu.'*znM1;
            lambdaTilnM1 = Azlmb.'*znM1;
            uTilbprnM1 = fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn);
            uTilbpprnM1 = fuppnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn);
            %Vetor posição completo, e suas derivadas temporais:
            uTilbnM1 = fubVc(uTilbrnM1,UcnM1);
            UpcnM1 = fUpc(tnM1);
            uTilbpnM1 = fubVc(uTilbprnM1,UpcnM1);
            UppcnM1 = fUppc(tnM1);
            uTilbppnM1 = fubVc(uTilbpprnM1,UppcnM1);
            %Armazenamento:
            ZZ(i,:) = [uTilbnM1.', uTilbpnM1.', uTilbppnM1.', lambdaTilnM1.', UcnM1.', UpcnM1.', UppcnM1.'];
            %Atualização:
            uTilbrn = uTilbrnM1;
            uTilbprn = uTilbprnM1;
            uTilbpprn = uTilbpprnM1;
            lambdaTiln = lambdaTilnM1;
        end
    case 9
        uTilbr0 = zeros(Ng-2,1);
        Uc0 = fUc(tn);
        uTilb0 = fubVc(uTilbr0,Uc0);
        uTilbpr0 = zeros(Ng-2,1);
        Upc0 = fUpc(tn);
        uTilbp0 = fubVc(uTilbpr0,Upc0);
        uTilbppr0 = zeros(Ng-2,1);
        Uppc0 = fUppc(tn);
        uTilbpp0 = fubVc(uTilbppr0,Uppc0); %Válido apenas para a próxima expressão:
        lambdaTil0 = zeros(rEstr,1);
        uTilbpp0 = M1gEF\(-CgEF*uTilbp0 - KgEF*uTilb0 - fcEFb(uTilb0,uTilbp0,uTilbpp0) -...
            kNmk*GamaR.'*lambdaTil0);
        ZZ = zeros(nT,3*Ng + 8 + 3*2);
        ZZ(1,:) = [uTilb0.', uTilbp0.', uTilbpp0.', lambdaTil0.',Uc0.',Upc0.',Uppc0.'];
        uTilbrn = Aur.'*uTilb0;
        uTilbprn = Aur.'*uTilbp0;
        uTilbpprn = Aur.'*uTilbpp0;
        lambdaTiln = lambdaTil0;
        Ucn = Uc0;
        Upcn = Upc0;
        Uppcn = Uppc0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            if ind_LVC
                %Limitação das variáveis cíclicas:
                UcnM1 = fUc(tnM1) - 2*pi*[1;0]*floor([1,0]*fUc(tnM1)/(2*pi));
                uTilbrn = uTilbrn - 2*pi*IdLVC*double(([1,0]*UcnM1)<([1,0]*Ucn));
            else
                %Não limitação das variáveis cíclicas:
                UcnM1 = fUc(tnM1);
            end
            uTilbrnM10 = uTilbrn + dtnM1*uTilbprn + 0.5*dtnM1^2*uTilbpprn;
            lambdaTilnM10 = lambdaTiln;
            znM10 = [uTilbrnM10;lambdaTilnM10];
            Eta0 = znM10;
            Psi = @(eta)(Psic(dtnM1,eta,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn));
            dPsidEtaT = @(eta)(dPsicdznM1T(dtnM1,eta,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            uTilbrnM1 = Azu.'*znM1;
            lambdaTilnM1 = Azlmb.'*znM1;
            uTilbprnM1 = fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn);
            uTilbpprnM1 = fuppnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn);
            %Vetor posição completo, e suas derivadas temporais:
            uTilbnM1 = fubVc(uTilbrnM1,UcnM1);
            %UpcnM1 = fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn);
            UpcnM1 = fUpc(tnM1);
            uTilbpnM1 = fubVc(uTilbprnM1,UpcnM1);
            UppcnM1 = fuppnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn);
            uTilbppnM1 = fubVc(uTilbpprnM1,UppcnM1);
            %Armazenamento:
            ZZ(i,:) = [uTilbnM1.', uTilbpnM1.', uTilbppnM1.', lambdaTilnM1.', UcnM1.', UpcnM1.', UppcnM1.'];
            %Atualização:
            uTilbrn = uTilbrnM1;
            uTilbprn = uTilbprnM1;
            uTilbpprn = uTilbpprnM1;
            lambdaTiln = lambdaTilnM1;
            Ucn = UcnM1;
            Upcn = UpcnM1;
            Uppcn = UppcnM1;
        end
    case 10
        uTilb0 = zeros(6,1);
        uTilc0 = fuTilc(uTilb0);
        uTilbp0 = zeros(6,1);
        uTilcp0 = fuTilc(uTilbp0);
        Z0 = [uTilc0; uTilcp0];
        mid0 = Ami*Z0;
        Zd0 = fcZdHat(t0,mid0,Z0);
        fiBd0 = Kfi^-1*Kcd(t0,mid0,Zd0);
        mipd0 = fcMidpHat(t0,mid0,Zd0);
        fipBd0 = 0;
        mippd0 = Ami*AbHat*fcZpdHat(t0,mid0,Zd0);
        %fippBd0 = -Kfi*fipBd0 + DKcDXdT(t0,mid0,Zd0)*fcXpdHat(t0,mid0,Zd0);
        Uc0 = Uccu(t0,uTilc0,uTilcp0,fiBd0,mid0,Zd0);
        tetaP0 = tetaPuc(uTilc0);
        uTilbpp0 = zeros(6,1);
        uTilbpp0 = M1\(B6(uTilb0)*[Uc0; tetaP0] - C*uTilbp0 - K*uTilb0 - fTilc(uTilb0,uTilbp0,uTilbpp0));
        xTil0 = [uTilb0; fiBd0; mid0];
        xTilp0 = [uTilbp0; fipBd0; mipd0];
        fippBd0 = fFippBdx(t0,xTil0,xTilp0,Zd0);
        xTilpp0 = [uTilbpp0; fippBd0; mippd0];
        X0 = fcXhatx(xTil0,xTilp0);
        Upc0 = Upccx(t0,xTil0,xTilp0,Zd0);
        ZZ = zeros(nT,3*(6+1+1)+4+1+1+1+4);
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', Zd0.', Uc0, Upc0, tetaP0, X0.'];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        Zdn = Zd0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            Eta0 = xTilnM10;
            Psi = @(eta)([KnmkEta(dtnM1)*eta + fTilcEta(dtnM1,eta,xTiln,xTilpn,xTilppn) - B6Eta(eta)*[UccEta(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,Zdn); tetaPEta(dtnM1,eta,xTiln,xTilpn,xTilppn)] + M1x*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn);
                GamaFiBdEta(dtnM1,eta,xTiln,xTilpn,xTilppn) + Kfix*eta - KcdEta(tnM1,eta,Zdn);
                fmidpEta(dtnM1,eta,xTiln,xTilpn,xTilppn) - Ami*AbHat*fcZdHatEta(tnM1,eta,Zdn) - (Ami*Bhat)*UccMidEta(tnM1,eta,Zdn)]);
            dPsidEtaT = @(eta)([KnmkEta(dtnM1) + DfTilcEtaDetaT(dtnM1,eta,xTiln,xTilpn,xTilppn) - dB6UccTetaPEtadEtaT(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zdn);
                dGamaEtadEtaT(dtnM1,eta,xTiln,xTilpn,xTilppn) + Kfix - dKcEtadEtaT(tnM1,eta,Zdn);
                fdmidpdEtaT(dtnM1) - Ami*AbHat*dfcZdHatEtadEtaT(tnM1,eta,Zdn) - (Ami*Bhat)*dUccMidEtadEtaT(tnM1,eta,Zdn)]);
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            ZdnM1 = fcZdHatx(tnM1,xTilnM1,Zdn);
            XnM1 = fcXhatx(xTilnM1,xTilpnM1);
            %Cálculo de fippBdnM1 por fórmula direta:
            fippBdnM1 = fFippBdx(tnM1,xTilnM1,xTilpnM1,Zdn);
            xTilppnM1(6+1,1) = fippBdnM1;
            %Entradas:
            UcnM1 = Uccx(tnM1,xTilnM1,xTilpnM1,Zdn);
            UpcnM1 = Upccx(tnM1,xTilnM1,xTilpnM1,Zdn);
            tetaPnM1 = tetaPx(xTilpnM1);
            if ind_LVC
                %Limitação das variáveis cíclicas teta2 e teta3:
                xTilnM1 = xTilnM1 - 2*pi*[1; zeros(4,1); 1; 0; q]*double(UcnM1>=2*pi);
                ZdnM1 = ZdnM1 - 2*pi*[ones(2,1); zeros(2,1)]*double(UcnM1>=2*pi);
                XnM1 = XnM1 - 2*pi*[zeros(3,1); q]*double(UcnM1>=2*pi);
                UcnM1 = UcnM1 - 2*pi*double(UcnM1>=2*pi);
            end
            ZZ(i,:) = [xTilnM1.', xTilpnM1.', xTilppnM1.', ZdnM1.', UcnM1, UpcnM1, tetaPnM1, XnM1.'];
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            Zdn = ZdnM1;
        end
    case 11
        uTilbr0 = zeros(Ng-2,1);
        uTilc0 = fuTilc(uTilbr0);
        uTilbpr0 = zeros(Ng-2,1);
        uTilcp0 = fuTilc(uTilbpr0);
        Z0 = [uTilc0; uTilcp0];
        mid0 = Ami*Z0;
        Zd0 = fcZdHat(t0,mid0,Z0);
        fiBd0 = Kfi^-1*Kcd(t0,mid0,Zd0);
        Uc0 = fUccr(t0,uTilbr0,uTilbpr0,fiBd0,mid0,Zd0);
        uTilb0 = fubVc(uTilbr0,Uc0);
        Upc0 = zeros(2,1);
        uTilbp0 = fubVc(uTilbpr0,Upc0);
        uTilbppr0 = zeros(Ng-2,1);
        Uppc0 = zeros(2,1);
        uTilbpp0 = fubVc(uTilbppr0,Uppc0); %Válido apenas para a próxima expressão:
        lambdaTil0 = zeros(rEstr,1);
        uTilbpp0 = M1gEF\(-CgEF*uTilbp0 - KgEF*uTilb0 - fcEFb(uTilb0,uTilbp0,uTilbpp0) -...
            kNmk*GamaR.'*lambdaTil0);
        mipd0 = fcMidpHat(t0,mid0,Zd0);
        fipBd0 = 0;
        mippd0 = Ami*AbHat*fcZpdHat(t0,mid0,Zd0);
        fFippBd = repositorio.fFippBd;
        fippBd0 = fFippBd(t0,fipBd0,mid0,Zd0);
        PsiDifcHat = repositorio.PsiDifcHat;
        X0 = PsiDifcHat(Z0);
        ZZ = zeros(nT,3*Ng+3*1+3*1+rEstr+2+2+2+4+4);
        ZZ(1,:) = [uTilb0.', uTilbp0.', uTilbpp0.', fiBd0, fipBd0, fippBd0, mid0, mipd0, mippd0, lambdaTil0.', Uc0.', Upc0.', Uppc0.', Zd0.', X0.'];
        xTil0 = [uTilbr0; fiBd0; mid0];
        xTilp0 = [uTilbpr0; fipBd0; mipd0];
        xTilpp0 = [uTilbppr0; fippBd0; mippd0];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        lambdaTiln = lambdaTil0;
        Ucn = Uc0;
        Upcn = Upc0;
        Uppcn = Uppc0;
        Zdn = Zd0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*(dtnM1*xTilppn)*dtnM1;
            lambdaTilnM10 = lambdaTiln;
            znM10 = [xTilnM10;lambdaTilnM10];
            Eta0 = znM10;
            Psi = @(eta)(Psic(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn));
            dPsidEtaT = @(eta)(dPsicdznM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            znM1 = Eta_conv;
            xTilnM1 = Azx.'*znM1;
            lambdaTilnM1 = Azlmb.'*znM1;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            ZdnM1 = fcZdHatx(tnM1,xTilnM1,Zdn);
            XnM1 = fcXhatx(xTilnM1,xTilpnM1);
            %Camada limite dinâmica:
            fiBdnM1 = AxFiBd.'*xTilnM1;
            fipBdnM1 = AxFiBd.'*xTilpnM1;
            %Cálculo de fippBdnM1 por fórmula direta:
            fippBdnM1 = fFippBdx(tnM1,xTilnM1,xTilpnM1,Zdn);
            %Vetor posição completo, e suas derivadas temporais:
            UcnM1 = fUccxnM1(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn);
            uTilrnM1 = Axur.'*xTilnM1;
            uTilbnM1 = fubVc(uTilrnM1,UcnM1);
            %UpcnM1 = fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn);
            xTilppnM1(Ng-2+1,1) = fippBdnM1;
            UpcnM1 = fUpccx(tnM1,xTilnM1,xTilpnM1,xTilppnM1,Zdn);
            uTilprnM1 = Axur.'*xTilpnM1;
            uTilbpnM1 = fubVc(uTilprnM1,UpcnM1);
            UppcnM1 = fuppnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn);
            uTilpprnM1 = Axur.'*xTilppnM1;
            uTilbppnM1 = fubVc(uTilpprnM1,UppcnM1);
            %Dinâmica interna:
            midnM1 = AxMi.'*xTilnM1;
            mipdnM1 = AxMi.'*xTilpnM1;
            mippdnM1 = AxMi.'*xTilppnM1;
            if ind_LVC
                %Limitação das variáveis cíclicas teta2 e teta3:
                xTilnM1 = xTilnM1 - 2*pi*[1; zeros(4,1); 1; 0; q]*double([1,0]*UcnM1>=2*pi);
                ZdnM1 = ZdnM1 - 2*pi*[ones(2,1); zeros(2,1)]*double([1,0]*UcnM1>=2*pi);
                XnM1 = XnM1 - 2*pi*[zeros(3,1); q]*double([1,0]*UcnM1>=2*pi);
                UcnM1 = UcnM1 - 2*pi*double([1,0]*UcnM1>=2*pi);
            end
            %Armazenamento:
            ZZ(i,:) = [uTilbnM1.', uTilbpnM1.', uTilbppnM1.', fiBdnM1, fipBdnM1, fippBdnM1, midnM1, mipdnM1, mippdnM1, lambdaTilnM1.', UcnM1.', UpcnM1.', UppcnM1.', ZdnM1.', XnM1.'];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            lambdaTiln = lambdaTilnM1;
            Ucn = UcnM1;
            Upcn = UpcnM1;
            Uppcn = UppcnM1;
            Zdn = ZdnM1;
        end
    case 12
        Z0 = zeros(4,1);
        mid0 = Ami*Z0;
        Zd0 = fcZdHat(t0,mid0,Z0);
        mipd0 = fcMidpHat(tn,mid0,Z0);
        mippd0 = Ami*AbHat*fcZpdHat(t0,mid0,Z0);
        ZZ = zeros(nT,3*1 + 4);        
        xTil0 = mid0; 
        xTilp0 = mipd0; 
        xTilpp0 = mippd0;
        ZZ(1,:) = [xTil0.', xTilp0.', xTilpp0.', Zd0.'];
        xTiln = xTil0;
        xTilpn = xTilp0;
        xTilppn = xTilpp0;
        Zdn = Zd0;
        for i = 2:nT
            tnM1 = T(i,1)
            dtnM1 = dt;
            xTilnM10 = xTiln + dtnM1*xTilpn + 0.5*dtnM1^2*xTilppn;
            Eta0 = xTilnM10;
            Psi = @(eta)(fmidpEta(dtnM1,eta,xTiln,xTilpn,xTilppn) -...
                Ami*AbHat*fcZdHatEta(tnM1,eta,Zdn) -...
                (Ami*Bhat)*UccMidEta(tnM1,eta,Zdn));
            dPsidEtaT = @(eta)(fdmidpdEtaT(dtnM1) -...
                Ami*AbHat*dfcZdHatEtadEtaT(tnM1,eta,Zdn) -...
                (Ami*Bhat)*dUccMidEtadEtaT(tnM1,eta,Zdn));
            [ Eta_conv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
            xTilnM1 = Eta_conv;
            xTilpnM1 = fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            xTilppnM1 = fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn);
            ZdnM1 = fcZdHatx(tnM1,xTilnM1,Z0);
            if ind_LVC
                %Limitação das variáveis cíclicas teta1, teta2 e teta3:
                xTilnM1 = xTilnM1 - 2*pi*[0; q]*double([1,0,0,0]*ZdnM1 >= 2*pi);
                ZdnM1 = ZdnM1 - 2*pi*[ones(2,1); zeros(2,1)]*double([1,0,0,0]*ZdnM1 >= 2*pi);
            end
            %Registro:
            ZZ(i,:) = [ xTilnM1.', xTilpnM1.', xTilppnM1.', ZdnM1.'];
            %Atualização:
            xTiln = xTilnM1;
            xTilpn = xTilpnM1;
            xTilppn = xTilppnM1;
            Zdn = ZdnM1;
        end
    otherwise
        disp('other value')
end
argumentos.ZZ = ZZ;
end