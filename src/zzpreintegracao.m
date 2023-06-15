function [ argumentos, repositorio ] = zzpreintegracao( argumentos, repositorio )
simulacao = argumentos.simulacao;
gamaNmk = argumentos.gama_newmark;
betaNmk = argumentos.beta_newmark;
fi1Nmk = @(dtnM1,un,upn,uppn)(-gamaNmk*(betaNmk*dtnM1)^-1*un +...
    (1 - gamaNmk*betaNmk^-1)*upn +...
    (1 - gamaNmk*(2*betaNmk)^-1)*dtnM1*uppn);
repositorio.fi1Nmk = fi1Nmk;
fi2Nmk = @(dtnM1,un,upn,uppn)(-(betaNmk*dtnM1^2)^-1*un -...
    (betaNmk*dtnM1)^-1*upn +...
    (1 - (2*betaNmk)^-1)*uppn);
repositorio.fi2Nmk = fi2Nmk;
fupnM1 = @(dtnM1,unM1,un,upn,uppn)(gamaNmk*(betaNmk*dtnM1)^-1*unM1 +...
    fi1Nmk(dtnM1,un,upn,uppn));
repositorio.fupnM1 = fupnM1;
fuppnM1 = @(dtnM1,unM1,un,upn,uppn)((betaNmk*dtnM1^2)^-1*unM1 +...
    fi2Nmk(dtnM1,un,upn,uppn));
repositorio.fuppnM1 = fuppnM1;
fdupduTnmk = @(dtnM1)(gamaNmk*(betaNmk*dtnM1)^-1);
repositorio.fdupduTnmk = fdupduTnmk;
fduppduTnmk = @(dtnM1)((betaNmk*dtnM1^2)^-1);
repositorio.fduppduTnmk = fduppduTnmk;
switch simulacao
    case 1
        %Matriz Knmk:
        M = argumentos.M;
        C = argumentos.C;
        K = argumentos.K;
        Knmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*M + gamaNmk*(betaNmk*dtnM1)^-1*C + K);
        repositorio.Knmk = Knmk;
        %Vetor fTilcEta:
        fTilc = repositorio.fTilc; %fTilc = @(ubr,ubpr)
        fTilcEta = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            fTilc(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)));
        repositorio.fTilcEta = fTilcEta;
        %Matriz dfTilc0EtadEtaT:
        dfTilcdubT = repositorio.dfTilcdubT;
        dfTilcdubpT = repositorio.dfTilcdubpT;
        duTilpnM1duTilnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        DfTilcEtaDEtaT = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            dfTilcdubT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)) +...
            dfTilcdubpT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn))*duTilpnM1duTilnM1(dtnM1));
        repositorio.DfTilcEtaDEtaT = DfTilcEtaDEtaT;
        %Matriz dfTildc0EtadEtaT:
        dfTildcdubT = repositorio.dfTildcdubT;
        dfTildcdubpT = repositorio.dfTildcdubpT;
        duTilpnM1duTilnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        DfTildcEtaDEtaT = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            dfTildcdubT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)) +...
            dfTildcdubpT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn))*duTilpnM1duTilnM1(dtnM1));
        repositorio.DfTildcEtaDEtaT = DfTildcEtaDEtaT;
    case 2
        %Matriz Knmk:
        M = argumentos.M;
        C = argumentos.C;
        K = argumentos.K;
        Knmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*M + gamaNmk*(betaNmk*dtnM1)^-1*C + K);
        repositorio.Knmk = Knmk;
        %Vetor fTilcEta:
        fTilc = repositorio.fTilc; %fTilc = @(ubr,ubpr)
        fTilcEta = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            fTilc(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)));
        repositorio.fTilcEta = fTilcEta;
        %Matriz dfTilc0EtadEtaT:
        dfTilcduT = repositorio.dfTilcduT;
        dfTilcdupT = repositorio.dfTilcdupT;
        duTilpnM1duTilnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        DfTilcEtaDEtaT = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            dfTilcduT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)) +...
            dfTilcdupT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn))*duTilpnM1duTilnM1(dtnM1));
        repositorio.DfTilcEtaDEtaT = DfTilcEtaDEtaT;
        %Matriz dfTildc0EtadEtaT:
        dfTildcduT = repositorio.dfTildcduT;
        dfTildcdupT = repositorio.dfTildcdupT;
        duTilpnM1duTilnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        DfTildcEtaDEtaT = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            dfTildcduT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)) +...
            dfTildcdupT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn))*duTilpnM1duTilnM1(dtnM1));
        repositorio.DfTildcEtaDEtaT = DfTildcEtaDEtaT;
    case 3
        %Matriz Knmk:
        M = argumentos.M;
        C = argumentos.C;
        K = argumentos.K;
        Knmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*M + gamaNmk*(betaNmk*dtnM1)^-1*C + K);
        repositorio.Knmk = Knmk;
        %Vetor fTilEta:
        fTilc = repositorio.fTilc; %fTilc = @(ubr,ubpr,ubppr)
        fTilcEta = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            fTilc(eta,...
            fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn),...
            fuppnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)));
        repositorio.fTilcEta = fTilcEta;
        %Matriz DfTilcEtaDEtaT:
        duTilpnM1duTilnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        duTilppnM1duTilnM1 = @(dtnM1)(fduppduTnmk(dtnM1));
        dfTilcdubT = repositorio.dfTilcdubT; %dfTilcduT = @(ubr,ubpr,ubppr)
        dfTilcdubpT = repositorio.dfTilcdubpT; %dfTilcdupT = @(ubr,ubpr)
        dfTildubppT = repositorio.dfTildubppT; %dfTilcduppT = @(ubr)
        DfTilcEtaDetaT = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            dfTilcdubT(eta,...
            fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn),...
            fuppnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)) +...
            dfTilcdubpT(eta,fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn))*...
            duTilpnM1duTilnM1(dtnM1) +...
            dfTildubppT(eta)*duTilppnM1duTilnM1(dtnM1));
        repositorio.DfTilcEtaDetaT = DfTilcEtaDetaT;
        %Matriz DfTildcEtaDEtaT:
        dfTildcdubT = repositorio.dfTildcdubT;
        dfTildcdubpT = repositorio.dfTildcdubpT;
        DfTildcEtaDetaT = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            dfTildcdubT(eta,...
            fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn),...
            fuppnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)) +...
            dfTildcdubpT(eta,...
            fupnM1(dtnM1,eta,uTilbn,uTilbpn,uTilbppn))*...
            duTilpnM1duTilnM1(dtnM1) +...
            dfTildubppT(eta)*duTilppnM1duTilnM1(dtnM1));
        repositorio.DfTildcEtaDetaT = DfTildcEtaDetaT;
    case 4
        Ng = argumentos.Ng;
        r = argumentos.rEstr;
        %Matriz Knmk:
        M = argumentos.M;
        C = argumentos.C;
        K = argumentos.K;
        B = argumentos.B;
        W0 = argumentos.W0;
        Knmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*M + gamaNmk*(betaNmk*dtnM1)^-1*C + K);
        repositorio.Knmk = Knmk;
        %Vetor fTilEta:
        fTilc = repositorio.fTilc; %fTilc = @(ubr,ubpr,ubppr)
        fTilcuTilb = @(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn)(...
            fTilc(uTilbnM1,...
            fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn),...
            fuppnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn)));
        repositorio.fTilcuTilb = fTilcuTilb;
        %Matriz DfTilcEtaDEtaT:
        duTilpnM1duTilnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        duTilppnM1duTilnM1 = @(dtnM1)(fduppduTnmk(dtnM1));
        dfTilcdubT = repositorio.dfTilcdubT; %dfTilcduT = @(ubr,ubpr,ubppr)
        dfTilcdubpT = repositorio.dfTilcdubpT; %dfTilcdupT = @(ubr,ubpr)
        dfTildubppT = repositorio.dfTildubppT; %dfTilcduppT = @(ubr)
        DfTilcDuTilbnM1T = @(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn)(...
            dfTilcdubT(uTilbnM1,...
            fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn),...
            fuppnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn)) +...
            dfTilcdubpT(uTilbnM1,fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn))*...
            duTilpnM1duTilnM1(dtnM1) +...
            dfTildubppT(uTilbnM1)*duTilppnM1duTilnM1(dtnM1));
        repositorio.DfTilcDuTilbnM1T = DfTilcDuTilbnM1T;
        %Matriz DfTildcEtaDEtaT:
        dfTildcdubT = repositorio.dfTildcdubT;
        dfTildcdubpT = repositorio.dfTildcdubpT;
        DfTildcDuTilbnM1T = @(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn)(...
            dfTildcdubT(uTilbnM1,...
            fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn),...
            fuppnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn)) +...
            dfTildcdubpT(uTilbnM1,fupnM1(dtnM1,uTilbnM1,uTilbn,uTilbpn,uTilbppn))*...
            duTilpnM1duTilnM1(dtnM1) +...
            dfTildubppT(uTilbnM1)*duTilppnM1duTilnM1(dtnM1));
        repositorio.DfTildcDuTilbnM1T = DfTildcDuTilbnM1T;
        %Sistema em z = [ubT, lambdaT]T:
        kNmk = argumentos.k_newmark;
        IdNgMr = eye(Ng+r);
        Auz = IdNgMr(:,1:Ng);
        argumentos.Auz = Auz;
        Almbz = IdNgMr(:,Ng+1:Ng+r);
        argumentos.Almbz = Almbz;
        GamaR = argumentos.GamaR;
        Kz = @(dtnM1)([Knmk(dtnM1), kNmk*GamaR.'; kNmk*GamaR, zeros(r,r)]);
        repositorio.Kz = Kz;
        fTilcEta = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            [fTilcuTilb(dtnM1,Auz.'*eta,uTilbn,uTilbpn,uTilbppn); 
             zeros(r,1)]);
        repositorio.fTilcEta = fTilcEta;
        Bz = [B; zeros(r,2)];
        argumentos.Bz = Bz;
        W0z = [W0; zeros(r,1)];
        argumentos.W0z = W0z;
        DfTilcEtaDEtaT = @(dtnM1,eta,uTilbn,uTilbpn,uTilbppn)(...
            [DfTilcDuTilbnM1T(dtnM1,Auz.'*eta,uTilbn,uTilbpn,uTilbppn), zeros(Ng,r); 
             zeros(r,Ng), zeros(r,r)]);
        repositorio.DfTilcEtaDEtaT = DfTilcEtaDEtaT;
    case 5
        nMid = argumentos.nMid;
        Idx = eye(nMid+2);
        Afix = Idx(:,1:2);
        argumentos.Afix = Afix;
        Amix = Idx(:,2+1:2+nMid);
        argumentos.Amix = Amix;
        GamaFiBd = repositorio.GamaFiBd; %GamaFiBd = @(fipBdTil)
        Klcd = repositorio.Klcd; %Klcd = @(t,mid)
        dGamaFidfipBdT = repositorio.dGamaFidfipBdT; %dGamaFidfipBdT = @(fipBdTil)
        dKlcddmidT = repositorio.dKlcddmidT; %dKlcddmidT = @(t,mid)
        
        GamaFiBdx = @(xTilp)(GamaFiBd(Afix.'*xTilp));
        Klcdx = @(t,xTil)(Klcd(t,Amix.'*xTil));
        dGamaFidfipBdTx = @(xTilp)(dGamaFidfipBdT(Afix.'*xTilp));
        dKlcddmidTx = @(t,xTil)(dKlcddmidT(t,Amix.'*xTil));
        
        GamaFiBdxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        repositorio.GamaFiBdxnM1 = GamaFiBdxnM1;
        KlcdxnM1 = @(tnM1,xTilnM1)(Klcdx(tnM1,xTilnM1));
        repositorio.KlcdxnM1 = KlcdxnM1;
        dGamaFidfipBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFidfipBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        repositorio.dGamaFidfipBdTxnM1 = dGamaFidfipBdTxnM1;
        dKlcddmidTxnM1 = @(tnM1,xTilnM1)(dKlcddmidTx(tnM1,xTilnM1));
        repositorio.dKlcddmidTxnM1 = dKlcddmidTxnM1;
    case 6
        nMid = argumentos.nMid;
        Idx = eye(nMid+2);
        Afix = Idx(:,1:2);
        argumentos.Afix = Afix;
        Amix = Idx(:,2+1:2+nMid);
        argumentos.Amix = Amix;
        GamaFiBd = repositorio.GamaFiBdFN; %GamaFiBd = @(fipBdTil)
        Klcd = repositorio.Klcd; %Klcd = @(t,mid)
        dGamaFidfipBdT = repositorio.dGamaFiFNdfipBdT; %dGamaFidfipBdT = @(fipBdTil)
        dKlcddmidT = repositorio.dKlcddmidT; %dKlcddmidT = @(t,mid)
        
        GamaFiBdx = @(xTilp)(GamaFiBd(Afix.'*xTilp));
        Klcdx = @(t,xTil)(Klcd(t,Amix.'*xTil));
        dGamaFidfipBdTx = @(xTilp)(dGamaFidfipBdT(Afix.'*xTilp));
        dKlcddmidTx = @(t,xTil)(dKlcddmidT(t,Amix.'*xTil));
        
        GamaFiBdxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            GamaFiBdx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        repositorio.GamaFiBdxnM1 = GamaFiBdxnM1;
        KlcdxnM1 = @(tnM1,xTilnM1)(Klcdx(tnM1,xTilnM1));
        repositorio.KlcdxnM1 = KlcdxnM1;
        dGamaFidfipBdTxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamaFidfipBdTx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        repositorio.dGamaFidfipBdTxnM1 = dGamaFidfipBdTxnM1;
        dKlcddmidTxnM1 = @(tnM1,xTilnM1)(dKlcddmidTx(tnM1,xTilnM1));
        repositorio.dKlcddmidTxnM1 = dKlcddmidTxnM1;
    case 7
        Id8 = eye(8);
        Axu = Id8(:,1:6);
        argumentos.Axu = Axu;
        AxFiBd = Id8(:,7);
        argumentos.AxFiBd = AxFiBd;
        AxMi = Id8(:,8);
        argumentos.AxMi = AxMi;
        
        %PLANTA:
        %Matriz Knmk:
        M1 = argumentos.M1;
        C = argumentos.C;
        K = argumentos.K;
        Knmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*M1 + gamaNmk*(betaNmk*dtnM1)^-1*C + K);
        repositorio.Knmk = Knmk;
        KnmkEta = @(dtnM1)(Knmk(dtnM1)*Axu.');
        repositorio.KnmkEta = KnmkEta;
        M1x = M1*Axu.';
        argumentos.M1x = M1x;
        Cx = C*Axu.';
        argumentos.Cx = Cx;
        %Vetor fTilEta:
        fTilc = repositorio.fTilc;
        fTilcx = @(xTil,xTilp,xTilpp)(fTilc(Axu.'*xTil,Axu.'*xTilp,Axu.'*xTilpp));
        fTilcEta = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            fTilcx(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),...
            fuppnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)));
        repositorio.fTilcEta = fTilcEta;
        %Vetor B6Eta:
        B6 = repositorio.B6;
        B6x = @(xTil)(B6(Axu.'*xTil));
        B6Eta = @(eta)(B6x(eta));
        repositorio.B6Eta = B6Eta;
        %Matriz DfTilcEtaDEtaT:
        duTilbdxTilT = Axu.';
        duTilbpdxTilpT = Axu.';
        duTilbppdxTilppT = Axu.';
        dxTilpnM1dxTilnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        dxTilppnM1dxTilnM1 = @(dtnM1)(fduppduTnmk(dtnM1));
        dfTilcduTilbT = repositorio.dfTilcdubT;
        dfTilcduTilbpT = repositorio.dfTilcdubpT;
        dfTilduTilbppT = repositorio.dfTildubppT;
        dfTilcxdxTilT = @(xTil,xTilp,xTilpp)(dfTilcduTilbT(Axu.'*xTil,Axu.'*xTilp,Axu.'*xTilpp)*duTilbdxTilT);
        dfTilcxdxTilpT = @(xTil,xTilp,xTilpp)(dfTilcduTilbpT(Axu.'*xTil,Axu.'*xTilp)*duTilbpdxTilpT);
        dfTilxdxTilppT = @(xTil,xTilp,xTilpp)(dfTilduTilbppT(Axu.'*xTil)*duTilbppdxTilppT);
        DfTilcEtaDetaT = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            dfTilcxdxTilT(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),...
            fuppnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)) +...
            dfTilcxdxTilpT(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn))*...
            dxTilpnM1dxTilnM1(dtnM1) +...
            dfTilxdxTilppT(eta)*dxTilppnM1dxTilnM1(dtnM1));
        repositorio.DfTilcEtaDetaT = DfTilcEtaDetaT;
        %Matriz DfTildcEtaDEtaT:
        dfTildcduTilbT = repositorio.dfTildcdubT;
        dfTildcduTilbpT = repositorio.dfTildcdubpT;
        dfTildcxdxTilT = @(xTil,xTilp,xTilpp)(dfTildcduTilbT(Axu.'*xTil,Axu.'*xTilp,Axu.'*xTilpp)*duTilbdxTilT);
        dfTildcxdxTilpT = @(xTil,xTilp,xTilpp)(dfTildcduTilbpT(Axu.'*xTil,Axu.'*xTilp)*duTilbpdxTilpT);
        DfTildcEtaDetaT = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            dfTildcxdxTilT(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),...
            fuppnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)) +...
            dfTildcxdxTilpT(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn))*...
            dxTilpnM1dxTilnM1(dtnM1) +...
            dfTilxdxTilppT(eta)*dxTilppnM1dxTilnM1(dtnM1));
        repositorio.DfTildcEtaDetaT = DfTildcEtaDetaT;
        
        %LEI DE CONTROLE:
        %Variáveis do vetor posição (lei de controle): uTilc = [teta2; teta3]
        fuTilc = repositorio.fuTilc;
        fuTilcx = @(xTil)(fuTilc(Axu.'*xTil));
        %Vetores UccEta e UcdcEta:
        Uccu = repositorio.Uccu;
        Uccx = @(t,xTil,xTilp,Zd0)(Uccu(t,fuTilcx(xTil),fuTilcx(xTilp),AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        repositorio.Uccx = Uccx;
        UccEta = @(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            Uccx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        repositorio.UccEta = UccEta;
        Ucdcu = repositorio.Ucdcu;
        Ucdcx = @(t,xTil,xTilp,Zd0)(Ucdcu(t,fuTilcx(xTil),fuTilcx(xTilp),AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        repositorio.Ucdcx = Ucdcx;
        UcdcEta = @(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            Ucdcx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        repositorio.UcdcEta = UcdcEta;
        %Vetor Upccx e Updccx:
        Upcc = repositorio.fUpcc; 
        Upccx = @(t,xTil,xTilp,Z0)(Upcc(t,fuTilcx(xTil),AxFiBd.'*xTil,AxMi.'*xTil,fuTilcx(xTilp),AxFiBd.'*xTilp,AxMi.'*xTilp,Z0));
        repositorio.Upccx = Upccx;
        Upcdc = repositorio.fUpcdc;
        Upcdcx = @(t,xTil,xTilp,Z0)(Upcdc(t,fuTilcx(xTil),AxFiBd.'*xTil,AxMi.'*xTil,fuTilcx(xTilp),AxFiBd.'*xTilp,AxMi.'*xTilp,Z0));
        repositorio.Upcdcx = Upcdcx;
        
        %Função tetaPEta:
        tetaPub = repositorio.tetaPub;
        tetaPx = @(xTilp)(tetaPub(Axu.'*xTilp));
        repositorio.tetaPx = tetaPx;
        tetaPEta = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(tetaPx(fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)));
        repositorio.tetaPEta = tetaPEta;
        
        %Vetor GamaFiEta:
        GamaFiBd = repositorio.GamaFiBd;
        GamaFiBdx = @(xTilp)(GamaFiBd(AxFiBd.'*xTilp));
        GamaFiBdEta = @(dtnM1,eta,xTilbn,xTilbpn,xTilbppn)(GamaFiBdx(fupnM1(dtnM1,eta,xTilbn,xTilbpn,xTilbppn)));
        repositorio.GamaFiBdEta = GamaFiBdEta;
        %Matriz Kfix:
        Kfi = argumentos.Kfi;
        Kfix = Kfi*AxFiBd.';
        argumentos.Kfix = Kfix;
        %Vetores KcdEta e KdcdEta:
        Kcd = repositorio.Kcd;
        Kcdx = @(t,xTil,Zd0)(Kcd(t,AxMi.'*xTil,Zd0));
        KcdEta = @(tnM1,eta,Zd0)(Kcdx(tnM1,eta,Zd0));
        repositorio.KcdEta = KcdEta;
        Kdcd = repositorio.Kdcd;
        Kdcdx = @(t,xTil,Z0)(Kdcd(t,AxMi.'*xTil,Z0));
        KdcdEta = @(tnM1,eta,Z0)(Kdcdx(tnM1,eta,Z0));
        repositorio.KdcdEta = KdcdEta;
        %Vetores fcZdHatEta e fdcZdHatEta:
        fcZdHat = repositorio.fcZdHat;
        fcZdHatx = @(t,xTil,Zd0)(fcZdHat(t,AxMi.'*xTil,Zd0));
        repositorio.fcZdHatx = fcZdHatx;
        fcZdHatEta = @(tnM1,eta,Zd0)(fcZdHatx(tnM1,eta,Zd0));
        repositorio.fcZdHatEta = fcZdHatEta;
        fdcZdHat = repositorio.fdcZdHat;
        fdcZdHatx = @(t,xTil,Z0)(fdcZdHat(t,AxMi.'*xTil,Z0));
        fdcZdHatEta = @(tnM1,eta,Z0)(fdcZdHatx(tnM1,eta,Z0));
        repositorio.fdcZdHatEta = fdcZdHatEta;
        %Expressões Uccd e Ucdcd:
        UccMid = repositorio.UccMid;
        UccMidx = @(t,xTil,Z0)(UccMid(t,AxMi.'*xTil,Z0));
        UccMidEta = @(tnM1,eta,Zd0)(UccMidx(tnM1,eta,Zd0));
        repositorio.UccMidEta = UccMidEta;
        UcdcMid = repositorio.UcdcMid;
        UcdcMidx = @(t,xTil,Z0)(UcdcMid(t,AxMi.'*xTil,Z0));
        UcdcMidEta = @(tnM1,eta,Z0)(UcdcMidx(tnM1,eta,Z0));
        repositorio.UcdcMidEta = UcdcMidEta;
        %Vetor midpEta:
        fmidpx = @(xTilp)(AxMi.'*xTilp);
        fmidpEta = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(fmidpx(fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)));
        repositorio.fmidpEta = fmidpEta;
        %Vetor fZx:
        fZ = repositorio.fZ;
        fZx = @(xTil,xTilp)(fZ(fuTilcx(xTil),fuTilcx(xTilp)));
        repositorio.fZx = fZx;
        fcXhat = repositorio.fcXhat;
        fcXhatx = @(xTil,xTilp)(fcXhat(fuTilcx(xTil),fuTilcx(xTilp)));
        repositorio.fcXhatx = fcXhatx;
        
        %DERIVADAS:
        dxTilpnM1dEtaT = @(dtnM1)(fdupduTnmk(dtnM1));
        duTilbpdxTilpT = Axu.';
        
        %Matriz dGamaEtadEtaT:
        dGamadfipBd = repositorio.dGamadfipBd;
        dFipBddxTilpT = AxFiBd.';
        dGamaFiBdxdxTilpT = @(xTilp)(dGamadfipBd(AxFiBd.'*xTilp)*dFipBddxTilpT);
        dGamaEtadEtaT = @(dtnM1,eta,xTilbn,xTilbpn,xTilbppn)(...
            dGamaFiBdxdxTilpT(fupnM1(dtnM1,eta,xTilbn,xTilbpn,xTilbppn))*...
            dxTilpnM1dEtaT(dtnM1));
        repositorio.dGamaEtadEtaT = dGamaEtadEtaT;
        %Matrizes dKcddEtaT e dKdcddEtaT:
        DKcDXdT = repositorio.DKcDXdT;
        DKcDXdTx = @(t,xTil,Zd0)(DKcDXdT(t,AxMi.'*xTil,Zd0));
        Ad = argumentos.Ad;
        dXdmi = Ad;
        dmidEtaT = AxMi.';
        dKcEtadEtaT = @(tnM1,eta,Zd0)(DKcDXdTx(tnM1,eta,Zd0)*dXdmi*dmidEtaT);
        repositorio.dKcEtadEtaT = dKcEtadEtaT;
        DKdcDXdT = repositorio.DKdcDXdT;
        DKdcDXdTx = @(t,xTil,Z0)(DKdcDXdT(t,AxMi.'*xTil,Z0));
        dKdcEtadEtaT = @(tnM1,eta,Z0)(DKdcDXdTx(tnM1,eta,Z0)*dXdmi*dmidEtaT);
        repositorio.dKdcEtadEtaT = dKdcEtadEtaT;
        %Matriz dmidpdEtaT:
        fdmidpdEtaT = @(dtnM1)(AxMi.'*dxTilpnM1dEtaT(dtnM1));
        repositorio.fdmidpdEtaT = fdmidpdEtaT;
        %Matrizes dfcZdHatEtadEtaT e dfdcZdHatdEtaT:
        dInvPsiDifcdXTHat = repositorio.dInvPsiDifcdXTHat;
        fXd = repositorio.fXd;
        dInvPsiDifcdXdTHat = @(t,mid,Z0)(dInvPsiDifcdXTHat(fXd(t,mid),Z0));
        dInvPsiDifcdXdTHatx = @(t,xTil,Z0)(dInvPsiDifcdXdTHat(t,AxMi.'*xTil,Z0));
        dInvPsiDifcdXdTHatEta = @(tnM1,eta,Z0)(dInvPsiDifcdXdTHatx(tnM1,eta,Z0));
        dfcZdHatEtadEtaT = @(tnM1,eta,Zd0)(dInvPsiDifcdXdTHatEta(tnM1,eta,Zd0)*dXdmi*dmidEtaT);
        repositorio.dfcZdHatEtadEtaT = dfcZdHatEtadEtaT;
        dInvPsiDifdcdXTHat = repositorio.dInvPsiDifdcdXTHat;
        dInvPsiDifdcdXdTHat = @(t,mid,Z0)(dInvPsiDifdcdXTHat(fXd(t,mid),Z0));
        dInvPsiDifdcdXdTHatx = @(t,xTil,Z0)(dInvPsiDifdcdXdTHat(t,AxMi.'*xTil,Z0));
        dInvPsiDifdcdXdTHatEta = @(tnM1,eta,Z0)(dInvPsiDifdcdXdTHatx(tnM1,eta,Z0));
        dfdcZdHatEtadEtaT = @(tnM1,eta,Z0)(dInvPsiDifdcdXdTHatEta(tnM1,eta,Z0)*dXdmi*dmidEtaT);
        repositorio.dfdcZdHatEtadEtaT = dfdcZdHatEtadEtaT;
        %Matrizes dUccMidEtadEtaT e dUcdcMidEtadEtaT:
        dUccMiddmid = repositorio.dUccMiddmid;
        dUccMiddmid = @(t,xTil,Zd0)(dUccMiddmid(t,AxMi.'*xTil,Zd0));
        dUccMidEtadEtaT = @(tnM1,eta,Zd0)(dUccMiddmid(tnM1,eta,Zd0)*dmidEtaT);
        repositorio.dUccMidEtadEtaT = dUccMidEtadEtaT;
        dUcdcMiddmid = repositorio.dUcdcMiddmid;
        dUcdcMiddmidx = @(t,xTil,Zd0)(dUcdcMiddmid(t,AxMi.'*xTil,Zd0));
        dUcdcMidEtadEtaT = @(tnM1,eta,Zd0)(dUcdcMiddmidx(tnM1,eta,Zd0));
        repositorio.dUcdcMidEtadEtaT = dUcdcMidEtadEtaT;
        
        %Matriz dB6UccTetaPEtadEtaT:
        dB6UccTetaPduTilbT = repositorio.dB6UccTetaPduTilbT;
        dB6UccTetaPduTilbpT = repositorio.dB6UccTetaPduTilbpT;
        dB6UccTetaPdxTilTx = @(t,xTil,xTilp,Zd0)(...
            dB6UccTetaPduTilbT(t,Axu.'*xTil,Axu.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0)*duTilbdxTilT);
        dB6UccTetaPdxTilpTx = @(t,xTil,xTilp,Zd0)(...
            dB6UccTetaPduTilbpT(t,Axu.'*xTil,Axu.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0)*duTilbpdxTilpT);
        dB6UccTetaPdxTilnM1T = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UccTetaPdxTilTx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        dB6UccTetaPdxTilpnM1T = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UccTetaPdxTilpTx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        dB6UccTetaPEtadEtaT = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UccTetaPdxTilnM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0) +...
            dB6UccTetaPdxTilpnM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)*dxTilpnM1dEtaT(dtnM1));
        repositorio.dB6UccTetaPEtadEtaT = dB6UccTetaPEtadEtaT;
        
        %Matriz dB6UcdcTetaPEtadEtaT:
        dB6UcdcTetaPduTilbT = repositorio.dB6UcdcTetaPduTilbT;
        dB6UcdcTetaPduTilbpT = repositorio.dB6UcdcTetaPduTilbpT;
        dB6UcdcTetaPdxTilTx = @(t,xTil,xTilp,Zd0)(...
            dB6UcdcTetaPduTilbT(t,Axu.'*xTil,Axu.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0)*duTilbdxTilT);
        dB6UcdcTetaPdxTilpTx = @(t,xTil,xTilp,Zd0)(...
            dB6UcdcTetaPduTilbpT(t,Axu.'*xTil,Axu.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0)*duTilbpdxTilpT);
        dB6UcdcTetaPdxTilnM1T = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UcdcTetaPdxTilTx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        dB6UcdcTetaPdxTilpnM1T = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UcdcTetaPdxTilpTx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        dB6UcdcTetaPEtadEtaT = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UcdcTetaPdxTilnM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0) +...
            dB6UcdcTetaPdxTilpnM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)*dxTilpnM1dEtaT(dtnM1));
        repositorio.dB6UcdcTetaPEtadEtaT = dB6UcdcTetaPEtadEtaT;
        
        %Função dtetaPEtadEtaT:
        dtetaPubduTilbpT = repositorio.dtetaPubduTilbpT;
        dtetaPxdxTilbpT = @(xTilp)(dtetaPubduTilbpT(Axu.'*xTilp)*duTilbpdxTilpT);
        repositorio.dtetaPxdxTilbpT = dtetaPxdxTilbpT;
        dtetaPEtadEtaT = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            dtetaPxdxTilbpT(fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn))*dxTilpnM1dEtaT(dtnM1));
        repositorio.dtetaPEtadEtaT = dtetaPEtadEtaT;
    case 8
        Ng = argumentos.Ng;
        IdNg = argumentos.IdNg;
        IdNgm2 = IdNg(1:Ng-2,1:Ng-2);
        
        %PLANTA:
        %Matriz Knmk:
        M1 = argumentos.M1;
        C = argumentos.C;
        K = argumentos.K;
        Knmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*M1 + gamaNmk*(betaNmk*dtnM1)^-1*C + K);
        repositorio.Knmk = Knmk;
        %Vetor fTilEta:
        fTilcr = repositorio.fTilcr; 
        fTilcrnM1 = @(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            fTilcr(uTilbrnM1,...
            fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
            fuppnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
            UcnM1,...
            fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn),...
            fuppnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn)));
        %Vetores kTilVcNmk e fiTiln:
        GamaRr = argumentos.GamaRr;
        GamaRvc = argumentos.GamaRvc;
        mTilVc = argumentos.mTilVc;
        cTilVc = argumentos.cTilVc;
        kTilVc = argumentos.kTilVc;
        kTilVcNmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*mTilVc + gamaNmk*(betaNmk*dtnM1)^-1*cTilVc + kTilVc);
        fiTiln = @(dtnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            M1*fi2Nmk(dtnM1,uTilbrn,uTilbprn,uTilbpprn) + C*fi1Nmk(dtnM1,uTilbrn,uTilbprn,uTilbpprn) +...
            mTilVc*fi2Nmk(dtnM1,Ucn,Upcn,Uppcn) + cTilVc*fi1Nmk(dtnM1,Ucn,Upcn,Uppcn));
        
        %SISTEMA EM znM1:
        %Matriz KnmkLmb:
        rEstr = argumentos.rEstr;
        kNmk = argumentos.k_newmark;
        Azu = [IdNgm2; zeros(rEstr,Ng-2)];
        argumentos.Azu = Azu;
        KnmkLmb = @(dtnM1)([Knmk(dtnM1), kNmk*GamaRr.'; 
            kNmk*GamaRr, zeros(rEstr)]);
        repositorio.KnmkLmb = KnmkLmb;
        %Vetor fcLmb:
        IdRestr = eye(rEstr);
        Azlmb = [zeros(Ng-2,rEstr); IdRestr];
        argumentos.Azlmb = Azlmb;
        fcLmb = @(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            [fTilcrnM1(dtnM1,Azu.'*znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn); 
            zeros(rEstr,1)]);
        repositorio.fcLmb = fcLmb;
        %Matriz Blmb:
        Blmb = @(dtnM1)([kTilVcNmk(dtnM1);
            kNmk*GamaRvc]);
        repositorio.Blmk = Blmb;
        %Vetor fiTilnLmb:
        fiTilnLmb = @(dtnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            [fiTiln(dtnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn);
            zeros(rEstr,1)]);
        
        %SISTEMA EM eta:
        Psic = @(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            KnmkLmb(dtnM1)*znM1 +...
            fcLmb(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn) +...
            Blmb(dtnM1)*UcnM1 +...
            fiTilnLmb(dtnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn));
        repositorio.Psic = Psic;
        
        %%
        %Derivadas:
        dfTilcrdubrT = repositorio.dfTilcrdubrT;
        dfTilcrdubprT = repositorio.dfTilcrdubprT;
        dfTilrdubpprT = repositorio.dfTilrdubpprT;
        dfTilcrdUcT = repositorio.dfTilcrdUcT;
        dfTilcrdUpcT = repositorio.dfTilcrdUpcT;
        dfTilrdUppcT = repositorio.dfTilrdUppcT;
        dubrpnM1dubrnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dubrppnM1dubrnM1T = @(dtnM1)(fduppduTnmk(dtnM1));
        
        %Matriz dfcLmbdEtaT:
        dfTilcrdubrnM1T = @(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            dfTilcrdubrT(uTilbrnM1,...
            fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
            fuppnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
            UcnM1,...
            fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn),...
            fuppnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn)) +...
            dfTilcrdubprT(uTilbrnM1,...
            fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
            UcnM1,...
            fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn))*dubrpnM1dubrnM1T(dtnM1) +...
            dfTilrdubpprT(uTilbrnM1,UcnM1)*dubrppnM1dubrnM1T(dtnM1));
        dfTilcrdUcnM1T = @(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            dfTilcrdUcT(uTilbrnM1,...
            fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
            fuppnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
            UcnM1,...
            fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn),...
            fuppnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn)) +...
            dfTilcrdUpcT(uTilbrnM1,...
            fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
            UcnM1,...
            fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn))*dubrpnM1dubrnM1T(dtnM1) +...
            dfTilrdUppcT(uTilbrnM1,UcnM1)*dubrppnM1dubrnM1T(dtnM1));
        dUcdubrT = zeros(2,Ng-2);
        DfTilcrDuTilbrT = @(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            dfTilcrdubrnM1T(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn) +...
            dfTilcrdUcnM1T(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)*dUcdubrT);
        duTilrdzTilT = Azu.';
        dfcLmbdznM1T = @(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            [DfTilcrDuTilbrT(dtnM1,Azu.'*znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)*duTilrdzTilT; 
            zeros(rEstr,Ng-2+rEstr)]);
        dUcnM1dznM1T = dUcdubrT*duTilrdzTilT;
        dPsicdznM1T = @(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            KnmkLmb(dtnM1) +...
            dfcLmbdznM1T(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn) +...
            Blmb(dtnM1)*...
            dUcnM1dznM1T);
        repositorio.dPsicdznM1T = dPsicdznM1T;
        %Obs.: eta == znM1 
    case 9
        Ng = argumentos.Ng;
        IdNg = argumentos.IdNg;
        IdNgm2 = IdNg(1:Ng-2,1:Ng-2);
        
        %PLANTA:
        %Matriz Knmk:
        M1 = argumentos.M1;
        C = argumentos.C;
        K = argumentos.K;
        Knmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*M1 + gamaNmk*(betaNmk*dtnM1)^-1*C + K);
        repositorio.Knmk = Knmk;
        %Vetor fTilEta:
        fTilcr = repositorio.fTilcr; 
        fTilcrnM1 = @(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            fTilcr(uTilbrnM1,...
                   fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
                   fuppnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
                   UcnM1,...
                   fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn),...
                   fuppnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn)));
        %Vetores kTilVcNmk e fiTiln:
        GamaRr = argumentos.GamaRr;
        GamaRvc = argumentos.GamaRvc;
        mTilVc = argumentos.mTilVc;
        cTilVc = argumentos.cTilVc;
        kTilVc = argumentos.kTilVc;
        kTilVcNmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*mTilVc + gamaNmk*(betaNmk*dtnM1)^-1*cTilVc + kTilVc);
        fiTiln = @(dtnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            M1*fi2Nmk(dtnM1,uTilbrn,uTilbprn,uTilbpprn) + C*fi1Nmk(dtnM1,uTilbrn,uTilbprn,uTilbpprn) +...
            mTilVc*fi2Nmk(dtnM1,Ucn,Upcn,Uppcn) + cTilVc*fi1Nmk(dtnM1,Ucn,Upcn,Uppcn));
        
        %SISTEMA EM znM1:
        %Matriz KnmkLmb:
        rEstr = argumentos.rEstr;
        kNmk = argumentos.k_newmark;
        Azu = [IdNgm2; zeros(rEstr,Ng-2)];
        argumentos.Azu = Azu;
        KnmkLmb = @(dtnM1)([Knmk(dtnM1), kNmk*GamaRr.'; 
            kNmk*GamaRr, zeros(rEstr)]);
        repositorio.KnmkLmb = KnmkLmb;
        %Vetor fcLmb:
        IdRestr = eye(rEstr);
        Azlmb = [zeros(Ng-2,rEstr); IdRestr];
        argumentos.Azlmb = Azlmb;
        fcLmb = @(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            [fTilcrnM1(dtnM1,Azu.'*znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn); 
            zeros(rEstr,1)]);
        repositorio.fcLmb = fcLmb;
        %Matriz Blmb:
        Blmb = @(dtnM1)([kTilVcNmk(dtnM1);
            kNmk*GamaRvc]);
        repositorio.Blmk = Blmb;
        %Vetor fiTilnLmb:
        fiTilnLmb = @(dtnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            [fiTiln(dtnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn);
            zeros(rEstr,1)]);
        
        %SISTEMA EM eta:
        Psic = @(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            KnmkLmb(dtnM1)*znM1 +...
            fcLmb(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn) +...
            Blmb(dtnM1)*UcnM1 +...
            fiTilnLmb(dtnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn));
        repositorio.Psic = Psic;
        
        %%
        %Derivadas:
        dfTilcrdubrT = repositorio.dfTilcrdubrT;
        dfTilcrdubprT = repositorio.dfTilcrdubprT;
        dfTilrdubpprT = repositorio.dfTilrdubpprT;
        dfTilcrdUcT = repositorio.dfTilcrdUcT;
        dfTilcrdUpcT = repositorio.dfTilcrdUpcT;
        dfTilrdUppcT = repositorio.dfTilrdUppcT;
        dubrpnM1dubrnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dubrppnM1dubrnM1T = @(dtnM1)(fduppduTnmk(dtnM1));
        
        %Matriz dfcLmbdEtaT:
        dfTilcrdubrnM1T = @(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            dfTilcrdubrT(uTilbrnM1,...
                         fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
                         fuppnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
                         UcnM1,...
                         fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn),...
                         fuppnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn)) +...
            dfTilcrdubprT(uTilbrnM1,...
                          fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
                          UcnM1,...
                          fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn))*dubrpnM1dubrnM1T(dtnM1) +...
            dfTilrdubpprT(uTilbrnM1,UcnM1)*dubrppnM1dubrnM1T(dtnM1));
        dfTilcrdUcnM1T = @(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            dfTilcrdUcT(uTilbrnM1,...
                        fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
                        fuppnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
                        UcnM1,...
                        fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn),...
                        fuppnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn)) +...
            dfTilcrdUpcT(uTilbrnM1,...
                         fupnM1(dtnM1,uTilbrnM1,uTilbrn,uTilbprn,uTilbpprn),...
                         UcnM1,...
                         fupnM1(dtnM1,UcnM1,Ucn,Upcn,Uppcn))*dubrpnM1dubrnM1T(dtnM1) +...
            dfTilrdUppcT(uTilbrnM1,UcnM1)*dubrppnM1dubrnM1T(dtnM1));
        dUcdubrT = zeros(2,Ng-2);
        DfTilcrDuTilbrT = @(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            dfTilcrdubrnM1T(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn) +...
            dfTilcrdUcnM1T(dtnM1,uTilbrnM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)*dUcdubrT);
        duTilrdzTilT = Azu.';
        dfcLmbdznM1T = @(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            [DfTilcrDuTilbrT(dtnM1,Azu.'*znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)*duTilrdzTilT; 
            zeros(rEstr,Ng-2+rEstr)]);
        dUcnM1dznM1T = dUcdubrT*duTilrdzTilT;
        dPsicdznM1T = @(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn)(...
            KnmkLmb(dtnM1) +...
            dfcLmbdznM1T(dtnM1,znM1,UcnM1,uTilbrn,uTilbprn,uTilbpprn,Ucn,Upcn,Uppcn) +...
            Blmb(dtnM1)*...
            dUcnM1dznM1T);
        repositorio.dPsicdznM1T = dPsicdznM1T;
        %Obs.: eta == znM1
    case 10
        Id8 = eye(8);
        Axu = Id8(:,1:6);
        argumentos.Axu = Axu;
        AxFiBd = Id8(:,7);
        argumentos.AxFiBd = AxFiBd;
        AxMi = Id8(:,8);
        argumentos.AxMi = AxMi;
        
        %PLANTA:
        %Matriz Knmk:
        M1 = argumentos.M1;
        C = argumentos.C;
        K = argumentos.K;
        Knmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*M1 + gamaNmk*(betaNmk*dtnM1)^-1*C + K);
        repositorio.Knmk = Knmk;
        KnmkEta = @(dtnM1)(Knmk(dtnM1)*Axu.');
        repositorio.KnmkEta = KnmkEta;
        M1x = M1*Axu.';
        argumentos.M1x = M1x;
        Cx = C*Axu.';
        argumentos.Cx = Cx;
        %Vetor fTilEta:
        fTilc = repositorio.fTilc;
        fTilcx = @(xTil,xTilp,xTilpp)(fTilc(Axu.'*xTil,Axu.'*xTilp,Axu.'*xTilpp));
        fTilcEta = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            fTilcx(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),...
            fuppnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)));
        repositorio.fTilcEta = fTilcEta;
        %Vetor B6Eta:
        B6 = repositorio.B6;
        B6x = @(xTil)(B6(Axu.'*xTil));
        B6Eta = @(eta)(B6x(eta));
        repositorio.B6Eta = B6Eta;
        %Matriz DfTilcEtaDEtaT:
        duTilbdxTilT = Axu.';
        duTilbpdxTilpT = Axu.';
        duTilbppdxTilppT = Axu.';
        dxTilpnM1dxTilnM1 = @(dtnM1)(fdupduTnmk(dtnM1));
        dxTilppnM1dxTilnM1 = @(dtnM1)(fduppduTnmk(dtnM1));
        dfTilcduTilbT = repositorio.dfTilcdubT;
        dfTilcduTilbpT = repositorio.dfTilcdubpT;
        dfTilduTilbppT = repositorio.dfTildubppT;
        dfTilcxdxTilT = @(xTil,xTilp,xTilpp)(dfTilcduTilbT(Axu.'*xTil,Axu.'*xTilp,Axu.'*xTilpp)*duTilbdxTilT);
        dfTilcxdxTilpT = @(xTil,xTilp,xTilpp)(dfTilcduTilbpT(Axu.'*xTil,Axu.'*xTilp)*duTilbpdxTilpT);
        dfTilxdxTilppT = @(xTil,xTilp,xTilpp)(dfTilduTilbppT(Axu.'*xTil)*duTilbppdxTilppT);
        DfTilcEtaDetaT = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            dfTilcxdxTilT(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),...
            fuppnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)) +...
            dfTilcxdxTilpT(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn))*...
            dxTilpnM1dxTilnM1(dtnM1) +...
            dfTilxdxTilppT(eta)*dxTilppnM1dxTilnM1(dtnM1));
        repositorio.DfTilcEtaDetaT = DfTilcEtaDetaT;
        %Matriz DfTildcEtaDEtaT:
        dfTildcduTilbT = repositorio.dfTildcdubT;
        dfTildcduTilbpT = repositorio.dfTildcdubpT;
        dfTildcxdxTilT = @(xTil,xTilp,xTilpp)(dfTildcduTilbT(Axu.'*xTil,Axu.'*xTilp,Axu.'*xTilpp)*duTilbdxTilT);
        dfTildcxdxTilpT = @(xTil,xTilp,xTilpp)(dfTildcduTilbpT(Axu.'*xTil,Axu.'*xTilp)*duTilbpdxTilpT);
        DfTildcEtaDetaT = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            dfTildcxdxTilT(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),...
            fuppnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)) +...
            dfTildcxdxTilpT(eta,...
            fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn))*...
            dxTilpnM1dxTilnM1(dtnM1) +...
            dfTilxdxTilppT(eta)*dxTilppnM1dxTilnM1(dtnM1));
        repositorio.DfTildcEtaDetaT = DfTildcEtaDetaT;
        
        %LEI DE CONTROLE:
        %Variáveis do vetor posição (lei de controle): uTilc = [teta2; teta3]
        fuTilc = repositorio.fuTilc;
        fuTilcx = @(xTil)(fuTilc(Axu.'*xTil));
        %Vetores UccEta e UcdcEta:
        Uccu = repositorio.Uccu;
        Uccx = @(t,xTil,xTilp,Zd0)(Uccu(t,fuTilcx(xTil),fuTilcx(xTilp),AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        repositorio.Uccx = Uccx;
        UccEta = @(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            Uccx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        repositorio.UccEta = UccEta;
        Ucdcu = repositorio.Ucdcu;
        Ucdcx = @(t,xTil,xTilp,Zd0)(Ucdcu(t,fuTilcx(xTil),fuTilcx(xTilp),AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        repositorio.Ucdcx = Ucdcx;
        UcdcEta = @(tnM1,dtnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            Ucdcx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        repositorio.UcdcEta = UcdcEta;
        %Vetor Upccx e Updccx:
        Upcc = repositorio.fUpcc; 
        Upccx = @(t,xTil,xTilp,Z0)(Upcc(t,fuTilcx(xTil),AxFiBd.'*xTil,AxMi.'*xTil,fuTilcx(xTilp),AxFiBd.'*xTilp,AxMi.'*xTilp,Z0));
        repositorio.Upccx = Upccx;
        Upcdc = repositorio.fUpcdc;
        Upcdcx = @(t,xTil,xTilp,Z0)(Upcdc(t,fuTilcx(xTil),AxFiBd.'*xTil,AxMi.'*xTil,fuTilcx(xTilp),AxFiBd.'*xTilp,AxMi.'*xTilp,Z0));
        repositorio.Upcdcx = Upcdcx;
        
        %Função tetaPEta:
        tetaPub = repositorio.tetaPub;
        tetaPx = @(xTilp)(tetaPub(Axu.'*xTilp));
        repositorio.tetaPx = tetaPx;
        tetaPEta = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(tetaPx(fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)));
        repositorio.tetaPEta = tetaPEta;
        
        %Vetor GamaFiEta:
        GamaFiBd = repositorio.GamaFiBd;
        GamaFiBdx = @(xTilp)(GamaFiBd(AxFiBd.'*xTilp));
        GamaFiBdEta = @(dtnM1,eta,xTilbn,xTilbpn,xTilbppn)(GamaFiBdx(fupnM1(dtnM1,eta,xTilbn,xTilbpn,xTilbppn)));
        repositorio.GamaFiBdEta = GamaFiBdEta;
        fFippBd = repositorio.fFippBd;
        fFippBdx = @(t,xTil,xTilp,Zd0)(fFippBd(t,AxFiBd.'*xTilp,AxMi.'*xTil,Zd0));
        repositorio.fFippBdx = fFippBdx;
        
        %Matriz Kfix:
        Kfi = argumentos.Kfi;
        Kfix = Kfi*AxFiBd.';
        argumentos.Kfix = Kfix;
        %Vetores KcdEta e KdcdEta:
        Kcd = repositorio.Kcd;
        Kcdx = @(t,xTil,Zd0)(Kcd(t,AxMi.'*xTil,Zd0));
        KcdEta = @(tnM1,eta,Zd0)(Kcdx(tnM1,eta,Zd0));
        repositorio.KcdEta = KcdEta;
        Kdcd = repositorio.Kdcd;
        Kdcdx = @(t,xTil,Z0)(Kdcd(t,AxMi.'*xTil,Z0));
        KdcdEta = @(tnM1,eta,Z0)(Kdcdx(tnM1,eta,Z0));
        repositorio.KdcdEta = KdcdEta;
        %Vetores fcZdHatEta e fdcZdHatEta:
        fcZdHat = repositorio.fcZdHat;
        fcZdHatx = @(t,xTil,Zd0)(fcZdHat(t,AxMi.'*xTil,Zd0));
        repositorio.fcZdHatx = fcZdHatx;
        fcZdHatEta = @(tnM1,eta,Zd0)(fcZdHatx(tnM1,eta,Zd0));
        repositorio.fcZdHatEta = fcZdHatEta;
        fdcZdHat = repositorio.fdcZdHat;
        fdcZdHatx = @(t,xTil,Z0)(fdcZdHat(t,AxMi.'*xTil,Z0));
        fdcZdHatEta = @(tnM1,eta,Z0)(fdcZdHatx(tnM1,eta,Z0));
        repositorio.fdcZdHatEta = fdcZdHatEta;
        %Expressões Uccd e Ucdcd:
        UccMid = repositorio.UccMid;
        UccMidx = @(t,xTil,Z0)(UccMid(t,AxMi.'*xTil,Z0));
        UccMidEta = @(tnM1,eta,Zd0)(UccMidx(tnM1,eta,Zd0));
        repositorio.UccMidEta = UccMidEta;
        UcdcMid = repositorio.UcdcMid;
        UcdcMidx = @(t,xTil,Z0)(UcdcMid(t,AxMi.'*xTil,Z0));
        UcdcMidEta = @(tnM1,eta,Z0)(UcdcMidx(tnM1,eta,Z0));
        repositorio.UcdcMidEta = UcdcMidEta;
        %Vetor midpEta:
        fmidpx = @(xTilp)(AxMi.'*xTilp);
        fmidpEta = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(fmidpx(fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)));
        repositorio.fmidpEta = fmidpEta;
        %Vetor fZx:
        fZ = repositorio.fZ;
        fZx = @(xTil,xTilp)(fZ(fuTilcx(xTil),fuTilcx(xTilp)));
        repositorio.fZx = fZx;
        fcXhat = repositorio.fcXhat;
        fcXhatx = @(xTil,xTilp)(fcXhat(fuTilcx(xTil),fuTilcx(xTilp)));
        repositorio.fcXhatx = fcXhatx;
        
        %DERIVADAS:
        dxTilpnM1dEtaT = @(dtnM1)(fdupduTnmk(dtnM1));
        duTilbpdxTilpT = Axu.';
        
        %Matriz dGamaEtadEtaT:
        dGamadfipBd = repositorio.dGamadfipBd;
        dFipBddxTilpT = AxFiBd.';
        dGamaFiBdxdxTilpT = @(xTilp)(dGamadfipBd(AxFiBd.'*xTilp)*dFipBddxTilpT);
        dGamaEtadEtaT = @(dtnM1,eta,xTilbn,xTilbpn,xTilbppn)(...
            dGamaFiBdxdxTilpT(fupnM1(dtnM1,eta,xTilbn,xTilbpn,xTilbppn))*...
            dxTilpnM1dEtaT(dtnM1));
        repositorio.dGamaEtadEtaT = dGamaEtadEtaT;
        %Matrizes dKcddEtaT e dKdcddEtaT:
        DKcDXdT = repositorio.DKcDXdT;
        DKcDXdTx = @(t,xTil,Zd0)(DKcDXdT(t,AxMi.'*xTil,Zd0));
        Ad = argumentos.Ad;
        dXdmi = Ad;
        dmidEtaT = AxMi.';
        dKcEtadEtaT = @(tnM1,eta,Zd0)(DKcDXdTx(tnM1,eta,Zd0)*dXdmi*dmidEtaT);
        repositorio.dKcEtadEtaT = dKcEtadEtaT;
        DKdcDXdT = repositorio.DKdcDXdT;
        DKdcDXdTx = @(t,xTil,Z0)(DKdcDXdT(t,AxMi.'*xTil,Z0));
        dKdcEtadEtaT = @(tnM1,eta,Z0)(DKdcDXdTx(tnM1,eta,Z0)*dXdmi*dmidEtaT);
        repositorio.dKdcEtadEtaT = dKdcEtadEtaT;
        %Matriz dmidpdEtaT:
        fdmidpdEtaT = @(dtnM1)(AxMi.'*dxTilpnM1dEtaT(dtnM1));
        repositorio.fdmidpdEtaT = fdmidpdEtaT;
        %Matrizes dfcZdHatEtadEtaT e dfdcZdHatdEtaT:
        dInvPsiDifcdXTHat = repositorio.dInvPsiDifcdXTHat;
        fXd = repositorio.fXd;
        dInvPsiDifcdXdTHat = @(t,mid,Z0)(dInvPsiDifcdXTHat(fXd(t,mid),Z0));
        dInvPsiDifcdXdTHatx = @(t,xTil,Z0)(dInvPsiDifcdXdTHat(t,AxMi.'*xTil,Z0));
        dInvPsiDifcdXdTHatEta = @(tnM1,eta,Z0)(dInvPsiDifcdXdTHatx(tnM1,eta,Z0));
        dfcZdHatEtadEtaT = @(tnM1,eta,Zd0)(dInvPsiDifcdXdTHatEta(tnM1,eta,Zd0)*dXdmi*dmidEtaT);
        repositorio.dfcZdHatEtadEtaT = dfcZdHatEtadEtaT;
        dInvPsiDifdcdXTHat = repositorio.dInvPsiDifdcdXTHat;
        dInvPsiDifdcdXdTHat = @(t,mid,Z0)(dInvPsiDifdcdXTHat(fXd(t,mid),Z0));
        dInvPsiDifdcdXdTHatx = @(t,xTil,Z0)(dInvPsiDifdcdXdTHat(t,AxMi.'*xTil,Z0));
        dInvPsiDifdcdXdTHatEta = @(tnM1,eta,Z0)(dInvPsiDifdcdXdTHatx(tnM1,eta,Z0));
        dfdcZdHatEtadEtaT = @(tnM1,eta,Z0)(dInvPsiDifdcdXdTHatEta(tnM1,eta,Z0)*dXdmi*dmidEtaT);
        repositorio.dfdcZdHatEtadEtaT = dfdcZdHatEtadEtaT;
        %Matrizes dUccMidEtadEtaT e dUcdcMidEtadEtaT:
        dUccMiddmid = repositorio.dUccMiddmid;
        dUccMiddmid = @(t,xTil,Zd0)(dUccMiddmid(t,AxMi.'*xTil,Zd0));
        dUccMidEtadEtaT = @(tnM1,eta,Zd0)(dUccMiddmid(tnM1,eta,Zd0)*dmidEtaT);
        repositorio.dUccMidEtadEtaT = dUccMidEtadEtaT;
        dUcdcMiddmid = repositorio.dUcdcMiddmid;
        dUcdcMiddmidx = @(t,xTil,Zd0)(dUcdcMiddmid(t,AxMi.'*xTil,Zd0));
        dUcdcMidEtadEtaT = @(tnM1,eta,Zd0)(dUcdcMiddmidx(tnM1,eta,Zd0));
        repositorio.dUcdcMidEtadEtaT = dUcdcMidEtadEtaT;
        
        %Matriz dB6UccTetaPEtadEtaT:
        dB6UccTetaPduTilbT = repositorio.dB6UccTetaPduTilbT;
        dB6UccTetaPduTilbpT = repositorio.dB6UccTetaPduTilbpT;
        dB6UccTetaPdxTilTx = @(t,xTil,xTilp,Zd0)(...
            dB6UccTetaPduTilbT(t,Axu.'*xTil,Axu.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0)*duTilbdxTilT);
        dB6UccTetaPdxTilpTx = @(t,xTil,xTilp,Zd0)(...
            dB6UccTetaPduTilbpT(t,Axu.'*xTil,Axu.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0)*duTilbpdxTilpT);
        dB6UccTetaPdxTilnM1T = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UccTetaPdxTilTx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        dB6UccTetaPdxTilpnM1T = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UccTetaPdxTilpTx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        dB6UccTetaPEtadEtaT = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UccTetaPdxTilnM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0) +...
            dB6UccTetaPdxTilpnM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)*dxTilpnM1dEtaT(dtnM1));
        repositorio.dB6UccTetaPEtadEtaT = dB6UccTetaPEtadEtaT;
        
        %Matriz dB6UcdcTetaPEtadEtaT:
        dB6UcdcTetaPduTilbT = repositorio.dB6UcdcTetaPduTilbT;
        dB6UcdcTetaPduTilbpT = repositorio.dB6UcdcTetaPduTilbpT;
        dB6UcdcTetaPdxTilTx = @(t,xTil,xTilp,Zd0)(...
            dB6UcdcTetaPduTilbT(t,Axu.'*xTil,Axu.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0)*duTilbdxTilT);
        dB6UcdcTetaPdxTilpTx = @(t,xTil,xTilp,Zd0)(...
            dB6UcdcTetaPduTilbpT(t,Axu.'*xTil,Axu.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0)*duTilbpdxTilpT);
        dB6UcdcTetaPdxTilnM1T = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UcdcTetaPdxTilTx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        dB6UcdcTetaPdxTilpnM1T = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UcdcTetaPdxTilpTx(tnM1,eta,fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn),Zd0));
        dB6UcdcTetaPEtadEtaT = @(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)(...
            dB6UcdcTetaPdxTilnM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0) +...
            dB6UcdcTetaPdxTilpnM1T(dtnM1,tnM1,eta,xTiln,xTilpn,xTilppn,Zd0)*dxTilpnM1dEtaT(dtnM1));
        repositorio.dB6UcdcTetaPEtadEtaT = dB6UcdcTetaPEtadEtaT;
        
        %Função dtetaPEtadEtaT:
        dtetaPubduTilbpT = repositorio.dtetaPubduTilbpT;
        dtetaPxdxTilbpT = @(xTilp)(dtetaPubduTilbpT(Axu.'*xTilp)*duTilbpdxTilpT);
        repositorio.dtetaPxdxTilbpT = dtetaPxdxTilbpT;
        dtetaPEtadEtaT = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            dtetaPxdxTilbpT(fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn))*dxTilpnM1dEtaT(dtnM1));
        repositorio.dtetaPEtadEtaT = dtetaPEtadEtaT;
    case 11
        %%
        Ng = argumentos.Ng;
        IdNg = argumentos.IdNg;
        rEstr = argumentos.rEstr;
        IdNgMrEstr = eye(Ng+rEstr);
        Azx = IdNgMrEstr(:,1:Ng);
        argumentos.Azx = Azx;
        Azlmb = IdNgMrEstr(:,Ng+1:Ng+rEstr);
        argumentos.Azlmb = Azlmb;
        
        %Decomposi��o do vetor xTil (Ng x 1):
        %Obs.: 
        %(1) xTil = Axur*ur + AxFiBd*fiBd + AxMi*mid
        %(2) ur: Ng-2 x 1        
        Axur = IdNg(:,1:Ng-2);
        argumentos.Axur = Axur;
        AxFiBd = IdNg(:,Ng-1);
        argumentos.AxFiBd = AxFiBd;
        AxMi = IdNg(:,Ng);
        argumentos.AxMi = AxMi;
        
        %PLANTA:
        %Matriz Knmkx:
        M1 = argumentos.M1;
        C = argumentos.C;
        K = argumentos.K;
        M1x = M1*Axur.';
        Cx = C*Axur.';
        Kx = K*Axur.';
        Knmkx = @(dtnM1)((betaNmk*dtnM1^2)^-1*M1x + gamaNmk*(betaNmk*dtnM1)^-1*Cx + Kx);
        repositorio.Knmkx = Knmkx;
        
        %Vetores kTilVcNmk e fiTilxn:
        GamaRr = argumentos.GamaRr;
        GamaRx = GamaRr*Axur.';
        GamaRvc = argumentos.GamaRvc;
        mTilVc = argumentos.mTilVc;
        cTilVc = argumentos.cTilVc;
        kTilVc = argumentos.kTilVc;
        kTilVcNmk = @(dtnM1)((betaNmk*dtnM1^2)^-1*mTilVc + gamaNmk*(betaNmk*dtnM1)^-1*cTilVc + kTilVc);
        fiTilxn = @(dtnM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn)(...
            M1x*fi2Nmk(dtnM1,xTiln,xTilpn,xTilppn) + Cx*fi1Nmk(dtnM1,xTiln,xTilpn,xTilppn) +...
            mTilVc*fi2Nmk(dtnM1,Ucn,Upcn,Uppcn) + cTilVc*fi1Nmk(dtnM1,Ucn,Upcn,Uppcn));
        
        %Matriz KnmkLmb:
        kNmk = argumentos.k_newmark;
        KnmkLmb = @(dtnM1)([Knmkx(dtnM1), kNmk*GamaRr.'; 
            kNmk*GamaRx, zeros(rEstr)]);
        repositorio.KnmkLmb = KnmkLmb;
        
        %Vetor fTilxEta:
        fTilcr = repositorio.fTilcr;
        fUccr = repositorio.fUccr;
        fUccx = @(t,xTil,xTilp,Zd0)(...
            fUccr(t,Axur.'*xTil,Axur.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        fUpccr = repositorio.fUpccr;
        fUpccx = @(t,xTil,xTilp,xTilpp,Zd0)(...
            fUpccr(t,...
                   Axur.'*xTil,Axur.'*xTilp,Axur.'*xTilpp,...
                   AxFiBd.'*xTil,AxFiBd.'*xTilp,...
                   AxMi.'*xTil,AxMi.'*xTilp,...
                   Zd0));
        repositorio.fUpccx = fUpccx;
        fTilcx = @(xTil,xTilp,xTilpp,Uc,Upc,Uppc)(...
            fTilcr(Axur.'*xTil,Axur.'*xTilp,Axur.'*xTilpp,Uc,Upc,Uppc));
        fTilcxnM1 = @(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn)(...
            fTilcx(xTilnM1,...
                   fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                   fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),...
                   fupnM1(dtnM1,fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),Ucn,Upcn,Uppcn),...
                   fuppnM1(dtnM1,fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),Ucn,Upcn,Uppcn)));
        fTilcLmb = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn)(...
            fTilcxnM1(dtnM1,tnM1,Azx.'*znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn));
                       
        %Vetor fcLmb:
        fcLmb = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn)(...
            [fTilcLmb(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn);
            zeros(rEstr,1)]);
        repositorio.fcLmb = fcLmb;
        
        %Matriz Blmb:
        Blmb = @(dtnM1)([kTilVcNmk(dtnM1);
            kNmk*GamaRvc]);
        repositorio.Blmk = Blmb;
        
        %Vetor fiTilxnLmb:
        fiTilxnLmb = @(dtnM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn)(...
            [fiTilxn(dtnM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn);
            zeros(rEstr,1)]);
        
        %Vetor GamaFiBdLmb:
        GamaFiBd = repositorio.GamaFiBd;
        GamaFiBdx = @(xTilp)(GamaFiBd(AxFiBd.'*xTilp));
        GamaFiBdxnM1 = @(dtnM1,xTilbnM1,xTilbn,xTilbpn,xTilbppn)(...
            GamaFiBdx(fupnM1(dtnM1,xTilbnM1,xTilbn,xTilbpn,xTilbppn)));
        repositorio.GamaFiBdxnM1 = GamaFiBdxnM1;
        GamaFiBdLmb = @(dtnM1,znM1,xTilbn,xTilbpn,xTilbppn)(...
            GamaFiBdx(fupnM1(dtnM1,Azx.'*znM1,xTilbn,xTilbpn,xTilbppn)));
        fFippBd = repositorio.fFippBd;
        fFippBdx = @(t,xTil,xTilp,Zd0)(fFippBd(t,AxFiBd.'*xTilp,AxMi.'*xTil,Zd0));
        repositorio.fFippBdx = fFippBdx;
        
        %Matriz KfiLmb:
        Kfi = argumentos.Kfi;
        Kfix = Kfi*AxFiBd.';
        KfiLmb = Kfix*Azx.';
        fcZdHat = repositorio.fcZdHat;
        fcZdHatx = @(t,xTil,Zd0)(fcZdHat(t,AxMi.'*xTil,Zd0));
        repositorio.fcZdHatx = fcZdHatx;
        fcXhat = repositorio.fcXhat;
        fuTilc = repositorio.fuTilc;
        fuTilcx = @(xTil)(fuTilc(Axur.'*xTil));
        fcXhatx = @(xTil,xTilp)(fcXhat(fuTilcx(xTil),fuTilcx(xTilp)));
        repositorio.fcXhatx = fcXhatx;
        
        %Vetores KcdEta e KdcdEta:
        Kcd = repositorio.Kcd;
        Kcdx = @(t,xTil,Zd0)(Kcd(t,AxMi.'*xTil,Zd0));
        KcdxnM1 = @(tnM1,xTilnM1,Zdn)(Kcdx(tnM1,xTilnM1,Zdn));
        KcdLmb = @(tnM1,znM1,Zdn)(KcdxnM1(tnM1,Azx.'*znM1,Zdn));
        repositorio.KcdLmb = KcdLmb;
        Kdcd = repositorio.Kdcd;
        Kdcdx = @(t,xTil,Zd0)(Kdcd(t,AxMi.'*xTil,Zd0));
        KdcdxnM1 = @(tnM1,xTilnM1,Zdn)(Kdcdx(tnM1,xTilnM1,Zdn));
        KdcdLmb = @(tnM1,znM1,Zdn)(KdcdxnM1(tnM1,Azx.'*znM1,Zdn));
        repositorio.KdcdLmb = KdcdLmb;
        
        %Vetor midpEta:
        fmidpx = @(xTilp)(AxMi.'*xTilp);
        fmidpxnM1 = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(fmidpx(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)));
        fmidpLmb = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(fmidpxnM1(dtnM1,Azx.'*znM1,xTiln,xTilpn,xTilppn));
        repositorio.fmidpLmb = fmidpLmb;
        
        %Expressão da dinâmica interna:
        fcMidpHat = repositorio.fcMidpHat;
        fcMidpHatx = @(t,xTil,Zd0)(fcMidpHat(t,AxMi.'*xTil,Zd0));
        fcMidpHatxnM1 = @(tnM1,xTilnM1,Zdn)(fcMidpHatx(tnM1,xTilnM1,Zdn));
        fcMidpHatLmb = @(tnM1,znM1,Zdn)(fcMidpHatxnM1(tnM1,Azx.'*znM1,Zdn));
        repositorio.fcMidpHatLmb = fcMidpHatLmb;
        fdcMidpHat = repositorio.fdcMidpHat;
        fdcMidpHatx = @(t,xTil,Zd0)(fdcMidpHat(t,AxMi.'*xTil,Zd0));
        fdcMidpHatxnM1 = @(tnM1,xTilnM1,Zdn)(fdcMidpHatx(tnM1,xTilnM1,Zdn));
        fdcMidpHatLmb = @(tnM1,znM1,Zdn)(fdcMidpHatxnM1(tnM1,Azx.'*znM1,Zdn));
        repositorio.fdcMidpHatLmb = fdcMidpHatLmb;
        
        %SISTEMA EM znM1:
        fUccxnM1 = @(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn)(...
            fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn));
        repositorio.fUccxnM1 = fUccxnM1;
        fUccLmb = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Zdn)(fUccxnM1(dtnM1,tnM1,Azx.'*znM1,xTiln,xTilpn,xTilppn,Zdn));
        Psic = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn)(...
            [KnmkLmb(dtnM1)*znM1 + fcLmb(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn) + Blmb(dtnM1)*fUccLmb(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Zdn) + fiTilxnLmb(dtnM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn);
            GamaFiBdLmb(dtnM1,znM1,xTiln,xTilpn,xTilppn) + KfiLmb*znM1 - KcdLmb(tnM1,znM1,Zdn);
            fmidpLmb(dtnM1,znM1,xTiln,xTilpn,xTilppn) - fcMidpHatLmb(tnM1,znM1,Zdn)]);
        repositorio.Psic = Psic;
        
        %%
        %DERIVADAS:
        dfTilcrdubrT = repositorio.dfTilcrdubrT;
        dfTilcrdubprT = repositorio.dfTilcrdubprT;
        dfTilrdubpprT = repositorio.dfTilrdubpprT;
        dfTilcrdUcT = repositorio.dfTilcrdUcT;
        dfTilcrdUpcT = repositorio.dfTilcrdUpcT;
        dfTilrdUppcT = repositorio.dfTilrdUppcT;
        
        dfTilcrdubrTx = @(xTil,xTilp,xTilpp,Uc,Upc,Uppc)(...
            dfTilcrdubrT(Axur.'*xTil,Axur.'*xTilp,Axur.'*xTilpp,Uc,Upc,Uppc));
        dfTilcrdubprTx = @(xTil,xTilp,Uc,Upc)(...
            dfTilcrdubprT(Axur.'*xTil,Axur.'*xTilp,Uc,Upc));
        dfTilcrdubpprTx = @(xTil,Uc)(...
            dfTilrdubpprT(Axur.'*xTil,Uc));
        dfTilcrdUcTx = @(xTil,xTilp,xTilpp,Uc,Upc,Uppc)(...
            dfTilcrdUcT(Axur.'*xTil,Axur.'*xTilp,Axur.'*xTilpp,Uc,Upc,Uppc));
        dfTilcrdUpcTx = @(xTil,xTilp,Uc,Upc)(...
            dfTilcrdUpcT(Axur.'*xTil,Axur.'*xTilp,Uc,Upc));
        dfTilcrdUppcTx = @(xTil,Uc)(...
            dfTilrdUppcT(Axur.'*xTil,Uc));
        
        %LEI DE CONTROLE
        dubrdxT = Axur.';
        dubprdxpT = Axur.';
        dubpprdxppT = Axur.';
        dxpnM1dxnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dxppnM1dxnM1T = @(dtnM1)(fduppduTnmk(dtnM1));
        dUpcnM1dUcnM1T = @(dtnM1)(fdupduTnmk(dtnM1));
        dUppcnM1dUcnM1T = @(dtnM1)(fduppduTnmk(dtnM1));
        
        %Vari�veis do vetor posi��o (lei de controle): uTilc = [teta2; teta3]
        dUccrdurT = repositorio.fdUccrdurT;
        dUccrduprT = repositorio.fdUccrduprT;
        fdUccrdFiBd = repositorio.fdUccrdFiBd;
        fdUccrdmid = repositorio.fdUccrdmid;
        dUccrdurTx = @(t,xTil,xTilp,Zd0)(...
            dUccrdurT(t,Axur.'*xTil,Axur.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        dUccrduprTx = @(t,xTil,xTilp,Zd0)(...
            dUccrduprT(t,Axur.'*xTil,Axur.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        fdUccrdFiBdx = @(t,xTil,xTilp,Zd0)(...
            fdUccrdFiBd(t,Axur.'*xTil,Axur.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        fdUccrdmidx = @(t,xTil,xTilp,Zd0)(...
            fdUccrdmid(t,Axur.'*xTil,Axur.'*xTilp,AxFiBd.'*xTil,AxMi.'*xTil,Zd0));
        
        dfiBddxT = AxFiBd.';
        dmiddxT = AxMi.';
        dUccdxT = @(t,xTil,xTilp,Zd0)(dUccrdurTx(t,xTil,xTilp,Zd0)*dubrdxT +...
            fdUccrdFiBdx(t,xTil,xTilp,Zd0)*dfiBddxT +...
            fdUccrdmidx(t,xTil,xTilp,Zd0)*dmiddxT);
        dUccdxpT = @(t,xTil,xTilp,Zd0)(dUccrduprTx(t,xTil,xTilp,Zd0)*dubprdxpT);
        
        dUccdxnM1T = @(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn)(...
            dUccdxT(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn));
        dUccdxpnM1T = @(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn)(...
            dUccdxpT(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn));
        DUccDxnM1T = @(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn)(...
            dUccdxnM1T(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn) +...
            dUccdxpnM1T(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn)*dxpnM1dxnM1T(dtnM1));
        
        dmiddxT = AxMi.';
        dmipddxpT = AxMi.';
        dfTilxdxnM1T = @(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn)(...
            dfTilcrdubrTx(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),...
                          fupnM1(dtnM1,fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),Ucn,Upcn,Uppcn),...
                          fuppnM1(dtnM1,fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),Ucn,Upcn,Uppcn))*...
            dubrdxT +...
            dfTilcrdubprTx(xTilnM1,...
                           fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                           fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),...
                           fupnM1(dtnM1,fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),Ucn,Upcn,Uppcn))*...
            dubprdxpT*dxpnM1dxnM1T(dtnM1) +...
            dfTilcrdubpprTx(xTilnM1,...
                            fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn))*...
            dubpprdxppT*dxppnM1dxnM1T(dtnM1) +...
            dfTilcrdUcTx(xTilnM1,...
                         fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fuppnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                         fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),...
                         fupnM1(dtnM1,fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),Ucn,Upcn,Uppcn),...
                         fuppnM1(dtnM1,fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),Ucn,Upcn,Uppcn))*...
            DUccDxnM1T(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn) +...
            dfTilcrdUpcTx(xTilnM1,...
                          fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),...
                          fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),...
                          fupnM1(dtnM1,fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn),Ucn,Upcn,Uppcn))*...
            dUpcnM1dUcnM1T(dtnM1)*...
            DUccDxnM1T(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn) +...
            dfTilcrdUppcTx(xTilnM1,...
                           fUccx(tnM1,xTilnM1,fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn),Zdn))*...
            dUppcnM1dUcnM1T(dtnM1)*...
            DUccDxnM1T(dtnM1,tnM1,xTilnM1,xTiln,xTilpn,xTilppn,Zdn));
        
        dxdzT = Azx.';
        dfcLmbdznM1T = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn)(...
            [dfTilxdxnM1T(dtnM1,tnM1,Azx.'*znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn)*dxdzT;
            zeros(rEstr,Ng+rEstr)]);
        dUccLmbdznM1T = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Zdn)(...
            DUccDxnM1T(dtnM1,tnM1,Azx.'*znM1,xTiln,xTilpn,xTilppn,Zdn)*dxdzT);
        
        %Vetor GamaFiBdLmb:
        dGamadfipBd = repositorio.dGamadfipBd;
        dGamadfipBdx = @(xTilp)(dGamadfipBd(AxFiBd.'*xTilp));
        dfipBddxpT = AxFiBd.';
        dGamadxTilpT = @(xTilp)(dGamadfipBdx(xTilp)*dfipBddxpT);
        dGamadxnM1T = @(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn)(...
            dGamadxTilpT(fupnM1(dtnM1,xTilnM1,xTiln,xTilpn,xTilppn))*dxpnM1dxnM1T(dtnM1));
        dGamaLmbdznM1T = @(dtnM1,znM1,xTiln,xTilpn,xTilppn)(...
            dGamadxnM1T(dtnM1,Azx.'*znM1,xTiln,xTilpn,xTilppn)*dxdzT);
        
        %Vetores dKcdzT e dKdcdzT:
        dKcddmid = repositorio.dKcddmid;
        dKcddmidx = @(t,xTil,Zd0)(dKcddmid(t,AxMi.'*xTil,Zd0));
        dKcddxT = @(t,xTil,Zd0)(dKcddmidx(t,xTil,Zd0)*dmiddxT);
        dKcddxnM1T = @(tnM1,xTilnM1,Zdn)(dKcddxT(tnM1,xTilnM1,Zdn));
        dKcdLmbdznM1T = @(tnM1,znM1,Zdn)(dKcddxnM1T(tnM1,Azx.'*znM1,Zdn)*dxdzT);
        
        %Vetor dmidpdzT:
        dmipddznM1T = @(dtnM1)(dmipddxpT*dxpnM1dxnM1T(dtnM1)*dxdzT);
        
        %Expressão da dinâmica interna:
        fcdMidpdmid = repositorio.fcdMidpdmid;
        fcdMidpdmidx = @(t,xTil,Zd0)(fcdMidpdmid(t,AxMi.'*xTil,Zd0));
        fcdMidpdxTilT = @(t,xTil,Zd0)(fcdMidpdmidx(t,xTil,Zd0)*dmiddxT);
        fcdMidpdxTilnM1T = @(tnM1,xTilnM1,Zdn)(fcdMidpdxTilT(tnM1,xTilnM1,Zdn));
        fcdMidpLmbdznM1T = @(tnM1,znM1,Zdn)(fcdMidpdxTilnM1T(tnM1,Azx.'*znM1,Zdn)*dxdzT);
        repositorio.fcdMidpLmbdznM1T = fcdMidpLmbdznM1T;
        fdcdMidpdmid = repositorio.fdcdMidpdmid;
        fdcdMidpdmidx = @(t,xTil,Zd0)(fdcdMidpdmid(t,AxMi.'*xTil,Zd0));
        fdcdMidpdxTilT = @(t,xTil,Zd0)(fdcdMidpdmidx(t,xTil,Zd0)*dmiddxT);
        fdcdMidpdxTilnM1T = @(tnM1,xTilnM1,Zdn)(fdcdMidpdxTilT(tnM1,xTilnM1,Zdn));
        fdcdMidpLmbdznM1T = @(tnM1,znM1,Zdn)(fdcdMidpdxTilnM1T(tnM1,Azx.'*znM1,Zdn)*dxdzT);
        repositorio.fdcdMidpLmbdznM1T = fdcdMidpLmbdznM1T;
        
        %SISTEMA EM znM1:
        dPsicdznM1T = @(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn)(...
            [KnmkLmb(dtnM1) + dfcLmbdznM1T(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Ucn,Upcn,Uppcn,Zdn) + Blmb(dtnM1)*dUccLmbdznM1T(dtnM1,tnM1,znM1,xTiln,xTilpn,xTilppn,Zdn);
            dGamaLmbdznM1T(dtnM1,znM1,xTiln,xTilpn,xTilppn) + KfiLmb - dKcdLmbdznM1T(tnM1,znM1,Zdn);
            dmipddznM1T(dtnM1) - fcdMidpLmbdznM1T(tnM1,znM1,Zdn)]);
        repositorio.dPsicdznM1T = dPsicdznM1T;
    case 12
        AxMi = 1;
        argumentos.AxMi = AxMi;
        %Vetor midpEta:
        fmidpx = @(xTilp)(AxMi.'*xTilp);
        fmidpEta = @(dtnM1,eta,xTiln,xTilpn,xTilppn)(...
            fmidpx(fupnM1(dtnM1,eta,xTiln,xTilpn,xTilppn)));
        repositorio.fmidpEta = fmidpEta;
        %Vetores fcZdHatEta e fdcZdHatEta:
        fcZdHat = repositorio.fcZdHat;
        fcZdHatx = @(t,xTil,Z0)(fcZdHat(t,AxMi.'*xTil,Z0));
        repositorio.fcZdHatx = fcZdHatx;
        fcZdHatEta = @(tnM1,eta,Z0)(fcZdHatx(tnM1,eta,Z0));
        repositorio.fcZdHatEta = fcZdHatEta;
        fdcZdHat = repositorio.fdcZdHat;
        fdcZdHatx = @(t,xTil,Z0)(fdcZdHat(t,AxMi.'*xTil,Z0));
        fdcZdHatEta = @(tnM1,eta,Z0)(fdcZdHatx(tnM1,eta,Z0));
        repositorio.fdcZdHatEta = fdcZdHatEta;
        %Express�es Uccd e Ucdcd:
        UccMid = repositorio.UccMid;
        UccMidx = @(t,xTil,Z0)(UccMid(t,AxMi.'*xTil,Z0));
        UccMidEta = @(tnM1,eta,Z0)(UccMidx(tnM1,eta,Z0));
        repositorio.UccMidEta = UccMidEta;
        UcdcMid = repositorio.UcdcMid;
        UcdcMidx = @(t,xTil,Z0)(UcdcMid(t,AxMi.'*xTil,Z0));
        UcdcMidEta = @(tnM1,eta,Z0)(UcdcMidx(tnM1,eta,Z0));
        repositorio.UcdcMidEta = UcdcMidEta;
        %Matriz dGamaEtadEtaT:
        dxTilpnM1dEtaT = @(dtnM1)(fdupduTnmk(dtnM1));
        %Matrizes dKcddEtaT e dKdcddEtaT:
        Ad = argumentos.Ad;
        dXdmi = Ad;
        dmidEtaT = AxMi.';
        %Matriz dmidpdEtaT:
        fdmidpdEtaT = @(dtnM1)(AxMi.'*dxTilpnM1dEtaT(dtnM1));
        repositorio.fdmidpdEtaT = fdmidpdEtaT;
        %Matrizes dfcZdHatEtadEtaT e dfdcZdHatdEtaT:
        dInvPsiDifcdXTHat = repositorio.dInvPsiDifcdXTHat;
        fXd = repositorio.fXd;
        dInvPsiDifcdXdTHat = @(t,mid,Z0)(dInvPsiDifcdXTHat(fXd(t,mid),Z0));
        dInvPsiDifcdXdTHatx = @(t,xTil,Z0)(dInvPsiDifcdXdTHat(t,AxMi.'*xTil,Z0));
        dInvPsiDifcdXdTHatEta = @(tnM1,eta,Z0)(dInvPsiDifcdXdTHatx(tnM1,eta,Z0));
        dfcZdHatEtadEtaT = @(tnM1,eta,Z0)(dInvPsiDifcdXdTHatEta(tnM1,eta,Z0)*dXdmi*dmidEtaT);
        repositorio.dfcZdHatEtadEtaT = dfcZdHatEtadEtaT;
        dInvPsiDifdcdXTHat = repositorio.dInvPsiDifdcdXTHat;
        dInvPsiDifdcdXdTHat = @(t,mid,Z0)(dInvPsiDifdcdXTHat(fXd(t,mid),Z0));
        dInvPsiDifdcdXdTHatx = @(t,xTil,Z0)(dInvPsiDifdcdXdTHat(t,AxMi.'*xTil,Z0));
        dInvPsiDifdcdXdTHatEta = @(tnM1,eta,Z0)(dInvPsiDifdcdXdTHatx(tnM1,eta,Z0));
        dfdcZdHatEtadEtaT = @(tnM1,eta,Z0)(dInvPsiDifdcdXdTHatEta(tnM1,eta,Z0)*dXdmi*dmidEtaT);
        repositorio.dfdcZdHatEtadEtaT = dfdcZdHatEtadEtaT;
        %Matrizes dUccMidEtadEtaT e dUcdcMidEtadEtaT:
        dUccMiddmid = repositorio.dUccMiddmid;
        dUccMiddmid = @(t,xTil,Z0)(dUccMiddmid(t,AxMi.'*xTil,Z0));
        dUccMidEtadEtaT = @(tnM1,eta,Z0)(dUccMiddmid(tnM1,eta,Z0)*dmidEtaT);
        repositorio.dUccMidEtadEtaT = dUccMidEtadEtaT;
        dUcdcMiddmid = repositorio.dUcdcMiddmid;
        dUcdcMiddmidx = @(t,xTil,Z0)(dUcdcMiddmid(t,AxMi.'*xTil,Z0));
        dUcdcMidEtadEtaT = @(tnM1,eta,Z0)(dUcdcMiddmidx(tnM1,eta,Z0));
        repositorio.dUcdcMidEtadEtaT = dUcdcMidEtadEtaT;
    otherwise
        disp('other value')
end
end