clear
clc

simulacao = 4;
opcao_m2 = 2;
opcao_sign = 1;
inc_dsignS = 1;
ind_LVC = 0;
[ argumentos ] = zzpreprocess( simulacao, opcao_m2, opcao_sign, inc_dsignS, ind_LVC );

mi0 = 0.01663;
mi1 = 0.7016;
mi2 = 0.7173;
beta1 = 2.0427;
beta2 = 1.9205;

miDcAtr = @(w)((mi0 + mi1*(1 - 2./(1+exp(beta1*abs(w)))) -...
    mi2*(1 - 2./(1+exp(beta2*abs(w))))).*sign(w));

dmiDcdw = @(w)((2*beta1*mi1*exp(beta1*w))/(exp(beta1*w) + 1)^2 -...
    (2*beta2*mi2*exp(beta2*w))/(exp(beta2*w) + 1)^2);
d2miDcdw2 = @(w)((2*beta1^2*mi1*exp(beta1*w))/(exp(beta1*w) + 1)^2 -...
    (4*beta1^2*mi1*exp(2*beta1*w))/(exp(beta1*w) + 1)^3 -...
    (2*beta2^2*mi2*exp(beta2*w))/(exp(beta2*w) + 1)^2 +...
    (4*beta2^2*mi2*exp(2*beta2*w))/(exp(beta2*w) + 1)^3);
Eta0 = 0.5;
Psi = @(eta)(dmiDcdw(eta));
dPsidEtaT = @(eta)(d2miDcdw2(eta));
[ EtaConv ] = f_Newton_Raphson( Eta0, Psi,...
    dPsidEtaT, argumentos );
wMax = EtaConv;
w0 = wMax;

miC0 = mi0 + mi1 - mi2;
miCmax = miDcAtr(wMax);
alfaLmb = 1.0;
lambdami = -1 + alfaLmb*miCmax/miC0;
nmi = 2;

w = -5:0.00001:5;
fiAtr = @(nu)(tanh(pi*nu) +...
    lambdami*((2*nmi)/(2*nmi - 1))*nu.*(ones(size(nu)) +...
    (1/(2*nmi - 1))*nu.^(2*nmi)).^-1);
dfiAtrdnu = @(nu)(-pi*(tanh(pi*nu).^2 - 1) +...
    (2*lambdami*nmi)./(2*nmi + nu.^(2*nmi) - 1) -...
    (4*lambdami*nmi^2*nu.^(2*nmi))./(2*nmi + nu.^(2*nmi) - 1).^2);
d2fiAtrdnu2 = @(nu)(2*pi^2*tanh(pi*nu).*(tanh(pi*nu).^2 - 1) -...
    (4*lambdami*nmi^2*nu.^(2*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^2 -...
    (8*lambdami*nmi^3*nu.^(2*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^2 +...
    (16*lambdami*nmi^3*nu.^(2*nmi).*nu.^(2*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^3);

miC = @(w)(miC0*fiAtr(w/w0));
dmiCdw = @(w)(miC0*dfiAtrdnu(w/w0)*(1/w0));
d2miCdw2 = @(w)(miC0*d2fiAtrdnu2(w/w0)*(1/w0^2));

L = 2.5; %m
d = 3.175*10^-3; %m
Jp = pi*d^4/32; %m^4
E = 200*10^9; %Pa
nuP = 0.3;
G = E/(2*(1 + nuP));

L1b = 0.026; %m
Dr1 = 0.25; %m
rho = 7850; %kg/m^3
Jp1 = pi*Dr1^4/32; %m^4
Ip1 = rho*Jp1*L1b; %Momento polar do rotor 1

k = G*Jp/L; %N/m
wN = sqrt(k/Ip1); %rad/s
qsiA = 0.01;
cA = 2*qsiA*wN;

Nr = 50; %N
R = 0.1; %m

%%
figure
plot(w,cA*w + miC(w)*Nr*R,w,cA + dmiCdw(w)*Nr*R)
grid

%%
Psi = @(eta)(cA + dmiCdw(eta)*Nr*R);
dPsidEtaT = @(eta)(d2miCdw2(eta)*Nr*R);
nRaiz = 0;
dw = 0.01*w0;
w1 = 0;
w2 = w1 + dw;
while nRaiz<2
    if Psi(w1)*Psi(w2)<0
        Eta0 = 0.5*(w1+w2);
        [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
        nRaiz = nRaiz + 1;
        if dPsidEtaT(EtaConv)>0
            wAstMin = EtaConv;
        else
            wAstMax = EtaConv;
        end
    end
    w1 = w2;
    w2 = w1 + dw;
end

Tc = cA*w;
Ts = miC(w)*Nr*R;
figure
p = plot(w,Ts,'k--',w,Tc,'k:',w,Tc+Ts,'k-',...
    [min(w) max(w)],[0 0],'k-',[0 0],[min(Tc+Ts)-0.01 max(Tc+Ts)+0.01],'k-');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
xlim([min(w) max(w)])
ylim([min(Tc+Ts)-0.015 max(Tc+Ts)+0.015])
xlabel('w (rad/s)')
ylabel('Tatr (Nm)')
legend('Atrito seco','Atrito viscoso','Atrito total','Location','northwest')
grid

Nr1 = 20; %N
Nr2 = 40; %N
Nr3 = 60; %N
Nr4 = 80; %N
Ts1 = miC(w)*Nr1*R;
Ts2 = miC(w)*Nr2*R;
Ts3 = miC(w)*Nr3*R;
Ts4 = miC(w)*Nr4*R;

%%
Nr = Nr1;
Psi = @(eta)(cA + dmiCdw(eta)*Nr*R);
dPsidEtaT = @(eta)(d2miCdw2(eta)*Nr*R);
nRaiz = 0;
dw = 0.01*w0;
w1 = 0;
w2 = w1 + dw;
while nRaiz<2
    if Psi(w1)*Psi(w2)<0
        Eta0 = 0.5*(w1+w2);
        [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
        nRaiz = nRaiz + 1;
        if dPsidEtaT(EtaConv)>0
            wAstMin1 = EtaConv;
        else
            wAstMax1 = EtaConv;
        end
    end
    w1 = w2;
    w2 = w1 + dw;
end

%%
Nr = Nr2;
Psi = @(eta)(cA + dmiCdw(eta)*Nr*R);
dPsidEtaT = @(eta)(d2miCdw2(eta)*Nr*R);
nRaiz = 0;
dw = 0.01*w0;
w1 = 0;
w2 = w1 + dw;
while nRaiz<2
    if Psi(w1)*Psi(w2)<0
        Eta0 = 0.5*(w1+w2);
        [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
        nRaiz = nRaiz + 1;
        if dPsidEtaT(EtaConv)>0
            wAstMin2 = EtaConv;
        else
            wAstMax2 = EtaConv;
        end
    end
    w1 = w2;
    w2 = w1 + dw;
end

%%
Nr = Nr3;
Psi = @(eta)(cA + dmiCdw(eta)*Nr*R);
dPsidEtaT = @(eta)(d2miCdw2(eta)*Nr*R);
nRaiz = 0;
dw = 0.01*w0;
w1 = 0;
w2 = w1 + dw;
while nRaiz<2
    if Psi(w1)*Psi(w2)<0
        Eta0 = 0.5*(w1+w2);
        [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
        nRaiz = nRaiz + 1;
        if dPsidEtaT(EtaConv)>0
            wAstMin3 = EtaConv;
        else
            wAstMax3 = EtaConv;
        end
    end
    w1 = w2;
    w2 = w1 + dw;
end

%%
Nr = Nr4;
Psi = @(eta)(cA + dmiCdw(eta)*Nr*R);
dPsidEtaT = @(eta)(d2miCdw2(eta)*Nr*R);
nRaiz = 0;
dw = 0.01*w0;
w1 = 0;
w2 = w1 + dw;
while nRaiz<2
    if Psi(w1)*Psi(w2)<0
        Eta0 = 0.5*(w1+w2);
        [ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
        nRaiz = nRaiz + 1;
        if dPsidEtaT(EtaConv)>0
            wAstMin4 = EtaConv;
        else
            wAstMax4 = EtaConv;
        end
    end
    w1 = w2;
    w2 = w1 + dw;
end

%%
figure
p = plot(w,Tc+Ts1,'k:',w,Tc+Ts2,'k-.',w,Tc+Ts3,'k--',w,Tc+Ts4,'k-',...
    [min(w) max(w)],[0 0],'k-',[0 0],[min(Tc+Ts4)-0.01 max(Tc+Ts4)+0.01],'k-',...
    wAstMin1,cA*wAstMin1 + miC(wAstMin1)*Nr1*R,'ro',...
    wAstMin2,cA*wAstMin2 + miC(wAstMin2)*Nr2*R,'ro',...
    wAstMin3,cA*wAstMin3 + miC(wAstMin3)*Nr3*R,'ro',...
    wAstMin4,cA*wAstMin4 + miC(wAstMin4)*Nr4*R,'ro',...
    wAstMin1,0,'rx',wAstMin2,0,'rx',wAstMin3,0,'rx',wAstMin4,0,'rx',...
    [wAstMin1 wAstMin1],[0 cA*wAstMin1 + miC(wAstMin1)*Nr1*R],'r:',...
    [wAstMin2 wAstMin2],[0 cA*wAstMin2 + miC(wAstMin2)*Nr2*R],'r:',...
    [wAstMin3 wAstMin3],[0 cA*wAstMin3 + miC(wAstMin3)*Nr3*R],'r:',...
    [wAstMin4 wAstMin4],[0 cA*wAstMin4 + miC(wAstMin4)*Nr4*R],'r:',...
    wAstMax1,cA*wAstMax1 + miC(wAstMax1)*Nr1*R,'bo',...
    wAstMax2,cA*wAstMax2 + miC(wAstMax2)*Nr2*R,'bo',...
    wAstMax3,cA*wAstMax3 + miC(wAstMax3)*Nr3*R,'bo',...
    wAstMax4,cA*wAstMax4 + miC(wAstMax4)*Nr4*R,'bo',...
    wAstMax1,0,'bx',wAstMax2,0,'bx',wAstMax3,0,'bx',wAstMax4,0,'bx',...
    [wAstMax1 wAstMax1],[0 cA*wAstMax1 + miC(wAstMax1)*Nr1*R],'b:',...
    [wAstMax2 wAstMax2],[0 cA*wAstMax2 + miC(wAstMax2)*Nr2*R],'b:',...
    [wAstMax3 wAstMax3],[0 cA*wAstMax3 + miC(wAstMax3)*Nr3*R],'b:',...
    [wAstMax4 wAstMax4],[0 cA*wAstMax4 + miC(wAstMax4)*Nr4*R],'b:');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
p(4).LineWidth = 2;
xlim([min(w) max(w)])
ylim([min(Tc+Ts4)-0.015 max(Tc+Ts4)+0.015])
xlabel('w (rad/s)')
ylabel('Tatr (Nm)')
legend('Nr = 20 N','Nr = 40 N','Nr = 60 N','Nr = 80 N','Location','northwest')
grid

