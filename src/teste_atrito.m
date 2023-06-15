clear
clc

simulacao = 4;
opcao_m2 = 2;
opcao_sign = 1;
inc_dsignS = 1;
ind_LVC = 0;
[ argumentos ] = zzpreprocess( simulacao, opcao_m2, opcao_sign, inc_dsignS, ind_LVC );
%[ argumentos, repositorio ] = zzpreprocess2v2( argumentos );

%Modelo de atrito contínuo na origem:
lambdami = argumentos.lambdami;
nmi = argumentos.nmi;

%fiAtrVt = repositorio.fiAtr_vt;
fiAtr = @(nu)(tanh(pi*nu) +...
    lambdami*((2*nmi)/(2*nmi - 1))*nu.*(ones(size(nu)) +...
    (1/(2*nmi - 1))*nu.^(2*nmi)).^-1);
nu = -5:0.00001:5;
fiC = fiAtr(nu);

dfiAtrdnu = @(nu)(-pi*(tanh(pi*nu).^2 - 1) +...
    (2*lambdami*nmi)./(2*nmi + nu.^(2*nmi) - 1) -...
    (4*lambdami*nmi^2*nu.^(2*nmi))./(2*nmi + nu.^(2*nmi) - 1).^2);
dfiCdnu = dfiAtrdnu(nu);

%d2fiAtrdnu2 = @(nu)(2*pi^2*tanh(pi*nu).*(tanh(pi*nu).^2 - ones(size(nu))) +...
%    (16*lambdami*nmi^3*abs(nu).^(4*nmi - 2).*sign(nu).^2)./(2*nmi*ones(size(nu)) + abs(nu).^(2*nmi) - ones(size(nu))).^3 -...
%    (4*lambdami*nmi^2*abs(nu).^(2*nmi - 2).*sign(nu).^2*(2*nmi - 1))./(2*nmi*ones(size(nu)) + abs(nu).^(2*nmi) - ones(size(nu))).^2);
%d2fiAtrdnu2 = @(nu)(2*pi^2*tanh(pi*nu).*(tanh(pi*nu).^2 - ones(size(nu))) -...
%    (8*lambdami*nmi^2*abs(nu).^(2*nmi - 1).*sign(nu))./(2*nmi*ones(size(nu)) + abs(nu).^(2*nmi) - ones(size(nu))).^2 +...
%    (16*lambdami*nmi^3.*nu.*abs(nu).^(4*nmi - 2))/(2*nmi*ones(size(nu)) + abs(nu).^(2*nmi) - ones(size(nu))).^3 -...
%    (4*lambdami*nmi^2.*nu.*abs(nu).^(2*nmi - 2)*(2*nmi - 1))/(2*nmi*ones(size(nu)) + abs(nu).^(2*nmi) - ones(size(nu))).^2);
d2fiAtrdnu2 = @(nu)(2*pi^2*tanh(pi*nu).*(tanh(pi*nu).^2 - 1) -...
    (4*lambdami*nmi^2*nu.^(2*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^2 -...
    (8*lambdami*nmi^3*nu.^(2*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^2 +...
    (16*lambdami*nmi^3*nu.^(2*nmi).*nu.^(2*nmi - 1))./(2*nmi + nu.^(2*nmi) - 1).^3);
d2fiCdnu2 = d2fiAtrdnu2(nu);

figure
p = plot(nu,fiC,'k-',...
    [min(nu) max(nu)],[0 0],'k-',[0 0],[min(fiC)-0.5 max(fiC)+0.5],'k-',...
    1,fiAtr(1),'ko',-1,fiAtr(-1),'ko',...
    [1 1],[0 fiAtr(1)],'k:',[0 1],[fiAtr(1) fiAtr(1)],'k:',...
    [-1 -1],[0 fiAtr(-1)],'k:',[0 -1],[fiAtr(-1) fiAtr(-1)],'k:',...
    [0 5],[1 1],'k:',[0 -5],[-1 -1],'k:',...
    0,fiAtr(1),'kx',0,fiAtr(-1),'kx');
p(1).LineWidth = 2;
p(6).LineWidth = 1.5;
p(7).LineWidth = 1.5;
p(8).LineWidth = 1.5;
p(9).LineWidth = 1.5;
p(10).LineWidth = 1.5;
p(11).LineWidth = 1.5;
xlim([min(nu) max(nu)])
ylim([min(fiC)-0.5 max(fiC)+0.5])
xlabel('nu')
ylabel('fi_{C}')
grid
%figure
%plot(nu,dfiCdnu)
%grid
figure
plot(nu,d2fiCdnu2)
grid
%figure
%plot(nu,fiC,'b',nu,dfiCdnu,'k',nu,d2fiCdnu2,'r')
%legend('fi_{C}','dfi_{C}dnu','d2fi_{C}dnu2')
%grid

%%
%Atrito seco descontínuo na origem
mi0 = 0.01663;
mi1 = 0.7016;
mi2 = 0.7173;
beta1 = 2.0427;
beta2 = 1.9205;
miDcAtr = @(w)((mi0 + mi1*(1 - 2./(1+exp(beta1*abs(w)))) -...
    mi2*(1 - 2./(1+exp(beta2*abs(w))))).*sign(w));

w = -5:0.00001:5;
miDc = miDcAtr(w);
figure
p = plot(w,miDc,'k-',...
    [min(w) max(w)],[0 0],'k-',[0 0],[min(miDc)-0.01 max(miDc)+0.01],'k-');
p(1).LineWidth = 2;
xlabel('w (rad/s)')
ylabel('mi_{DC}')
grid

%%
%Os dois modelos superpostos:
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

miC0 = mi0 + mi1 - mi2;
miCmax = miDcAtr(wMax);
alfaLmb = 1.0;
lambdami = -1 + alfaLmb*miCmax/miC0;
nmi = 2;
%fiAtr = @(nu)(tanh(pi*nu) +...
%    lambdami*((2*nmi)/(2*nmi - 1))*nu.*(ones(size(nu)) + (1/(2*nmi - 1))*nu.^(2*nmi)).^-1);
fiAtr = @(nu)(tanh(pi*nu) +...
    lambdami*((2*nmi)/(2*nmi - 1))*nu.*(ones(size(nu)) +...
    (1/(2*nmi - 1))*nu.^(2*nmi)).^-1);
w0 = wMax;
fiC = fiAtr(w/w0);
miC = miC0*fiC;
figure
p = plot(w,miDc,'k-',w,miC,'k:',...
    [min(w) max(w)],[0 0],'k-',[0 0],[min(min(miC,miDc))-0.01 max(max(miC,miDc))+0.01],'k-');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
xlim([min(w) max(w)])
ylim([min(min(miC,miDc))-0.015 max(max(miC,miDc))+0.015])
xlabel('w (rad/s)')
ylabel('mi')
legend('mi_{DC}','mi_{C}')
grid

%%
%Observação de amortecimento negativo para baixas velocidades
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

Tc = cA*w;
Nr = 50; %N
R = 0.1; %m
Ts = miC*Nr*R;
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
legend('Atrito seco','Atrito viscoso','Location','northwest')
grid

Nr1 = 20; %N
Nr2 = 40; %N
Nr3 = 60; %N
Nr4 = 80; %N
Ts2 = miC*Nr1*R;
Ts3 = miC*Nr2*R;
Ts4 = miC*Nr3*R;
Ts5 = miC*Nr4*R;
figure
p = plot(w,Tc+Ts2,'k:',w,Tc+Ts3,'k-.',w,Tc+Ts4,'k--',w,Tc+Ts5,'k-',...
    [min(w) max(w)],[0 0],'k-',[0 0],[min(Tc+Ts5)-0.01 max(Tc+Ts5)+0.01],'k-');
p(1).LineWidth = 1;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
p(4).LineWidth = 2;
xlim([min(w) max(w)])
ylim([min(Tc+Ts5)-0.015 max(Tc+Ts5)+0.015])
xlabel('w (rad/s)')
ylabel('Tatr (Nm)')
legend('Nr = 20 N','Nr = 40 N','Nr = 60 N','Nr = 80 N','Location','northwest')
grid

miC = @(w)(miC0*fiAtr(w/w0));
dmiCdw = @(w)(miC0*dfiAtrdnu(w/w0)*(1/w0));
d2miCdw2 = @(w)(miC0*d2fiAtrdnu2(w/w0)*(1/w0^2));

figure
plot(w,cA*w + miC(w)*Nr4*R,w,cA + dmiCdw(w)*Nr4*R)
grid

Eta0 = 1.3*w0;
Psi = @(eta)(cA + dmiCdw(eta)*Nr*R);
dPsidEtaT = @(eta)(d2miCdw2(eta)*Nr*R);
[ EtaConv ] = f_Newton_Raphson( Eta0, Psi, dPsidEtaT, argumentos );
wMin = EtaConv;




