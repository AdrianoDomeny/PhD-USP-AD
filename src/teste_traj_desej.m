clear
clc

simulacao = 11; opcao_m2 = 1; opcao_sign = 1; inc_dsignS = 1; ind_LVC = 0;
[ argumentos ] = zzpreprocess( simulacao, opcao_m2, opcao_sign, inc_dsignS, ind_LVC );
[ argumentos, repositorio ] = zzrepositoriofv6( argumentos );
%[ argumentos, repositorio ] = zzpreprocess2v2( argumentos );
M1 = argumentos.M1rt;
g = argumentos.g;

%{
syms nu alfa1 alfa2 alfa3
fi = tanh(alfa1*nu) - alfa2*nu/(alfa3^2*nu^2 + 1);
finu0 = subs(fi,nu,0);
dfidnu = diff(fi,nu);
dfidnu_nu0 = subs(dfidnu,nu,0);
%Resultado: dfidnu_nu0 = alfa1 - alfa2
%(i) dfidnu_nu0 == 0 => alfa1 = alfa2
syms alfa beta
alfa1 = alfa;
alfa2 = alfa;
alfa3 = beta;
fi = tanh(alfa1*nu) - alfa2*nu/(alfa3^2*nu^2 + 1);
dfidnu = diff(fi,nu);
d2fidnu2 = diff(dfidnu,nu);
d2fidnu2_nu0 = simplify(subs(d2fidnu2,nu,0));
%(ii) Resultado: d2fidnu2_nu0 == 0, independente de alfa e beta
d3fidnu3 = diff(d2fidnu2,nu);
d3fidnu3_nu0 = simplify(subs(d3fidnu3,nu,0));
%Resultado: d3fidnu3_nu0 == -2*alfa^3 + 6*alfa*beta^2
%(iii) 2*alfa*(3*beta^2 - alfa^2) == 0 => beta^2 = alfa^2/3
fi = tanh(alfa*nu) - (alfa*nu)/((alfa^2*nu^2)/3 + 1);
%{
fiRef = subs(fi,beta^2,alfa^2/3);
erro = simplify(fiRef- fi);
ad = 1;
%}
fiRef = subs(fi,beta^2,alfa^2/3);
dfidnu = (6*alfa^3*nu^2)/(alfa^2*nu^2 + 3)^2 - (3*alfa)/(alfa^2*nu^2 + 3) - alfa*(tanh(alfa*nu)^2 - 1);
%{
dfidnuRef = diff(fiRef,nu);
erro = simplify(dfidnuRef - dfidnu);
ad = 1;
%}
d2fidnu2 = 2*alfa^2*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1) + (18*alfa^3*nu)/(alfa^2*nu^2 + 3)^2 - (24*alfa^5*nu^3)/(alfa^2*nu^2 + 3)^3;
%{
d2fidnu2Ref = diff(fiRef,nu,2);
erro = simplify(d2fidnu2Ref - d2fidnu2);
ad = 1;
%}
d3fidnu3 = (18*alfa^3)/(alfa^2*nu^2 + 3)^2 - 2*alfa^3*(tanh(alfa*nu)^2 - 1)^2 - 4*alfa^3*tanh(alfa*nu)^2*(tanh(alfa*nu)^2 - 1) - (144*alfa^5*nu^2)/(alfa^2*nu^2 + 3)^3 + (144*alfa^7*nu^4)/(alfa^2*nu^2 + 3)^4;
%{
d3fidnu3Ref = diff(fiRef,nu,3);
erro = simplify(d3fidnu3Ref - d3fidnu3);
ad = 1;
%}
d4fidnu4 = 16*alfa^4*tanh(alfa*nu)*(tanh(alfa*nu)^2 - 1)^2 - (360*alfa^5*nu)/(alfa^2*nu^2 + 3)^3 + 8*alfa^4*tanh(alfa*nu)^3*(tanh(alfa*nu)^2 - 1) + (1440*alfa^7*nu^3)/(alfa^2*nu^2 + 3)^4 - (1152*alfa^9*nu^5)/(alfa^2*nu^2 + 3)^5;
%{
d4fidnu4Ref = diff(fiRef,nu,4);
erro = simplify(d4fidnu4Ref - d4fidnu4);
ad = 1;
%}
%}

fiOpt = @(nu,alfa)(tanh(alfa*nu) - (alfa*nu)./((alfa^2*nu.^2)/3 + ones(size(nu))));
CUSTO = @(alfa)((fiOpt(1,alfa) - 0.9)^2);
alfa0 = 28; %fi(1) == 0.9
alfa = fminsearch(CUSTO, alfa0);

fi = @(nu)(fiOpt(nu,alfa));
nu = 0:0.0001:10;
finu = fi(nu);
figure
p = plot(nu,finu,'k-',1,fi(1),'kx',[-0.2 1.2],...
    [0 0],'k-',[0 0],[-0.1 1.0],'k-',...
    [0 1],[fi(1) fi(1)],'k:',[1 1],[0 fi(1)],'k:');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(5).LineWidth = 2;
p(6).LineWidth = 2;
xlim([-0.1 1.2])
ylim([-0.1 1.0])
xlabel('nu')
ylabel('fi_{W}')
%title('nu versus finu')
grid

%wRef = 90; %RPM
wRef = 120; %RPM
t09Ref = 5; %s
y1d = @(t)(wRef*fi(t/t09Ref));
T = 0:0.001:100;
y1d_T = y1d(T);
figure
p = plot(T,y1d_T,'k-',t09Ref,y1d(t09Ref),'kx',...
    [-0.2*t09Ref 1.2*t09Ref],[0 0],'k-',...
    [0 0],[-0.1*(y1d(t09Ref)/fi(1)) 1.0*(y1d(t09Ref)/fi(1))],'k-',...
    [0 1*t09Ref],[fi(1)*(y1d(t09Ref)/fi(1)) fi(1)*(y1d(t09Ref)/fi(1))],'k:',...
    [1*t09Ref 1*t09Ref],[0 fi(1)*(y1d(t09Ref)/fi(1))],'k:');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(5).LineWidth = 2;
p(6).LineWidth = 2;
xlim([-0.1*t09Ref 1.2*t09Ref])
ylim([-0.1*(y1d(t09Ref)/fi(1)) 1.0*(y1d(t09Ref)/fi(1))])
xlabel('t (s)')
ylabel('w (RPM)')
%title('T versus y1d_T')
grid

%%
gama_TD01 = 0.1;
%Obs.:
%(1) Foram feitas simulações com gamaTD = 0.1, 1, 2 e 5
%(2) Quanto menor o valor de gamaTD, maior a concavidade da curva, e menor a
%declividade da curva no início da trajetória

wRef = wRef*pi/30; %rad/s
fU = 1;
Nref = fU*M1*g; %N
beta0 = 10^6; %rad/s
Psi = @(beta)(beta*exp(wRef/beta) - wRef - Nref/gama_TD01 - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
while fiBeta(beta0) <= beta0/2
    beta0 = beta0/2;
end
[ beta_TD01 ] = f_Newton_Raphson( beta0, Psi, dPsidbeta, argumentos );

%%
gama_TD1 = 1;
beta0 = 10^6; %rad/s
Psi = @(beta)(beta*exp(wRef/beta) - wRef - Nref/gama_TD1 - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
while fiBeta(beta0) <= beta0/2
    beta0 = beta0/2;
end
[ beta_TD1 ] = f_Newton_Raphson( beta0, Psi, dPsidbeta, argumentos );

%%
gama_TD2 = 2;
beta0 = 10^6; %rad/s
Psi = @(beta)(beta*exp(wRef/beta) - wRef - Nref/gama_TD2 - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
while fiBeta(beta0) <= beta0/2
    beta0 = beta0/2;
end
[ beta_TD2 ] = f_Newton_Raphson( beta0, Psi, dPsidbeta, argumentos );

%%
gama_TD5 = 5;
beta0 = 10^6; %rad/s
Psi = @(beta)(beta*exp(wRef/beta) - wRef - Nref/gama_TD5 - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
while fiBeta(beta0) <= beta0/2
    beta0 = beta0/2;
end
[ beta_TD5 ] = f_Newton_Raphson( beta0, Psi, dPsidbeta, argumentos );

%%
w = 0:0.0001:1.1*wRef;
Nr01 = gama_TD01*(beta_TD01*exp(w/beta_TD01) - w - beta_TD01);
Nr1 = gama_TD1*(beta_TD1*exp(w/beta_TD1) - w - beta_TD1);
Nr2 = gama_TD2*(beta_TD2*exp(w/beta_TD2) - w - beta_TD2);
Nr5 = gama_TD5*(beta_TD5*exp(w/beta_TD5) - w - beta_TD5);

figure
p = plot(w*(30/pi),Nr01,'k-',w*(30/pi),Nr1,'k--',w*(30/pi),Nr2,'k-.',...
    w*(30/pi),Nr5,'k:',wRef*(30/pi),Nref,'ko');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;
p(3).LineWidth = 1.5;
p(4).LineWidth = 1.5;
p(5).LineWidth = 1.5;
xlabel('teta3p (RPM)')
ylabel('Nr (N)')
grid
legend('gama_{TD} = 0,1','gama_{TD} = 1','gama_{TD} = 2','gama_{TD} = 5',...
    'Ponto de operação','Location','northwest')


%%
gama_TD = 0.1;
%Obs.:
%(1) Foram feitas simulações com gamaTD = 0.1, 1, 2 e 5
%(2) Quanto menor o valor de gamaTD, maior a concavidade da curva, e menor a
%declividade da curva no início da trajetória
%Velocidade angular em RPM:
wRef = 120; %RPM
%Velocidade angular em rad/s:
wRef = wRef*pi/30; %rad/s
%Fator de carga:
fU = 1;
%Carga de operação:
%Obs.: M1 é massa do rotor 1, e g a gravidade
Nref = fU*M1*g; %N
%Função
Psi = @(beta)(beta*exp(wRef/beta) - wRef - Nref/gama_TD - beta);
dPsidbeta = @(beta)(exp(wRef/beta) - (wRef/beta)*exp(wRef/beta) - 1);
fiBeta = @(beta0)(beta0 - Psi(beta0)/dPsidbeta(beta0));
%Primeiro chute inicial:
%Obs.: É importante que seja um valor alto!!!
beta0 = 10^6; %rad/s
%Melhor chute inicial:
%Obs.: a convergência de beta_TD é muito sensível ao chute inicial beta0
while fiBeta(beta0) <= beta0/2
    beta0 = beta0/2;
end
[ beta_TD ] = f_Newton_Raphson( beta0, Psi, dPsidbeta, argumentos );

N = gama_TD*(beta_TD*exp(w/beta_TD) - w - beta_TD);
dNdw = gama_TD*(exp(w/beta_TD) - ones(size(w)));
d2Ndw2 = (gama_TD/beta_TD)*exp(w/beta_TD);
d3Ndw3 = (gama_TD/beta_TD^2)*exp(w/beta_TD);
d4Ndw4 = (gama_TD/beta_TD^3)*exp(w/beta_TD);

y1d = @(t)(wRef*fi(t/t09Ref));
y1d_T = y1d(T);
y2d_T = gama_TD*(beta_TD*exp(y1d_T/beta_TD) - y1d_T - beta_TD*ones(size(y1d_T)));
figure
plot(T,y2d_T)
title('T versus y2d_T')
grid

figure
plot(w*(30/pi),N,'k',w*(30/pi),dNdw,'r',w*(30/pi),d2Ndw2,'b',w*(30/pi),d3Ndw3,'g',w*(30/pi),d4Ndw4,'y',wRef*(30/pi),Nref,'o')
xlabel('teta3p (RPM)')
grid
legend('Nr (N)','dNrdw (N/rad*s^{-1})','d2Nrdw2 (N/rad*s^{-2})','d3Nrdw3 (N/rad*s^{-3})','d4Nrdw4 (N/rad*s^{-4})','Ponto de operação')

beta = 0.5:0.0001:3;
psi = beta.*exp(wRef*ones(size(beta))./beta) - wRef*ones(size(beta)) - Nref/gama_TD*ones(size(beta)) - beta;
psi_TD = beta_TD*exp(wRef/beta_TD) - wRef - Nref/gama_TD - beta_TD;
figure
plot(beta,psi,'b',beta_TD,psi_TD,'rx',1,0,'ko')
legend('beta versus psi','beta_{TD}','chute inicial')
title('curva beta_{TD}')
grid

%%
G = argumentos.G;
J = argumentos.J;
L = argumentos.L;
I1 = argumentos.Ip1;
mic0 = argumentos.mic0;
w0 = argumentos.w0;
R3 = argumentos.R3;
dfiAtrdnu = repositorio.dfiAtrdnu;
dfiAtrdnuVt = repositorio.dfiAtrdnuVt;
d2fiAtrdnu2 = repositorio.d2fiAtrdnu2;
d2fiAtrdnu2Vt = repositorio.d2fiAtrdnu2Vt;

fU = 1;
Nref = fU*M1*g;
K = G*J/L;
cR = 2*sqrt(K*I1);
qsi = 0.03;
c3 = qsi*cR;
Eta0 = 10*w0;
Psi = @(eta)(c3 + mic0*(1/w0)*dfiAtrdnu(eta/w0)*Nref*R3);
dPsidEta = @(eta)(mic0*(1/w0^2)*d2fiAtrdnu2(eta/w0)*Nref*R3);
[ EtaSol ] = f_Newton_Raphson( Eta0, Psi, dPsidEta, argumentos );
w3ast = EtaSol;

PsiVt = @(eta)(c3*ones(size(eta)) + mic0*(1/w0)*dfiAtrdnuVt(eta/w0)*Nref*R3);
wTeste = 0:0.0001:1;
EtaTeste = wTeste;
PsiTeste = PsiVt(EtaTeste);
figure
plot(EtaTeste,PsiTeste,'b',w3ast,Psi(w3ast),'rx',Eta0,0,'ko')
legend('wTeste versus c3Mdmicdw','w3ast','chute inicial')
%title('curva beta_{TD}')
grid

%
incr = 10^-3;
fUvt = 1:-incr:incr;
fUvt = [fUvt, incr-incr^2:-incr^2:incr^2];
nF = length(fUvt);
w3astVt = zeros(1,nF);
for i=1:nF
    i
    Ncg = fUvt(i)*M1*g;
    Psi = @(eta)(c3 + mic0*(1/w0)*dfiAtrdnu(eta/w0)*Ncg*R3);
    dPsidEta = @(eta)(mic0*(1/w0^2)*d2fiAtrdnu2(eta/w0)*Ncg*R3);
    [ EtaSol ] = f_Newton_Raphson( Eta0, Psi, dPsidEta, argumentos );
    w3astVt(1,i) = EtaSol;
    Eta0 = EtaSol;
end
fiW = atan2(1,qsi);
wRefcr1 = w3astVt*sqrt(1+qsi^2)*(1/sqrt(1+qsi^2) - exp(-qsi*(5*pi/2 - fiW)))^-1;
alfaWref = 0.65;
wMax = alfaWref*wRef;
w0 = 2*wMax;
Nts2 = Nref*(w/wRef).*((w - w0)/(wRef - w0));
figure
p = plot(wRefcr1*30/pi,fUvt*M1*g,'k--',w*(30/pi),N,'k-',...
    w*(30/pi),Nts2,'k:',...
    wRef*(30/pi),Nref,'ko',...
    [0,0],[0,220],'k',[0,140],[0,0],'k',...
    [0,1]*wRefcr1(1,10)*30/pi,[1,1]*fUvt(1,10)*M1*g,'k',...
    [0,1]*wRefcr1(1,125)*30/pi,[1,1]*fUvt(1,125)*M1*g,'k',...
    [0,1]*wRefcr1(1,250)*30/pi,[1,1]*fUvt(1,250)*M1*g,'k',...
    [0,1]*wRefcr1(1,370)*30/pi,[1,1]*fUvt(1,370)*M1*g,'k',...
    [0,1]*wRefcr1(1,500)*30/pi,[1,1]*fUvt(1,500)*M1*g,'k',...
    [0,1]*wRefcr1(1,625)*30/pi,[1,1]*fUvt(1,625)*M1*g,'k',...
    [0,1]*wRefcr1(1,750)*30/pi,[1,1]*fUvt(1,750)*M1*g,'k',...
    [0,1]*wRefcr1(1,875)*30/pi,[1,1]*fUvt(1,875)*M1*g,'k');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
p(4).LineWidth = 2;
xlabel('Velocidade (RPM)')
ylabel('Carga de operação (N)')
xlim([0 140])
ylim([0 220])
%title('Mapa de operador')
legend('Limite da zona de estabilidade - tetap_{X,1}',...
    'Trajetória 1 - tetap_{X,N+1}','Trajetória 2 - tetap_{X,N+1}',...
    'Ponto de operação (w_{OP}, N_{OP})','Location','northwest')
grid
%}

%Sistema com dois rotores:
I1 = argumentos.Ip1;
I2b = argumentos.Ip2b;
dm = argumentos.dm;
l2 = argumentos.l2;
I2 = I2b + dm*l2^2;
K1 = argumentos.K10;
K2 = argumentos.K20;

I = diag([I2,I1]);
K = [K1 + K2, -K2; -K2, K2];
b = (K1 + K2)/I2 + K2/I1;
c = K1*K2/(I1*I2);
delta = b^2 - 4*c;
lambda1 = (1/2)*(b - sqrt(delta));
lambda2 = (1/2)*(b + sqrt(delta));
w1 = sqrt(lambda1);
w2 = sqrt(lambda2);
f1 = w1/(2*pi);
f2 = w2/(2*pi);

P1 = [1 - lambda1*(I1/K2); 1];
P2 = [1 - lambda2*(I1/K2); 1];
P1hat = P1/sqrt(P1.'*I*P1);
P2hat = P2/sqrt(P2.'*I*P2);
Phat = [P1hat, P2hat];
p11 = Phat(1,1);
p21 = Phat(2,1);
p12 = Phat(1,2);
p22 = Phat(2,2);

figure
plot([1,2],P1hat,'k',[1,2],P2hat,'r')
grid
legend('P1hat','P2hat')

T = 0:0.001:50;
qsi1 = 0.1;
qsi2 = 0.001;
fiW1 = atan2(1,qsi1);
fiW2 = atan2(1,qsi2);
fW3 = @(t)(wRef*K1*p21*...
    (p22/(w2^2*(qsi2^2 + 1)) + p11/(w1^2*(qsi1^2 + 1))) -...
    wRef*K1*p21*...
    (p11*exp(-qsi1*w1*t)/(w1^2*sqrt(qsi1^2 + 1))*sin(w1*t + fiW1) +...
    p22*exp(-qsi2*w2*t)/(w2^2*sqrt(qsi2^2 + 1))*sin(w2*t + fiW2)));
fW3vt = @(T)(wRef*K1*p21*ones(size(T))*...
    (p22/(w2^2*(qsi2^2 + 1)) + p11/(w1^2*(qsi1^2 + 1))) -...
    wRef*K1*p21*...
    ((p11*exp(-qsi1*w1*T)/(w1^2*sqrt(qsi1^2 + 1))).*sin(w1*T + fiW1*ones(size(T))) +...
    (p22*exp(-qsi2*w2*T)/(w2^2*sqrt(qsi2^2 + 1))).*sin(w2*T + fiW2*ones(size(T)))));
fW3vtf1 = @(T)(wRef*K1*p21*ones(size(T))*...
    (p11/(w1^2*(qsi1^2 + 1))) -...
    wRef*K1*p21*...
    ((p11*exp(-qsi1*w1*T)/(w1^2*sqrt(qsi1^2 + 1))).*sin(w1*T + fiW1*ones(size(T)))));
fW3vtf2 = @(T)(wRef*K1*p21*ones(size(T))*...
    (p22/(w2^2*(qsi2^2 + 1))) -...
    wRef*K1*p21*...
    ((p22*exp(-qsi2*w2*T)/(w2^2*sqrt(qsi2^2 + 1))).*sin(w2*T + fiW2*ones(size(T)))));
W3 = fW3vt(T);
W3f1 = fW3vtf1(T);
W3f2 = fW3vtf2(T);
figure
plot(T,W3*30/pi,'b',T,W3f1*30/pi,'k',T,W3f2*30/pi,'r')
xlabel('tempo (s)')
ylabel('w3 (RPM)')
legend('w3','w3f1','w3f2')
grid
figure
p = plot(T,W3f1*30/pi,'r-',T,W3f2*30/pi,'b-',T,W3*30/pi,'k-');
p(1).LineWidth = 2;
xlabel('tempo (s)')
ylabel('w_{r,1} (RPM)')
legend('1^{o} modo','2^{o} modo','sinal completo')
grid
%Conclusão:
%(1) O primeiro modo é dominante
%(2) O segundo modo quase não influi na saída final

wRefcr2 = w3astVt*sqrt(1+qsi^2)/(K1*p21)*(1/sqrt(1+qsi^2)*(p11/w1^2 + p22/w2^2) - (p11/w1^2)*exp(-qsi*(5*pi/2 - fiW)))^-1;
figure
plot(wRefcr2*30/pi,fUvt*M1*g,'g',w*(30/pi),N,'b',wRef*(30/pi),Nref,'ko',[0,0],[0,1]*max(N),'k',...
    [0,1]*wRefcr2(1,10)*30/pi,[1,1]*fUvt(1,10)*M1*g,'r',...
    [0,1]*wRefcr2(1,125)*30/pi,[1,1]*fUvt(1,125)*M1*g,'r',...
    [0,1]*wRefcr2(1,250)*30/pi,[1,1]*fUvt(1,250)*M1*g,'r',...
    [0,1]*wRefcr2(1,370)*30/pi,[1,1]*fUvt(1,370)*M1*g,'r',...
    [0,1]*wRefcr2(1,500)*30/pi,[1,1]*fUvt(1,500)*M1*g,'r',...
    [0,1]*wRefcr2(1,625)*30/pi,[1,1]*fUvt(1,625)*M1*g,'r',...
    [0,1]*wRefcr2(1,750)*30/pi,[1,1]*fUvt(1,750)*M1*g,'r',...
    [0,1]*wRefcr2(1,875)*30/pi,[1,1]*fUvt(1,875)*M1*g,'r')
xlabel('Velocidade de operação (RPM)')
ylabel('Carga de operação (N)')
title('Mapa de operador')
legend('Limite da zona de estabilidade','Trajetória desejada','Ponto de operação')
grid
%}



