clear
clc

fRmin = 0.6; %Hz
wRmin = 2*pi*fRmin;
fRmax = 119.4; %Hz
wRmax = 2*pi*fRmax;
%Coeficiente máximo de amortecimento
xiMax = 0.5;

%Coeficiente mínimo de amortecimento
xiMin = xiMax*(0.5*(sqrt(wRmin/wRmax) + sqrt(wRmax/wRmin)))^(-1);

%Frequência crítica (rad/s)
wAst = sqrt(wRmin*wRmax);

w = wRmin:10^-3:wRmax;
alfa = 2*xiMax*wRmin*wRmax/(wRmin + wRmax);
beta = 2*xiMax/(wRmin + wRmax);
f = (1/2)*(alfa./w + beta*w);

xOrigem = 3*10^0;
yOrigem = 0;
xQuadroMin = 2*10^0;
xQuadroMax = 1*10^3;
figure
p = semilogx(wRmin,yOrigem,'kx',wAst,yOrigem,'k*',wRmax,yOrigem,'ko',w,f,'k',...
    xOrigem,xiMin,'k*',xOrigem,xiMax,'kx',xOrigem,xiMax,'ko',...
    [wRmin wRmin],[yOrigem xiMax],'k:',[wRmax wRmax],[yOrigem xiMax],'k:',[wAst wAst],[yOrigem xiMin],'k:',...
    [xQuadroMin xQuadroMax],[yOrigem yOrigem],'k',[xOrigem xOrigem],[yOrigem-0.1 yOrigem+1],'k',...
    [xOrigem wRmax],[xiMax xiMax],'k:',[xOrigem wAst],[xiMin xiMin],'k:');
p(4).LineWidth = 2;
p(8).LineWidth = 2;
p(9).LineWidth = 2;
p(10).LineWidth = 2;
p(13).LineWidth = 2;
p(14).LineWidth = 2;
xlim([2*10^0 10^3])
ylim([-0.12 xiMax+0.12])
xlabel('w (rad/s)')
ylabel('xi')
legend('w_{Rmin}','w*','w_{Rmax}','Location','northeastoutside')
grid



