clear
clc

fRmin = 0.6; %Hz
wRmin = 2*pi*fRmin;
fRmax = 119.4; %Hz
wRmax = 2*pi*fRmax;

alfaW = 2;
%Frequência crítica (rad/s)
wAst0 = sqrt(wRmin*wRmax);
wAst = alfaW*wAst0;

deltaXiMax0 = wAst0/wRmax;
deltaXiMin0 = wRmin/wAst0;
deltaXiMax = wAst/wRmax;
deltaXiMin = wRmin/wAst;

deltaXi = 10^-2:10^-3:1.5;
wMin0 = wAst0*deltaXi;
wMax0 = wAst0./deltaXi;
wMin = wAst*deltaXi;
wMax = wAst./deltaXi;

epsADx = 10^-1;
epsADy = 10^2;
figure
p = plot(deltaXi,wMin,'k-',deltaXi,wMax,'k--',...
    deltaXiMax,wRmax,'rx',deltaXiMax,wAst*deltaXiMax,'rx',...
    deltaXiMin,wRmin,'bo',deltaXiMin,wAst/deltaXiMin,'bo',...
    deltaXi,wMin0,'k-',deltaXi,wMax0,'k--',...
    deltaXiMax0,wRmax,'kx',deltaXiMax0,wAst0*deltaXiMax0,'kx',...
    deltaXiMin0,wRmin,'ko',deltaXiMin0,wAst0/deltaXiMin0,'ko',...
    [0 0],[-epsADy wAst/deltaXiMin+epsADy],'k-',[-epsADx 1+epsADx],[0 0],'k-',...
    [0 deltaXiMax],[wRmax wRmax],'k:',[0 deltaXiMax0],[wRmin wRmin],'k:');
p(3).LineWidth = 2.5;
p(5).LineWidth = 2.5;
p(7).LineWidth = 1.5;
p(8).LineWidth = 1.5;
p(9).LineWidth = 2.5;
p(11).LineWidth = 1.5;
p(15).LineWidth = 1.5;
p(16).LineWidth = 1.5;
xlim([-epsADx 1+epsADx])
ylim([-epsADy wAst/deltaXiMin+epsADy])
xlabel('delta_{Xi}')
ylabel('w (rad/s)')
legend('gráfico delta_{Xi} x w_{m} -- alfa_{W} = 2','gráfico delta_{Xi} x w_{M} -- alfa_{W} = 2',...
    'ponto (delta_{Xi,M}, w_{R,max}) -- alfa_{W} = 2', 'ponto (delta_{Xi,M}, w* delta_{Xi,M}) -- alfa_{W} = 2',...
    'ponto (delta_{Xi,m}, w_{R,min}) -- alfa_{W} = 2', 'ponto (delta_{Xi,m}, w*/delta_{Xi,M}) -- alfa_{W} = 2',...
    'gráfico delta_{Xi} x w_{m} -- alfa_{W} = 1','gráfico delta_{Xi} x w_{M} -- alfa_{W} = 1',...
    'ponto (delta_{Xi,M}, w_{R,max}) -- alfa_{W} = 1', 'ponto (delta_{Xi,M}, w* delta_{Xi,M}) -- alfa_{W} = 1',...
    'ponto (delta_{Xi,m}, w_{R,min}) -- alfa_{W} = 1', 'ponto (delta_{Xi,m}, w*/delta_{Xi,M}) -- alfa_{W} = 1')
grid

figure
h = semilogy(deltaXi,wMin,'k-',deltaXi,wMax,'k--',...
    deltaXiMax,wRmax,'rx',deltaXiMax,wAst*deltaXiMax,'rx',...
    deltaXiMin,wRmin,'bo',deltaXiMin,wAst/deltaXiMin,'bo',...
    deltaXi,wMin0,'k-',deltaXi,wMax0,'k--',...
    deltaXiMax0,wRmax,'kx',deltaXiMax0,wAst0*deltaXiMax0,'kx',...
    deltaXiMin0,wRmin,'ko',deltaXiMin0,wAst0/deltaXiMin0,'ko',...
    [0 0],[10^(-1) 10^5],'k-',[-epsADx 1+epsADx],[10^0 10^0],'k-',...
    [0 deltaXiMax],[wRmax wRmax],'k:',[0 deltaXiMax0],[wRmin wRmin],'k:');
h(3).LineWidth = 2.5;
h(5).LineWidth = 2.5;
h(7).LineWidth = 1.5;
h(8).LineWidth = 1.5;
h(9).LineWidth = 2.5;
h(11).LineWidth = 1.5;
h(15).LineWidth = 1.5;
h(16).LineWidth = 1.5;
xlim([-epsADx 1+epsADx])
%ylim([10^-1 10^15])
xlabel('delta_{Xi}')
ylabel('w (rad/s)')
%{
legend('gráfico delta_{Xi} x w_{m} -- alfa_{W} = 2','gráfico delta_{Xi} x w_{M} -- alfa_{W} = 2',...
    'ponto (delta_{Xi,M}, w_{R,max}) -- alfa_{W} = 2', 'ponto (delta_{Xi,M}, w* delta_{Xi,M}) -- alfa_{W} = 2',...
    'ponto (delta_{Xi,m}, w_{R,min}) -- alfa_{W} = 2', 'ponto (delta_{Xi,m}, w*/delta_{Xi,M}) -- alfa_{W} = 2',...
    'gráfico delta_{Xi} x w_{m} -- alfa_{W} = 1','gráfico delta_{Xi} x w_{M} -- alfa_{W} = 1',...
    'ponto (delta_{Xi,M}, w_{R,max}) -- alfa_{W} = 1', 'ponto (delta_{Xi,M}, w* delta_{Xi,M}) -- alfa_{W} = 1',...
    'ponto (delta_{Xi,m}, w_{R,min}) -- alfa_{W} = 1', 'ponto (delta_{Xi,m}, w*/delta_{Xi,M}) -- alfa_{W} = 1')
%}
grid

figure
p = plot(deltaXi,wMin0,'k-',deltaXi,wMax0,'k--',...
    1,wAst0,'o',...
    [0 0],[-epsADy wAst/deltaXiMin+epsADy],'k-',[-epsADx 1+epsADx],[0 0],'k-',...
    0,wAst0,'kx',[0 1],[wAst0 wAst0],'k:',1,0,'kx',[1 1],[0 wAst0],'k:');
p(1).LineWidth = 1.5;
p(7).LineWidth = 1.5;
p(9).LineWidth = 1.5;
xlim([-epsADx 1+epsADx])
%ylim([-epsADy wAst/deltaXiMin+epsADy])
ylim([-epsADy 500+epsADy])
xlabel('delta_{Xi}')
ylabel('w (rad/s)')
legend('wMin = w* delta_{xi}','wMax = w*/delta_{xi}','ponto (1, w*)')
grid

figure
h = semilogy(deltaXi,wMin0,'k-',deltaXi,wMax0,'k--',...
    1,wAst0,'o',...
    [0 0],[10^(-1) 10^5],'k-',[-epsADx 1+epsADx],[10^0 10^0],'k-',...
    0,wAst0,'kx',[0 1],[wAst0 wAst0],'k:',1,1,'kx',[1 1],[1 wAst0],'k:');
h(1).LineWidth = 1.5;
h(2).LineWidth = 1.5;
h(7).LineWidth = 1.5;
h(9).LineWidth = 1.5;
xlim([-epsADx 1+epsADx])
ylim([10^-1 10^4])
xlabel('delta_{Xi}')
ylabel('w (rad/s)')
legend('wMin = w* delta_{xi}','wMax = w*/delta_{xi}','ponto (1, w*)')
grid
