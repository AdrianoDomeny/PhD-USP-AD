clear
clc

fg = 6*10^-3; %folga (m)
%Kip = 10^10; %N/m
Kip = 5.2500e+09; %Nm
%Obs.:
%(1) Kip ~ E*Af/Lf = 210*10^9*(50*10^-3*5*10^-3)/(10*10^-3)
lbd = 1 + 10^-3; %adm
%lbd = 2; %adm
delta0 = lbd*fg;
n = (1/pi)*(1 - fg/delta0)^-1;

%Teste 2: Satisfatório!
nu = (0:0.0000001:delta0).';
fi = Kip*(delta0 - fg)*(tanh(n*pi*((nu/delta0) - ones(size(nu)))) + ones(size(nu)));
fiRef = Kip*(nu - fg*ones(size(nu))).*double(nu>=fg);
figure
p = plot(nu*10^3,fiRef*10^(-3),'k-',nu*10^3,fi*10^(-3),'k:',...
    fg*10^3,0,'ko',delta0*10^3,0,'kx',...
    [delta0*10^3 delta0*10^3],[0 max(fi)*10^-3],'k:',[min(nu)*10^3 max(10^3)],[0 0],'k-');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
xlabel('u_{E} (mm)')
ylabel('F_{N} (kN)')
xlim([fg*(1-0.1)*10^3 fg*(1+0.05)*10^3])
ylim([min(fi*10^(-3))-0.5*10^2 max(fi*10^(-3))+0.5*10^2])
legend('modelo original','modelo suavizado','u_{E} = f_{r}','u_{E} = f_{r} + delta_{0}','Location','northwest')
grid

%{
x = -1:0.000001:1;
nS = 10;
fiSign = tanh(nS*pi*x);
figure
plot(x,fiSign,'k')
grid
%ok!
%}
%{
%Teste 1: Não satisfatório!

x = (0:0.0001:delta0).';
A = [delta0^2, delta0^3; 2*delta0, 3*delta0^2];
b = [Kip*(delta0 - fg); Kip];
aTil = A\b;
Fn = (x - fg*ones(size(x))).*double(x>=fg);
p = [x.^2, x.^3]*aTil;
figure
plot(x,p,'b',x,Fn,'r')
grid
%}