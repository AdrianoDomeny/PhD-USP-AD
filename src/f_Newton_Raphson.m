function [ Eta_conv ] = f_Newton_Raphson( Eta_0, Psi,...
    dPsi_dEtaT, argumentos )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
eta = Eta_0;
nEta = length(Eta_0);
cont = 0;
Nmax = argumentos.Nmax;
%
epsilon = argumentos.epsNewtRaph;
dEta = 2*(epsilon/sqrt(nEta))*ones(nEta,1);
while (sqrt(dEta'*dEta) > epsilon) && (cont < Nmax)
    A = -dPsi_dEtaT(eta);
    %cond(A)
    b = Psi(eta);
    dEta = A\b;
    eta = eta + dEta;
    cont = cont + 1;
    %ad = 1;
end
if cont>=5
    disp('cont>=5 NwtRph')
end
%}
%{
epsilon_2 = 10^-15;
eta1 = (1/sqrt(nEta))*ones(nEta,1);
eta2 = eta1;
dEta = 2*epsilon_2*eta1;
while (sqrt(dEta'*dEta)/(0.5*sqrt((eta1 + eta2).'*(eta1 + eta2))) > epsilon_2) && (cont < Nmax);
    A = -dPsi_dEtaT(eta);
    b = Psi(eta);
    dEta = A\b;
    eta1 = eta;
    eta = eta + dEta;
    eta2 = eta;
    cont = cont + 1
end
%}
Eta_conv = eta;
end