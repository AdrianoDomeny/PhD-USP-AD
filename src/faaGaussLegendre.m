function [ I ] = faaGaussLegendre( f, a, b, nGL )
%Integração:
%   f: integrando
%   a: limite inferior de integração
%   b: limite superior de integração
%   nGL: número de pontos para quadratura
m = (b-a)/2;
g = @(psi) f(a + m*(psi + 1))*m;
switch nGL
    case 2
        x2 = sqrt(3)/3;  
        w2 = 1;  
        I = w2(1)*g(x2(1)) + w2(1)*g(-x2(1));
    case 3
        x3 = [0 -sqrt(3/5) sqrt(3/5)];  
        w3 = [8/9 5/9 5/9];  
        I = w3(1)*g(x3(1)) + w3(2)*g(x3(2)) + w3(2)*g(-x3(2));
    case 4
        x4 = [sqrt((3-2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7)];  
        w4 = [(18+sqrt(30))/36 (18-sqrt(30))/36];  
        I = w4(1)*g(x4(1)) + w4(1)*g(-x4(1)) ...  
           + w4(2)*g(x4(2)) + w4(2)*g(-x4(2));
    case 5
        x5 = [0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))];  
        w5 = [128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];  
        I = w5(1)*g(x5(1)) + w5(2)*g(x5(2)) + w5(2)*g(-x5(2)) ...  
             + w5(3)*g(x5(3)) + w5(3)*g(-x5(3));
    otherwise
        disp('digite 2, 3, 4 ou 5 para nGL')
end
end