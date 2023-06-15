clear
clc

%Dinâmica interna:
syms gama71 gama81 gama91 gama101 gama111 gama121
syms gama72 gama82 gama92 gama102 gama112 gama122

Gama = [gama71, gama101; 
    gama72, gama102];
Beta = [gama81, gama91, gama111, gama121; 
    gama82, gama92, gama112, gama122];
detG = gama71*gama102 - gama72*gama101;
detBx7 = @(gamai1,gamai2)(det([gamai1, gama101; gamai2, gama102]));
detBx10 = @(gamai1,gamai2)(det([gama71, gamai1; gama72, gamai2]));
invGB = detG^-1*...
    [  detBx7(gama81,gama82),  detBx7(gama91,gama92),  detBx7(gama111,gama112),  detBx7(gama121,gama122);
      detBx10(gama81,gama82), detBx10(gama91,gama92), detBx10(gama111,gama112), detBx10(gama121,gama122)];
%
invGBref = Gama\Beta;
erro = simplify(invGBref - invGB);
ad = 1;
%ok?
%}