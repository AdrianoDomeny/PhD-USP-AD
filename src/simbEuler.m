clear
clc

syms tetaZ tetaY tetaX

Id3 = eye(3);

Tc01 = [ cos(tetaY), 0, sin(tetaY); 
                  0, 1,          0; 
        -sin(tetaY), 0, cos(tetaY)];
Tc10 = Tc01.';
Tc12 = [cos(tetaZ), -sin(tetaZ), 0; 
        sin(tetaZ),  cos(tetaZ), 0; 
                 0,           0, 1];
Tc21 = Tc12.';
Tc23 = [1,          0,           0; 
        0, cos(tetaX), -sin(tetaX); 
        0, sin(tetaX),  cos(tetaX)];
Tc32 = Tc23.';

uHat1cI = Id3(:,2);
uHat2cI = Tc01*Id3(:,3);
uHat3cI = Tc01*Tc12*Id3(:,1);

%%
syms fi teta alfa
Te01 = [1,       0,        0; 
        0, cos(fi), -sin(fi); 
        0, sin(fi),  cos(fi)];
Te10 = Te01.';
Te12 = [ cos(teta), 0, sin(teta); 
                 0, 1,         0; 
        -sin(teta), 0, cos(teta)];
Te21 = Te12.';
Te23 = [1,         0,          0; 
        0, cos(alfa), -sin(alfa); 
        0, sin(alfa),  cos(alfa)];
Te32 = Te23.';

uHat1eI = Id3(:,1);
uHat2eI = Te01*Id3(:,2);
uHat3eI = Te01*Te12*Id3(:,1);

%%
duHat1eIdFiT = jacobian(uHat1eI,[fi,teta,alfa]);
duHat2eIdFiT = jacobian(uHat2eI,[fi,teta,alfa]);
duHat3eIdFiT = jacobian(uHat3eI,[fi,teta,alfa]);

duHat1cIdTetaT = jacobian(uHat1cI,[tetaY,tetaZ,tetaX]);
duHat2cIdTetaT = jacobian(uHat2cI,[tetaY,tetaZ,tetaX]);
duHat3cIdTetaT = jacobian(uHat3cI,[tetaY,tetaZ,tetaX]);

ad = 1;




