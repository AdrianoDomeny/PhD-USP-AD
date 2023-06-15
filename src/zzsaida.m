function [ argumentos, repositorio ] = zzsaida( argumentos, repositorio )
simulacao = argumentos.simulacao;
N = argumentos.N;
N1 = argumentos.N1;
Ng = argumentos.Ng;
rEstr = argumentos.rEstr;
SigmaR = argumentos.SigmaR;
switch simulacao
    case 1
        %%
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        ubTilr = ZZ(:,1:Ng-rEstr);
        ubpTilr = ZZ(:,Ng-rEstr+1:2*(Ng-rEstr));
        ubppTilr = ZZ(:,2*(Ng-rEstr)+1:3*(Ng-rEstr));
        Uc = ZZ(:,3*(Ng-rEstr)+1:3*(Ng-rEstr)+2);
        ubTil = ubTilr*SigmaR.';
        ubpTil = ubpTilr*SigmaR.';
        ubppTil = ubppTilr*SigmaR.';
        Ttx = Uc(:,1);
        Fub = Uc(:,2);
        %Movimento longitudinal:
        ub1 = ubTil(:,6*1-5);
        ubN1M1 = ubTil(:,6*(N1+1)-5);
        ubNM1 = ubTil(:,6*(N+1)-5);
        ubp1 = ubpTil(:,6*1-5);
        ubpN1M1 = ubpTil(:,6*(N1+1)-5);
        ubpNM1 = ubpTil(:,6*(N+1)-5);
        ubpp1 = ubppTil(:,6*1-5);
        ubppN1M1 = ubppTil(:,6*(N1+1)-5);
        ubppNM1 = ubppTil(:,6*(N+1)-5);
        %Movimento de torção:
        tetaX1 = ubTil(:,6*1);
        tetaXN1M1 = ubTil(:,6*(N1+1));
        tetaXNM1 = ubTil(:,6*(N+1));
        tetaXp1 = ubpTil(:,6*1);
        tetaXpN1M1 = ubpTil(:,6*(N1+1));
        tetaXpNM1 = ubpTil(:,6*(N+1));
        tetaXpp1 = ubppTil(:,6*1);
        tetaXppN1M1 = ubppTil(:,6*(N1+1));
        tetaXppNM1 = ubppTil(:,6*(N+1));
        %Movimento lateral:
        vN1M1 = ubTil(:,6*(N1+1)-4);
        wN1M1 = ubTil(:,6*(N1+1)-2);
        rN1M1 = sqrt(vN1M1.^2 + wN1M1.^2);
        
        %Entradas:
        figure
        plot(T,Ttx,'k')
        xlabel('Tempo (s)')
        ylabel('Ttx (Nm)')
        grid
        figure
        plot(T,Fub,'k')
        xlabel('Tempo (s)')
        ylabel('Fub (N)')
        grid
        %Movimento longitudinal:
        figure
        plot(T,ubp1*10^3,'b',T,ubpN1M1*10^3,'k',T,ubpNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ubp (mm/s)')
        legend('ubp1','ubpN1M1','ubpNM1')
        grid
        %Movimento de torção:
        figure
        plot(T,tetaXp1*(30/pi),'b',T,tetaXpN1M1*(30/pi),'k',T,tetaXpNM1*(30/pi),'r')
        xlabel('Tempo (s)')
        ylabel('tetaXp (RPM)')
        legend('tetaXp1','tetaXpN1M1','tetaXpNM1')
        grid
        %Movimento lateral:
        figure
        plot(T,vN1M1*10^3,'b',T,wN1M1*10^3,'r',T,rN1M1*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('v_{N1+1}','w_{N1+1}','delta_{N1+1}')
        grid
        %Deslocamento longitudinal:
        figure
        plot(T,ub1*10^3,'b',T,ubN1M1*10^3,'k',T,ubNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        legend('ub1','ubN1M1','ubNM1')
        grid
        %Deslocamento de torção:
        figure
        plot(T,tetaX1,'b',T,tetaXN1M1,'k',T,tetaXNM1,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXp (rad)')
        legend('tetaXp_1','tetaXp_{N1+1}','tetaXp_{N+1}')
        grid
        figure
        plot(tetaXNM1 - tetaX1,tetaXpNM1,'k',tetaXN1M1 - tetaX1,tetaXpN1M1,'r')
        xlabel('dtetaX (rad)')
        ylabel('tetaXp (rad/s)')
        legend('(tetaXNM1 - tetaX1) X tetaXpNM1','(tetaXN1M1 - tetaX1) X tetaXpN1M1')
        grid
        %Aceleração longitudinal:
        figure
        plot(T,ubpp1*10^3,'b',T,ubppN1M1*10^3,'k',T,ubppNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        legend('ubpp1','ubppN1M1','ubppNM1')
        grid
        %Aceleração de torção:
        figure
        plot(T,tetaXpp1,'b',T,tetaXppN1M1,'k',T,tetaXppNM1,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXpp (rad)')
        legend('tetaXpp_1','tetaXpp_{N1+1}','tetaXpp_{N+1}')
        grid
        %Interações no contato de atrito:
        %(1) Força normal:
        Fbn = repositorio.FnbNM1; %FnbNM1 = @(ub_NM1)
        Normal = Fbn(ubNM1);
        figure
        plot(T,Normal,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        %(2) Torque de atrito:
        TcxNM1vt = repositorio.TcxNM1vt; %TcxNM1vt = @(tetaXp_NM1,ub_NM1)
        Tatr = TcxNM1vt(tetaXpNM1,ubNM1);
        figure
        plot(T,Tatr,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
    case 2
        %%
        T = argumentos.T;
        %Coeficientes temporais dos modos, e suas derivadas no tempo:
        ZZ = argumentos.ZZ;
        %Vetor deslocamento e suas derivadas:
        Pk = argumentos.Pk;
        sPk = size(Pk);
        npk = sPk(1,1);
        mpk = sPk(1,2);
        ZZdesloc = ZZ*...
            [Pk.', zeros(mpk,npk), zeros(mpk,npk), zeros(mpk,2); 
             zeros(mpk,npk), Pk.', zeros(mpk,npk), zeros(mpk,2); 
             zeros(mpk,npk), zeros(mpk,npk), Pk.', zeros(mpk,2); 
             zeros(2,npk), zeros(2,npk), zeros(2,npk), eye(2)];
        ubTilr = ZZdesloc(:,1:Ng-rEstr);
        ubpTilr = ZZdesloc(:,Ng-rEstr+1:2*(Ng-rEstr));
        ubppTilr = ZZdesloc(:,2*(Ng-rEstr)+1:3*(Ng-rEstr));
        Uc = ZZdesloc(:,3*(Ng-rEstr)+1:3*(Ng-rEstr)+2);
        ubTil = ubTilr*SigmaR.';
        ubpTil = ubpTilr*SigmaR.';
        ubppTil = ubppTilr*SigmaR.';
        Ttx = Uc(:,1);
        Fub = Uc(:,2);
        %Movimento longitudinal:
        ub1 = ubTil(:,6*1-5);
        ubN1M1 = ubTil(:,6*(N1+1)-5);
        ubNM1 = ubTil(:,6*(N+1)-5);
        ubp1 = ubpTil(:,6*1-5);
        ubpN1M1 = ubpTil(:,6*(N1+1)-5);
        ubpNM1 = ubpTil(:,6*(N+1)-5);
        ubpp1 = ubppTil(:,6*1-5);
        ubppN1M1 = ubppTil(:,6*(N1+1)-5);
        ubppNM1 = ubppTil(:,6*(N+1)-5);
        %Movimento de torção:
        tetaX1 = ubTil(:,6*1);
        tetaXN1M1 = ubTil(:,6*(N1+1));
        tetaXNM1 = ubTil(:,6*(N+1));
        tetaXp1 = ubpTil(:,6*1);
        tetaXpN1M1 = ubpTil(:,6*(N1+1));
        tetaXpNM1 = ubpTil(:,6*(N+1));
        tetaXpp1 = ubppTil(:,6*1);
        tetaXppN1M1 = ubppTil(:,6*(N1+1));
        tetaXppNM1 = ubppTil(:,6*(N+1));
        %Movimento lateral:
        vN1M1 = ubTil(:,6*(N1+1)-4);
        wN1M1 = ubTil(:,6*(N1+1)-2);
        rN1M1 = sqrt(vN1M1.^2 + wN1M1.^2);
        
        %Entradas:
        figure
        plot(T,Ttx,'k')
        xlabel('Tempo (s)')
        ylabel('Ttx (Nm)')
        grid
        figure
        plot(T,Fub,'k')
        xlabel('Tempo (s)')
        ylabel('Fub (N)')
        grid
        
        %Movimento longitudinal:
        figure
        plot(T,ubp1*10^3,'b',T,ubpN1M1*10^3,'k',T,ubpNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ubp (mm/s)')
        legend('ubp1','ubpN1M1','ubpNM1')
        grid
        %Movimento de torção:
        figure
        plot(T,tetaXp1*(30/pi),'b',T,tetaXpN1M1*(30/pi),'k',T,tetaXpNM1*(30/pi),'r')
        xlabel('Tempo (s)')
        ylabel('tetaXp (RPM)')
        legend('tetaXp1','tetaXpN1M1','tetaXpNM1')
        grid
        %Movimento lateral:
        figure
        plot(T,vN1M1*10^3,'b',T,wN1M1*10^3,'r',T,rN1M1*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('v_{N1+1}','w_{N1+1}','delta_{N1+1}')
        grid
        %Deslocamento longitudinal:
        figure
        plot(T,ub1*10^3,'b',T,ubN1M1*10^3,'k',T,ubNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        legend('ub1','ubN1M1','ubNM1')
        grid
        %Deslocamento de torção:
        figure
        plot(T,tetaX1,'b',T,tetaXN1M1,'k',T,tetaXNM1,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXp (rad)')
        legend('tetaXp_1','tetaXp_{N1+1}','tetaXp_{N+1}')
        grid
        figure
        plot(tetaXNM1 - tetaX1,tetaXpNM1,'k',tetaXN1M1 - tetaX1,tetaXpN1M1,'r')
        xlabel('dtetaX (rad)')
        ylabel('tetaXp (rad/s)')
        legend('(tetaXNM1 - tetaX1) X tetaXpNM1','(tetaXN1M1 - tetaX1) X tetaXpN1M1')
        grid
        %Aceleração longitudinal:
        figure
        plot(T,ubpp1*10^3,'b',T,ubppN1M1*10^3,'k',T,ubppNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        legend('ubpp1','ubppN1M1','ubppNM1')
        grid
        %Aceleração de torção:
        figure
        plot(T,tetaXpp1,'b',T,tetaXppN1M1,'k',T,tetaXppNM1,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXpp (rad)')
        legend('tetaXpp_1','tetaXpp_{N1+1}','tetaXpp_{N+1}')
        grid
        %Interações no contato de atrito:
        %(1) Força normal:
        Fbn = repositorio.FnbNM1; %FnbNM1 = @(ub_NM1)
        Normal = Fbn(ubNM1);
        figure
        plot(T,Normal,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        %(2) Torque de atrito:
        TcxNM1vt = repositorio.TcxNM1vt; %TcxNM1vt = @(tetaXp_NM1,ub_NM1)
        Tatr = TcxNM1vt(tetaXpNM1,ubNM1);
        figure
        plot(T,Tatr,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
    case 3
        %%
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        ubTilr = ZZ(:,1:Ng-rEstr);
        ubpTilr = ZZ(:,Ng-rEstr+1:2*(Ng-rEstr));
        ubppTilr = ZZ(:,2*(Ng-rEstr)+1:3*(Ng-rEstr));
        Uc = ZZ(:,3*(Ng-rEstr)+1:3*(Ng-rEstr)+2);
        ubTil = ubTilr*SigmaR.';
        ubpTil = ubpTilr*SigmaR.';
        ubppTil = ubppTilr*SigmaR.';
        Ttx = Uc(:,1);
        Fub = Uc(:,2);
        %Movimento longitudinal:
        ub1 = ubTil(:,6*1-5);
        ubN1M1 = ubTil(:,6*(N1+1)-5);
        ubNM1 = ubTil(:,6*(N+1)-5);
        ubp1 = ubpTil(:,6*1-5);
        ubpN1M1 = ubpTil(:,6*(N1+1)-5);
        ubpNM1 = ubpTil(:,6*(N+1)-5);
        ubpp1 = ubppTil(:,6*1-5);
        ubppN1M1 = ubppTil(:,6*(N1+1)-5);
        ubppNM1 = ubppTil(:,6*(N+1)-5);
        %Movimento de torção:
        tetaX1 = ubTil(:,6*1);
        tetaXN1M1 = ubTil(:,6*(N1+1));
        tetaXNM1 = ubTil(:,6*(N+1));
        tetaXp1 = ubpTil(:,6*1);
        tetaXpN1M1 = ubpTil(:,6*(N1+1));
        tetaXpNM1 = ubpTil(:,6*(N+1));
        tetaXpp1 = ubppTil(:,6*1);
        tetaXppN1M1 = ubppTil(:,6*(N1+1));
        tetaXppNM1 = ubppTil(:,6*(N+1));
        %Movimento lateral:
        vN1M1 = ubTil(:,6*(N1+1)-4);
        wN1M1 = ubTil(:,6*(N1+1)-2);
        rN1M1 = sqrt(vN1M1.^2 + wN1M1.^2);
        
        %Entradas:
        figure
        plot(T,Ttx,'k')
        xlabel('Tempo (s)')
        ylabel('Ttx (Nm)')
        grid
        figure
        plot(T,Fub,'k')
        xlabel('Tempo (s)')
        ylabel('Fub (N)')
        grid
        %Movimento longitudinal:
        figure
        plot(T,ubp1*10^3,'b',T,ubpN1M1*10^3,'k',T,ubpNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ubp (mm/s)')
        legend('ubp1','ubpN1M1','ubpNM1')
        grid
        %Movimento de torção:
        figure
        plot(T,tetaXp1*(30/pi),'b',T,tetaXpN1M1*(30/pi),'k',T,tetaXpNM1*(30/pi),'r')
        xlabel('Tempo (s)')
        ylabel('tetaXp (RPM)')
        legend('tetaXp1','tetaXpN1M1','tetaXpNM1')
        grid
        %Movimento lateral:
        figure
        plot(T,vN1M1*10^3,'b',T,wN1M1*10^3,'r',T,rN1M1*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('v_{N1+1}','w_{N1+1}','delta_{N1+1}')
        grid
        %Deslocamento longitudinal:
        figure
        plot(T,ub1*10^3,'b',T,ubN1M1*10^3,'k',T,ubNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        legend('ub1','ubN1M1','ubNM1')
        grid
        %Deslocamento de torção:
        figure
        plot(T,tetaX1,'b',T,tetaXN1M1,'k',T,tetaXNM1,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXp (rad)')
        legend('tetaXp_1','tetaXp_{N1+1}','tetaXp_{N+1}')
        grid
        figure
        plot(tetaXNM1 - tetaX1,tetaXpNM1,'k',tetaXN1M1 - tetaX1,tetaXpN1M1,'r')
        xlabel('dtetaX (rad)')
        ylabel('tetaXp (rad/s)')
        legend('(tetaXNM1 - tetaX1) X tetaXpNM1','(tetaXN1M1 - tetaX1) X tetaXpN1M1')
        grid
        %Aceleração longitudinal:
        figure
        plot(T,ubpp1*10^3,'b',T,ubppN1M1*10^3,'k',T,ubppNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        legend('ubpp1','ubppN1M1','ubppNM1')
        grid
        %Aceleração de torção:
        figure
        plot(T,tetaXpp1,'b',T,tetaXppN1M1,'k',T,tetaXppNM1,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXpp (rad)')
        legend('tetaXpp_1','tetaXpp_{N1+1}','tetaXpp_{N+1}')
        grid
        %Interações no contato de atrito:
        %(1) Força normal:
        Fbn = repositorio.FnbNM1; %FnbNM1 = @(ub_NM1)
        Normal = Fbn(ubNM1);
        figure
        plot(T,Normal,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        %(2) Torque de atrito:
        TcxNM1vt = repositorio.TcxNM1vt; %TcxNM1vt = @(tetaXp_NM1,ub_NM1)
        Tatr = TcxNM1vt(tetaXpNM1,ubNM1);
        figure
        plot(T,Tatr,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
        
        %Rubbing:
        yL = vN1M1;
        zL = wN1M1;
        uL = rN1M1;
        figure
        plot(T,yL*10^3,'b',T,zL*10^3,'k',T,uL*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('yL','zL','uL')
        grid
        fg = argumentos.fg;
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(yL*10^3,zL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
    case 4
        %%
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        ubTil = ZZ(:,1:Ng);
        ubpTil = ZZ(:,Ng+1:2*Ng);
        ubppTil = ZZ(:,2*Ng+1:3*Ng);
        Uc = ZZ(:,3*Ng+1:3*Ng+2);
        lambdaR = ZZ(:,3*Ng+3:3*Ng+3+rEstr);
        %Entradas:
        Ttx = Uc(:,1);
        Fub = Uc(:,2);
        %Reações:
        lambdaRv1 = lambdaR(:,1);
        lambdaRtz1 = lambdaR(:,2);
        lambdaRw1 = lambdaR(:,3);
        lambdaRty1 = lambdaR(:,4);
        lambdaRvNM1 = lambdaR(:,5);
        lambdaRtzNM1 = lambdaR(:,6);
        lambdaRwNM1 = lambdaR(:,7);
        lambdaRtyNM1 = lambdaR(:,8);
        %Movimento longitudinal:
        ub1 = ubTil(:,6*1-5);
        ubN1M1 = ubTil(:,6*(N1+1)-5);
        ubNM1 = ubTil(:,6*(N+1)-5);
        ubp1 = ubpTil(:,6*1-5);
        ubpN1M1 = ubpTil(:,6*(N1+1)-5);
        ubpNM1 = ubpTil(:,6*(N+1)-5);
        ubpp1 = ubppTil(:,6*1-5);
        ubppN1M1 = ubppTil(:,6*(N1+1)-5);
        ubppNM1 = ubppTil(:,6*(N+1)-5);
        %Movimento de torção:
        tetaX1 = ubTil(:,6*1);
        tetaXN1M1 = ubTil(:,6*(N1+1));
        tetaXNM1 = ubTil(:,6*(N+1));
        tetaXp1 = ubpTil(:,6*1);
        tetaXpN1M1 = ubpTil(:,6*(N1+1));
        tetaXpNM1 = ubpTil(:,6*(N+1));
        tetaXpp1 = ubppTil(:,6*1);
        tetaXppN1M1 = ubppTil(:,6*(N1+1));
        tetaXppNM1 = ubppTil(:,6*(N+1));
        %Movimento lateral:
        vN1M1 = ubTil(:,6*(N1+1)-4);
        wN1M1 = ubTil(:,6*(N1+1)-2);
        rN1M1 = sqrt(vN1M1.^2 + wN1M1.^2);
        
        %Entradas:
        figure
        plot(T,Ttx,'k')
        xlabel('Tempo (s)')
        ylabel('Ttx (Nm)')
        grid
        figure
        plot(T,Fub,'k')
        xlabel('Tempo (s)')
        ylabel('Fub (N)')
        grid
        %Movimento longitudinal:
        figure
        plot(T,ubp1*10^3,'b',T,ubpN1M1*10^3,'k',T,ubpNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ubp (mm/s)')
        legend('ubp1','ubpN1M1','ubpNM1')
        grid
        %Movimento de torção:
        figure
        plot(T,tetaXp1*(30/pi),'b',T,tetaXpN1M1*(30/pi),'k',T,tetaXpNM1*(30/pi),'r')
        xlabel('Tempo (s)')
        ylabel('tetaXp (RPM)')
        legend('tetaXp1','tetaXpN1M1','tetaXpNM1')
        grid
        %Movimento lateral:
        figure
        plot(T,vN1M1*10^3,'b',T,wN1M1*10^3,'r',T,rN1M1*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('v_{N1+1}','w_{N1+1}','delta_{N1+1}')
        grid
        %Deslocamento longitudinal:
        figure
        plot(T,ub1*10^3,'b',T,ubN1M1*10^3,'k',T,ubNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub (mm)')
        legend('ub1','ubN1M1','ubNM1')
        grid
        %Deslocamento de torção:
        figure
        plot(T,tetaX1,'b',T,tetaXN1M1,'k',T,tetaXNM1,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXp (rad)')
        legend('tetaXp_1','tetaXp_{N1+1}','tetaXp_{N+1}')
        grid
        figure
        plot(tetaXNM1 - tetaX1,tetaXpNM1,'k',tetaXN1M1 - tetaX1,tetaXpN1M1,'r')
        xlabel('dtetaX (rad)')
        ylabel('tetaXp (rad/s)')
        legend('(tetaXNM1 - tetaX1) X tetaXpNM1','(tetaXN1M1 - tetaX1) X tetaXpN1M1')
        grid
        %Aceleração longitudinal:
        figure
        plot(T,ubpp1*10^3,'b',T,ubppN1M1*10^3,'k',T,ubppNM1*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ubpp (mm/s2)')
        legend('ubpp1','ubppN1M1','ubppNM1')
        grid
        %Aceleração de torção:
        figure
        plot(T,tetaXpp1,'b',T,tetaXppN1M1,'k',T,tetaXppNM1,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXpp (rad)')
        legend('tetaXpp_1','tetaXpp_{N1+1}','tetaXpp_{N+1}')
        grid
        %Interações no contato de atrito:
        %(1) Força normal:
        Fbn = repositorio.FnbNM1; %FnbNM1 = @(ub_NM1)
        Normal = Fbn(ubNM1);
        figure
        plot(T,Normal,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        %(2) Torque de atrito:
        TcxNM1vt = repositorio.TcxNM1vt; %TcxNM1vt = @(tetaXp_NM1,ub_NM1)
        Tatr = TcxNM1vt(tetaXpNM1,ubNM1);
        figure
        plot(T,Tatr,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
        %Rubbing:
        yL = vN1M1;
        zL = wN1M1;
        uL = rN1M1;
        figure
        plot(T,yL*10^3,'b',T,zL*10^3,'k',T,uL*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('yL','zL','uL')
        grid
        fg = argumentos.fg;
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(yL*10^3,zL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        %Reações:
        figure
        plot(T,lambdaRv1,'b',T,lambdaRw1,'k',T,lambdaRvNM1,'r',T,lambdaRwNM1,'g')
        xlabel('Tempo (s)')
        ylabel('lambdaR (N)')
        legend('lambdaRv1','lambdaRw1','lambdaRvNM1','lambdaRwNM1')
        grid
        figure
        plot(T,lambdaRtz1,'b',T,lambdaRty1,'k',T,lambdaRtzNM1,'r',T,lambdaRtyNM1,'g')
        xlabel('Tempo (s)')
        ylabel('lambdaR (Nm)')
        legend('lambdaRtz1','lambdaRty1','lambdaRtzNM1','lambdaRtyNM1')
        grid
    case 5
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        nMid = argumentos.nMid;
        Afix = argumentos.Afix;
        %{
        Amix = argumentos.Amix;
        %}
        jPu2 = argumentos.jPu2;
        Pk = argumentos.Pk;
        
        nx = nMid+2;
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        ydTil = ZZ(:,3*nx+1:3*nx+3);
        Zd = ZZ(:,3*nx+3+1:3*nx+3+2*jPu2);
        fiBdTil = xTil*Afix;
        fipBdTil = xTilp*Afix;
        fippBdTil = xTilpp*Afix;
        %{
        midTil = xTil*Amix;
        mipdTil = xTilp*Amix;
        mippdTil = xTilpp*Amix;
        %}
        uTilrkd = Zd(:,1:jPu2);
        uTilprkd = Zd(:,jPu2+1:2*jPu2);
        uTilrd = uTilrkd*Pk.';
        uTilprd = uTilprkd*Pk.';
        uTilbd = uTilrd*SigmaR.';
        uTilbpd = uTilprd*SigmaR.';
        
        %Camada limite:
        figure
        plot(T,fiBdTil(:,1),'k',T,-fiBdTil(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('fiBd1')
        grid
        figure
        plot(T,fipBdTil(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('fipBd1')
        grid
        figure
        plot(T,fippBdTil(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('fippBd1')
        grid
        figure
        plot(T,fiBdTil(:,2),'k',T,-fiBdTil(:,2),'k')
        xlabel('Tempo (s)')
        ylabel('fiBd2')
        grid
        figure
        plot(T,fipBdTil(:,2),'k')
        xlabel('Tempo (s)')
        ylabel('fipBd2')
        grid
        figure
        plot(T,fippBdTil(:,2),'k')
        xlabel('Tempo (s)')
        ylabel('fippBd2')
        grid
        
        %Trajetória desejada:
        figure
        plot(T,ydTil(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('tetaXp1 (rad/s)')
        grid
        figure
        plot(T,ydTil(:,2)*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('u1bd (mm)')
        grid
        figure
        plot(ydTil(:,1),ydTil(:,2)*10^3,'k')
        xlabel('tetaXp1 (rad/s)')
        ylabel('u1bd (mm)')
        grid
        
        %Vetor de estados desejado:
        ub1d = uTilbd(:,6*(1)-5);
        ubN1M1d = uTilbd(:,6*(N1+1)-5);
        ubNM1d = uTilbd(:,6*(N+1)-5);
        tetaX1d = uTilbd(:,6*(1));
        tetaXN1M1d = uTilbd(:,6*(N1+1));
        tetaXNM1d = uTilbd(:,6*(N+1));
        tetaXp1d = uTilbpd(:,6*(1));
        tetaXpN1M1d = uTilbpd(:,6*(N1+1));
        tetaXpNM1d = uTilbpd(:,6*(N+1));
        figure
        plot(T,ub1d*10^3,'k',T,ubN1M1d*10^3,'b',T,ubNM1d*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub_{d} (mm)')
        legend('ub_{1,d}','ub_{N1+1,d}','ub_{N+1,d}')
        grid
        figure
        plot(T,tetaX1d,'k',T,tetaXN1M1d,'b',T,tetaXNM1d,'r')
        xlabel('Tempo (s)')
        ylabel('tetaX (rad)')
        legend('tetaX_{1,d}','tetaX_{N1+1,d}','tetaX_{N+1,d}')
        grid
        figure
        plot(T,tetaXp1d*30/pi,'k',T,tetaXpN1M1d*30/pi,'b',T,tetaXpNM1d*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXpd (RPM)')
        legend('tetaXp_{1,d}','tetaXp_{N1+1,d}','tetaXp_{N+1,d}')
        grid
    case 6
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        nMid = argumentos.nMid;
        Afix = argumentos.Afix;
        %{
        Amix = argumentos.Amix;
        %}
        jPu2 = argumentos.jPu2;
        Pk = argumentos.Pk;
        
        nx = nMid+2;
        xTil = ZZ(:,1:nx);
        xTilp = ZZ(:,nx+1:2*nx);
        xTilpp = ZZ(:,2*nx+1:3*nx);
        ydTil = ZZ(:,3*nx+1:3*nx+3);
        Zd = ZZ(:,3*nx+3+1:3*nx+3+2*jPu2);
        fiBdTil = xTil*Afix;
        fipBdTil = xTilp*Afix;
        fippBdTil = xTilpp*Afix;
        %{
        midTil = xTil*Amix;
        mipdTil = xTilp*Amix;
        mippdTil = xTilpp*Amix;
        %}
        uTilrkd = Zd(:,1:jPu2);
        uTilprkd = Zd(:,jPu2+1:2*jPu2);
        uTilrd = uTilrkd*Pk.';
        uTilprd = uTilprkd*Pk.';
        uTilbd = uTilrd*SigmaR.';
        uTilbpd = uTilprd*SigmaR.';
        
        %Camada limite:
        figure
        plot(T,fiBdTil(:,1),'k',T,-fiBdTil(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('fiBd1')
        grid
        figure
        plot(T,fipBdTil(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('fipBd1')
        grid
        figure
        plot(T,fippBdTil(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('fippBd1')
        grid
        figure
        plot(T,fiBdTil(:,2),'k',T,-fiBdTil(:,2),'k')
        xlabel('Tempo (s)')
        ylabel('fiBd2')
        grid
        figure
        plot(T,fipBdTil(:,2),'k')
        xlabel('Tempo (s)')
        ylabel('fipBd2')
        grid
        figure
        plot(T,fippBdTil(:,2),'k')
        xlabel('Tempo (s)')
        ylabel('fippBd2')
        grid
        
        %Trajetória desejada:
        figure
        plot(T,ydTil(:,1),'k')
        xlabel('Tempo (s)')
        ylabel('tetaXp1 (rad/s)')
        grid
        figure
        plot(T,ydTil(:,2)*10^3,'k')
        xlabel('Tempo (s)')
        ylabel('u1bd (mm)')
        grid
        figure
        plot(ydTil(:,1),ydTil(:,2)*10^3,'k')
        xlabel('tetaXp1 (rad/s)')
        ylabel('u1bd (mm)')
        grid
        
        %Vetor de estados desejado:
        ub1d = uTilbd(:,6*(1)-5);
        ubN1M1d = uTilbd(:,6*(N1+1)-5);
        ubNM1d = uTilbd(:,6*(N+1)-5);
        tetaX1d = uTilbd(:,6*(1));
        tetaXN1M1d = uTilbd(:,6*(N1+1));
        tetaXNM1d = uTilbd(:,6*(N+1));
        tetaXp1d = uTilbpd(:,6*(1));
        tetaXpN1M1d = uTilbpd(:,6*(N1+1));
        tetaXpNM1d = uTilbpd(:,6*(N+1));
        figure
        plot(T,ub1d*10^3,'k',T,ubN1M1d*10^3,'b',T,ubNM1d*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('ub_{d} (mm)')
        legend('ub_{1,d}','ub_{N1+1,d}','ub_{N+1,d}')
        grid
        figure
        plot(T,tetaX1d,'k',T,tetaXN1M1d,'b',T,tetaXNM1d,'r')
        xlabel('Tempo (s)')
        ylabel('tetaX (rad)')
        legend('tetaX_{1,d}','tetaX_{N1+1,d}','tetaX_{N+1,d}')
        grid
        figure
        plot(T,tetaXp1d*30/pi,'k',T,tetaXpN1M1d*30/pi,'b',T,tetaXpNM1d*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetaXpd (RPM)')
        legend('tetaXp_{1,d}','tetaXp_{N1+1,d}','tetaXp_{N+1,d}')
        grid
    case 7
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        sVt = repositorio.sVt;
        
        xTil = ZZ(:,1:8);
        xTilp = ZZ(:,9:16);
        xTilpp = ZZ(:,17:24);
        Zd = ZZ(:,25:28);
        Uc = ZZ(:,29);
        Upc = ZZ(:,30);
        tetaP = ZZ(:,31);
        X = ZZ(:,32:35);
        
        uTilb = xTil(:,1:6);
        uTilbp = xTilp(:,1:6);
        uTilbpp = xTilpp(:,1:6);
        fiBd = xTil(:,7);
        fipBd = xTilp(:,7);
        fippBd = xTilpp(:,7);
        mid = xTil(:,8);
        mipd = xTilp(:,8);
        mippd = xTilpp(:,8);
        teta1 = Uc;
        teta2 = uTilb(:,1);
        teta3 = uTilb(:,6);
        teta1p = Upc;
        teta2p = uTilbp(:,1);
        teta3p = uTilbp(:,6);
        teta2pp = uTilbpp(:,1);
        teta3pp = uTilbpp(:,6);
        
        figure
        plot(T,teta1,'b',T,teta2,'k',T,teta3,'r')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        legend('teta1','teta2','teta3')
        grid
        figure
        plot(T,teta1p*30/pi,'k',T,teta2p*30/pi,'b',T,teta3p*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        legend('teta1p','teta2p','teta3p')
        grid
        
        teta2d = Zd(:,1);
        teta3d = Zd(:,2);
        figure
        plot(T,teta2d,'k',T,teta3d,'r')
        xlabel('Tempo (s)')
        ylabel('tetad (rad)')
        legend('teta2d','teta3d')
        grid
        teta2pd = Zd(:,3); 
        teta3pd = Zd(:,4);
        figure
        plot(T,teta2p*30/pi,'k',T,teta2pd*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('teta2p','teta2pd')
        grid
        figure
        plot(T,teta3p*30/pi,'k',T,teta3pd*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('teta3p','teta3pd')
        grid
        figure
        plot(T,teta2pd*30/pi,'k',T,teta3pd*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('teta2pd','teta3pd')
        grid
        figure
        plot(T,teta2pp,'k',T,teta3pp,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (rad/s2)')
        legend('teta2pp','teta3pp')
        grid
        
        figure
        plot(T,fiBd,'b',T,fipBd,'r')
        xlabel('Tempo (s)')
        ylabel('fipBd')        
        legend('fiBd','fipBd')
        grid
        figure
        plot(T,fippBd,'k')
        xlabel('Tempo (s)')
        ylabel('fippBd1')
        grid
        
        s = sVt(T,X);
        figure
        plot(T,fiBd,'k',T,-fiBd,'k',T,s,'r')
        xlabel('Tempo (s)')
        ylabel('s, fiBd')
        grid
        figure
        plot(T,mid,'b',T,mipd,'k',T,mippd,'r')
        xlabel('Tempo (s)')
        ylabel('mid, mipd, mippd')
        legend('mid (rad)','mipd (rad/s)','mippd (rad/s2)')
        grid
        mi = X(:,4);
        figure
        plot(T,mid,'b',T,mi,'k')
        xlabel('Tempo (s)')
        ylabel('mid, mi')
        legend('mid (rad)','mi (rad)')
        grid
        
        betaS = argumentos.betaS;
        u1 = tetaP/betaS;
        u2b = uTilb(:,3);
        u3b = uTilb(:,4);
        figure
        plot(T,u1*10^6,'b',T,u2b*10^6,'k',T,u3b*10^6,'r')
        xlabel('Tempo (s)')
        ylabel('u')
        legend('u1 (10^{-6} m)','u2b (10^{-6} m)','u3b (10^{-6} m)')
        grid
        
        figure
        plot(teta3 - teta1,teta3p)
        xlabel('teta3-teta1 (rad)')
        ylabel('teta3p (rad/s)')
        grid
        Normalb = repositorio.NormalbVt;
        Tatr3c = repositorio.Tatr3cVt;
        Nb = Normalb(u3b);
        Tatr3 = Tatr3c(teta3p,u3b);
        figure
        plot(T,Nb,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        figure
        plot(T,Tatr3,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
        
        %Rubbing:
        yL = uTilb(:,2);
        zL = uTilb(:,3);
        uL = sqrt(yL.^2 + zL.^2);
        figure
        plot(T,yL*10^3,'b',T,zL*10^3,'k',T,uL*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('yL','zL','uL')
        grid
        fg = argumentos.fg;
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(yL*10^3,zL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
    case 8
        N = argumentos.N;
        N1 = argumentos.N1;
        Ng = argumentos.Ng;
        rEstr = argumentos.rEstr;
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        kNmk = argumentos.k_newmark;
        
        uTilb = ZZ(:,1:Ng);
        uTilbp = ZZ(:,Ng+1:2*Ng);
        uTilbpp = ZZ(:,2*Ng+1:3*Ng);
        lambdaTil = kNmk*ZZ(:,3*Ng+1:3*Ng+rEstr);
        Uc = ZZ(:,3*Ng+rEstr+1:3*Ng+rEstr+2);
        Upc = ZZ(:,3*Ng+rEstr+3:3*Ng+rEstr+4);
        Uppc = ZZ(:,3*Ng+rEstr+5:3*Ng+rEstr+6);
        tetaP = Uc(:,2);
        tetaPp = Upc(:,2);
        tetaPpp = Uppc(:,2);
        
        tetaX1 = uTilb(:,6*1);
        tetaXN1M1 = uTilb(:,6*(N1+1));
        tetaXNM1 = uTilb(:,6*(N+1));
        tetaX1p = uTilbp(:,6*1);
        tetaXN1M1p = uTilbp(:,6*(N1+1));
        tetaXNM1p = uTilbp(:,6*(N+1));
        tetaX1pp = uTilbpp(:,6*1);
        tetaXN1M1pp = uTilbpp(:,6*(N1+1));
        tetaXNM1pp = uTilbpp(:,6*(N+1));
        
        lambdaV1 = lambdaTil(:,1);
        lambdaTetaZ1 = lambdaTil(:,2);
        lambdaW1 = lambdaTil(:,3);
        lambdaTetaY1 = lambdaTil(:,4);
        lambdaVNM1 = lambdaTil(:,5);
        lambdaTetaZNM1 = lambdaTil(:,6);
        lambdaWNM1 = lambdaTil(:,7);
        lambdaTetaYNM1 = lambdaTil(:,8);
                
        u1b = uTilb(:,1);
        uNM1b = uTilb(:,6*(N+1)-5);
        uNM1bp = uTilbp(:,6*(N+1)-5);        
        
        %Atuação:
        figure
        plot(T,tetaX1,'b',T,tetaX1p,'k',T,tetaX1pp,'r')
        xlabel('Tempo (s)')
        ylabel('tetaX1, tetaX1p, tetaX1pp')
        legend('tetaX1 (rad)','tetaX1p (rad/s)','tetaX1pp (rad/s2)')
        grid
        figure
        plot(T,tetaP,'b',T,tetaPp,'k',T,tetaPpp,'r')
        xlabel('Tempo (s)')
        ylabel('tetaP, tetaPp, tetaPpp')
        legend('tetaP (rad)','tetaPp (rad/s)','tetaPpp (rad/s2)')
        grid
        
        %Posição angular:
        figure
        plot(T,tetaX1,'b',T,tetaXN1M1,'k',T,tetaXNM1,'r')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        legend('tetaX1','tetaXN1M1','tetaXNM1')
        grid
        %Velocidade angular:
        figure
        plot(T,tetaX1p*30/pi,'k',T,tetaXN1M1p*30/pi,'b',T,tetaXNM1p*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        legend('tetaX1p','tetaXN1M1p','tetaXNM1p')
        grid
        %Aceleração angular:
        figure
        plot(T,tetaX1pp,'k',T,tetaXN1M1pp,'b',T,tetaXNM1pp,'r')
        xlabel('Tempo (s)')
        ylabel('tetapp')
        legend('tetaX1pp (rad/s2)','tetaXN1M1pp (rad/s2)','tetaXNM1pp (rad/s2)')
        grid
        %Posição longitudinal:
        ubN1M1 = uTilb(:,6*(N1+1)-5);
        ubNM1 = uTilb(:,6*(N+1)-5);
        figure
        plot(T,u1b*10^6,'b',T,ubN1M1*10^6,'k',T,ubNM1*10^6,'r')
        xlabel('Tempo (s)')
        ylabel('ub')
        legend('u1b (10^{-6} m)','ubN1M1 (10^{-6} m)','ubNM1 (10^{-6} m)')
        grid
        %Reações relativas às restrições:
        %(1) Lineares:
        figure
        plot(T,lambdaV1,'b',T,lambdaW1,'k',T,lambdaVNM1,'r',T,lambdaWNM1,'g')
        xlabel('tempo (s)')
        ylabel('lambda (N)')
        legend('lambdaV1','lambdaW1','lambdaVNM1','lambdaWNM1')
        grid
        %(2) Angulares:
        figure
        plot(T,lambdaTetaZ1,'b',T,lambdaTetaY1,'k',T,lambdaTetaZNM1,'r',T,lambdaTetaYNM1,'g')
        xlabel('tempo (s)')
        ylabel('lambda (Nm)')
        legend('lambdaTetaZ1','lambdaTetaY1','lambdaTetaZNM1','lambdaTetaYNM1')
        grid
        
        %Outros gráficos:
        figure
        plot(tetaXNM1 - tetaX1,tetaXNM1p)
        xlabel('tetaXNM1-tetaX1 (rad)')
        ylabel('tetaXNM1p (rad/s)')
        grid
        Normalb = repositorio.F_Nb_NM1;
        Tatr3c = repositorio.TcxNM1vt;
        Nb = Normalb(uNM1b,uNM1bp);
        Tatr3 = Tatr3c(tetaXNM1p,uNM1b,uNM1bp);
        figure
        plot(T,Nb,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        figure
        plot(T,Tatr3,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
        
        %Rubbing:
        vN1M1 = uTilb(:,6*(N1+1)-4);
        yL = vN1M1;
        wN1M1 = uTilb(:,6*(N1+1)-2);
        zL = wN1M1;
        uL = sqrt(yL.^2 + zL.^2);
        figure
        plot(T,yL*10^3,'b',T,zL*10^3,'k',T,uL*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('yL','zL','uL')
        grid
        fg = argumentos.fg;
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(yL*10^3,zL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        
        %Precessão e nutação do rotor intermediário:
        tetaYN1M1 = uTilb(:,6*(N1+1)-1);
        tetaZN1M1 = uTilb(:,6*(N1+1)-3);
        tetaYN1M1p = uTilbp(:,6*(N1+1)-1);
        tetaZN1M1p = uTilbp(:,6*(N1+1)-3);        
        fi = atan2(tetaZN1M1p.*cos(tetaYN1M1),tetaYN1M1p);
        fip = tetaZN1M1p.*sin(tetaYN1M1);
        figure
        plot(T,fi,'b')
        xlabel('Tempo (s)')
        ylabel('Precessão (rad)')
        title('Precessão - posição angular')
        grid
        figure
        plot(T,fip,'r')
        xlabel('Tempo (s)')
        ylabel('Precessão (rad/s)')
        title('Precessão - velocidade angular')
        grid
        teta = acos(cos(tetaYN1M1).*cos(tetaZN1M1));
        %tetap = sqrt(tetaYN1M1p.^2 + tetaZN1M1p.^2.*(cos(tetaYN1M1)).^2);
        pFi = @(fi)(double((cos(fi)>0) | (cos(fi)<0)));
        qFi = @(fi)(double(~((cos(fi)>0) | (cos(fi)<0))));
        tetap = tetaYN1M1p./cos(fi).*pFi(fi) + tetaZN1M1p.*cos(tetaYN1M1)./sin(fi).*qFi(fi);
        figure
        plot(T,teta)
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        title('Nutação - posição angular')
        grid
        figure
        plot(T,tetap)
        xlabel('Tempo (s)')
        ylabel('tetap (rad/s)')
        title('Nutação - velocidade angular')
        grid
    case 9
        N = argumentos.N;
        N1 = argumentos.N1;
        Ng = argumentos.Ng;
        rEstr = argumentos.rEstr;
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        kNmk = argumentos.k_newmark;
        
        uTilb = ZZ(:,1:Ng);
        uTilbp = ZZ(:,Ng+1:2*Ng);
        uTilbpp = ZZ(:,2*Ng+1:3*Ng);
        lambdaTil = kNmk*ZZ(:,3*Ng+1:3*Ng+rEstr);
        Uc = ZZ(:,3*Ng+rEstr+1:3*Ng+rEstr+2);
        Upc = ZZ(:,3*Ng+rEstr+3:3*Ng+rEstr+4);
        Uppc = ZZ(:,3*Ng+rEstr+5:3*Ng+rEstr+6);
        tetaP = Uc(:,2);
        tetaPp = Upc(:,2);
        tetaPpp = Uppc(:,2);
        
        tetaX1 = uTilb(:,6*1);
        tetaXN1M1 = uTilb(:,6*(N1+1));
        tetaXNM1 = uTilb(:,6*(N+1));
        tetaX1p = uTilbp(:,6*1);
        tetaXN1M1p = uTilbp(:,6*(N1+1));
        tetaXNM1p = uTilbp(:,6*(N+1));
        tetaX1pp = uTilbpp(:,6*1);
        tetaXN1M1pp = uTilbpp(:,6*(N1+1));
        tetaXNM1pp = uTilbpp(:,6*(N+1));
        
        lambdaV1 = lambdaTil(:,1);
        lambdaTetaZ1 = lambdaTil(:,2);
        lambdaW1 = lambdaTil(:,3);
        lambdaTetaY1 = lambdaTil(:,4);
        lambdaVNM1 = lambdaTil(:,5);
        lambdaTetaZNM1 = lambdaTil(:,6);
        lambdaWNM1 = lambdaTil(:,7);
        lambdaTetaYNM1 = lambdaTil(:,8);
                
        u1b = uTilb(:,1);
        uNM1b = uTilb(:,6*(N+1)-5);
        uNM1bp = uTilbp(:,6*(N+1)-5);
        
        %Atuação:
        figure
        plot(T,tetaX1,'b',T,tetaX1p,'k')
        xlabel('Tempo (s)')
        ylabel('tetaX1, tetaX1p')
        legend('tetaX1 (rad)', 'tetaX1p (rad/s)')
        grid
        figure
        plot(T,tetaX1pp,'r')
        xlabel('Tempo (s)')
        ylabel('tetaX1pp')
        legend('tetaX1pp (rad/s2)')
        grid
        figure
        plot(T,tetaP,'b',T,tetaPp,'k')
        xlabel('Tempo (s)')
        ylabel('tetaP, tetaPp')
        legend('tetaP (rad)', 'tetaPp (rad/s)')
        grid
        figure
        plot(T,tetaPpp,'r')
        xlabel('Tempo (s)')
        ylabel('tetaPpp')
        legend('tetaPpp (rad/s2)')
        grid
        
        %Posição angular:
        figure
        plot(T,tetaX1,'b',T,tetaXN1M1,'k',T,tetaXNM1,'r')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        legend('tetaX1','tetaXN1M1','tetaXNM1')
        grid
        %Velocidade angular:
        figure
        plot(T,tetaX1p*30/pi,'k',T,tetaXN1M1p*30/pi,'b',T,tetaXNM1p*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        legend('tetaX1p','tetaXN1M1p','tetaXNM1p')
        grid
        %Aceleração angular:
        figure
        plot(T,tetaX1pp,'k')
        xlabel('Tempo (s)')
        ylabel('tetapp')
        legend('tetaX1pp (rad/s2)')
        grid
        figure
        plot(T,tetaXN1M1pp,'b',T,tetaXNM1pp,'r')
        xlabel('Tempo (s)')
        ylabel('tetapp (rad/s2)')
        legend('tetaXN1M1pp (rad/s2)','tetaXNM1pp (rad/s2)')
        grid
        %Posição longitudinal:
        ubN1M1 = uTilb(:,6*(N1+1)-5);
        ubNM1 = uTilb(:,6*(N+1)-5);
        figure
        plot(T,u1b*10^6,'b',T,ubN1M1*10^6,'k',T,ubNM1*10^6,'r')
        xlabel('Tempo (s)')
        ylabel('ub')
        legend('u1b (10^{-6} m)','ubN1M1 (10^{-6} m)','ubNM1 (10^{-6} m)')
        grid
        %Reações relativas às restrições:
        %(1) Lineares:
        figure
        plot(T,lambdaV1,'b',T,lambdaW1,'k',T,lambdaVNM1,'r',T,lambdaWNM1,'g')
        xlabel('tempo (s)')
        ylabel('lambda (N)')
        legend('lambdaV1','lambdaW1','lambdaVNM1','lambdaWNM1')
        grid
        %(2) Angulares:
        figure
        plot(T,lambdaTetaZ1,'b',T,lambdaTetaY1,'k',T,lambdaTetaZNM1,'r',T,lambdaTetaYNM1,'g')
        xlabel('tempo (s)')
        ylabel('lambda (Nm)')
        legend('lambdaTetaZ1','lambdaTetaY1','lambdaTetaZNM1','lambdaTetaYNM1')
        grid
        
        %Outros gráficos:
        figure
        plot(tetaXNM1 - tetaX1,tetaXNM1p)
        xlabel('tetaXNM1-tetaX1 (rad)')
        ylabel('tetaXNM1p (rad/s)')
        grid
        Normalb = repositorio.F_Nb_NM1;
        Tatr3c = repositorio.TcxNM1vt;
        Nb = Normalb(uNM1b,uNM1bp);
        Tatr3 = Tatr3c(tetaXNM1p,uNM1b,uNM1bp);
        figure
        plot(T,Nb,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        figure
        plot(T,Tatr3,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
        
        %Rubbing:
        vN1M1 = uTilb(:,6*(N1+1)-4);
        yL = vN1M1;
        wN1M1 = uTilb(:,6*(N1+1)-2);
        zL = wN1M1;
        uL = sqrt(yL.^2 + zL.^2);
        figure
        plot(T,yL*10^3,'b',T,zL*10^3,'k',T,uL*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('yL','zL','uL')
        grid
        fg = argumentos.fg;
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(yL*10^3,zL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
    case 10
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        sVt = repositorio.sVt;
        
        xTil = ZZ(:,1:8);
        xTilp = ZZ(:,9:16);
        xTilpp = ZZ(:,17:24);
        Zd = ZZ(:,25:28);
        Uc = ZZ(:,29);
        Upc = ZZ(:,30);
        tetaP = ZZ(:,31);
        X = ZZ(:,32:35);
        
        uTilb = xTil(:,1:6);
        uTilbp = xTilp(:,1:6);
        uTilbpp = xTilpp(:,1:6);
        fiBd = xTil(:,7);
        fipBd = xTilp(:,7);
        fippBd = xTilpp(:,7);
        mid = xTil(:,8);
        mipd = xTilp(:,8);
        mippd = xTilpp(:,8);
        teta1 = Uc;
        teta2 = uTilb(:,1);
        teta3 = uTilb(:,6);
        teta1p = Upc;
        teta2p = uTilbp(:,1);
        teta3p = uTilbp(:,6);
        teta2pp = uTilbpp(:,1);
        teta3pp = uTilbpp(:,6);
        
        figure
        plot(T,teta1,'b',T,teta2,'k',T,teta3,'r')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        legend('teta1','teta2','teta3')
        grid
        figure
        plot(T,teta1p*30/pi,'k',T,teta2p*30/pi,'b',T,teta3p*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        legend('teta1p','teta2p','teta3p')
        grid
        
        teta2d = Zd(:,1);
        teta3d = Zd(:,2);
        figure
        plot(T,teta2d,'k',T,teta3d,'r')
        xlabel('Tempo (s)')
        ylabel('tetad (rad)')
        legend('teta2d','teta3d')
        grid
        teta2pd = Zd(:,3); 
        teta3pd = Zd(:,4);
        figure
        plot(T,teta2p*30/pi,'k',T,teta2pd*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('teta2p','teta2pd')
        grid
        figure
        plot(T,teta3p*30/pi,'k',T,teta3pd*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('teta3p','teta3pd')
        grid
        figure
        plot(T,teta2pd*30/pi,'k',T,teta3pd*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('teta2pd','teta3pd')
        grid
        figure
        plot(T,teta2pp,'k',T,teta3pp,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (rad/s2)')
        legend('teta2pp','teta3pp')
        grid
        
        figure
        plot(T,fiBd,'b',T,fipBd,'r')
        xlabel('Tempo (s)')
        ylabel('fipBd')        
        legend('fiBd','fipBd')
        grid
        figure
        plot(T,fippBd,'k')
        xlabel('Tempo (s)')
        ylabel('fippBd')
        grid
        
        s = sVt(T,X);
        figure
        plot(T,fiBd,'k',T,-fiBd,'k',T,s,'r')
        xlabel('Tempo (s)')
        ylabel('s, fiBd')
        grid
        figure
        plot(T,mid,'b',T,mipd,'k',T,mippd,'r')
        xlabel('Tempo (s)')
        ylabel('mid, mipd, mippd')
        legend('mid (rad)','mipd (rad/s)','mippd (rad/s2)')
        grid
        mi = X(:,4);
        figure
        plot(T,mid,'b',T,mi,'k')
        xlabel('Tempo (s)')
        ylabel('mid, mi')
        legend('mid (rad)','mi (rad)')
        grid
        
        betaS = argumentos.betaS;
        u1 = tetaP/betaS;
        u2b = uTilb(:,3);
        u3b = uTilb(:,4);
        figure
        plot(T,u1*10^6,'b',T,u2b*10^6,'k',T,u3b*10^6,'r')
        xlabel('Tempo (s)')
        ylabel('u')
        legend('u1 (10^{-6} m)','u2b (10^{-6} m)','u3b (10^{-6} m)')
        grid
        
        figure
        plot(teta3 - teta1,teta3p)
        xlabel('teta3-teta1 (rad)')
        ylabel('teta3p (rad/s)')
        grid
        Normalb = repositorio.NormalbVt;
        Tatr3c = repositorio.Tatr3cVt;
        Nb = Normalb(u3b);
        Tatr3 = Tatr3c(teta3p,u3b);
        figure
        plot(T,Nb,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        figure
        plot(T,Tatr3,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
        
        %Rubbing:
        yL = uTilb(:,2);
        zL = uTilb(:,3);
        uL = sqrt(yL.^2 + zL.^2);
        figure
        plot(T,yL*10^3,'b',T,zL*10^3,'k',T,uL*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('yL','zL','uL')
        grid
        fg = argumentos.fg;
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(yL*10^3,zL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
    case 11
        N = argumentos.N;
        N1 = argumentos.N1;
        Ng = argumentos.Ng;
        rEstr = argumentos.rEstr;
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        kNmk = argumentos.k_newmark;
        sVt = repositorio.sVt;
        
        uTilb = ZZ(:,1:Ng);
        uTilbp = ZZ(:,Ng+1:2*Ng);
        uTilbpp = ZZ(:,2*Ng+1:3*Ng);
        fiBd = ZZ(:,3*Ng+1);
        fipBd = ZZ(:,3*Ng+2);
        fippBd = ZZ(:,3*Ng+3);
        mid = ZZ(:,3*Ng+4);
        mipd = ZZ(:,3*Ng+5);
        mippd = ZZ(:,3*Ng+6);
        lambdaTil = kNmk*ZZ(:,3*Ng+6+1:3*Ng+6+rEstr);
        Uc = ZZ(:,3*Ng+rEstr+7:3*Ng+rEstr+8);
        Upc = ZZ(:,3*Ng+rEstr+9:3*Ng+rEstr+10);
        Uppc = ZZ(:,3*Ng+rEstr+11:3*Ng+rEstr+12);
        Zd = ZZ(:,3*Ng+rEstr+13:3*Ng+rEstr+16);
        X = ZZ(:,3*Ng+rEstr+17:3*Ng+rEstr+20);
        
        tetaP = Uc(:,2);
        tetaPp = Upc(:,2);
        tetaPpp = Uppc(:,2);
        
        tetaX1 = uTilb(:,6*1);
        tetaXN1M1 = uTilb(:,6*(N1+1));
        tetaXNM1 = uTilb(:,6*(N+1));
        tetaX1p = uTilbp(:,6*1);
        tetaXN1M1p = uTilbp(:,6*(N1+1));
        tetaXNM1p = uTilbp(:,6*(N+1));
        tetaX1pp = uTilbpp(:,6*1);
        tetaXN1M1pp = uTilbpp(:,6*(N1+1));
        tetaXNM1pp = uTilbpp(:,6*(N+1));
        
        lambdaV1 = lambdaTil(:,1);
        lambdaTetaZ1 = lambdaTil(:,2);
        lambdaW1 = lambdaTil(:,3);
        lambdaTetaY1 = lambdaTil(:,4);
        lambdaVNM1 = lambdaTil(:,5);
        lambdaTetaZNM1 = lambdaTil(:,6);
        lambdaWNM1 = lambdaTil(:,7);
        lambdaTetaYNM1 = lambdaTil(:,8);
                
        u1b = uTilb(:,1);
        uNM1b = uTilb(:,6*(N+1)-5);
        uNM1bp = uTilbp(:,6*(N+1)-5);
        
        %Atuação:
        figure
        plot(T,tetaX1,'b',T,tetaX1p,'k')
        xlabel('Tempo (s)')
        ylabel('tetaX1, tetaX1p')
        legend('tetaX1 (rad)', 'tetaX1p (rad/s)')
        grid
        figure
        plot(T,tetaX1pp,'r')
        xlabel('Tempo (s)')
        ylabel('tetaX1pp')
        legend('tetaX1pp (rad/s2)')
        grid
        figure
        plot(T,tetaP,'b',T,tetaPp,'k')
        xlabel('Tempo (s)')
        ylabel('tetaP, tetaPp')
        legend('tetaP (rad)', 'tetaPp (rad/s)')
        grid
        figure
        plot(T,tetaPpp,'r')
        xlabel('Tempo (s)')
        ylabel('tetaPpp')
        legend('tetaPpp (rad/s2)')
        grid
        
        %Posição angular:
        figure
        plot(T,tetaX1,'b',T,tetaXN1M1,'k',T,tetaXNM1,'r')
        xlabel('Tempo (s)')
        ylabel('teta (rad)')
        legend('tetaX1','tetaXN1M1','tetaXNM1')
        grid
        %Velocidade angular:
        figure
        plot(T,tetaX1p*30/pi,'k',T,tetaXN1M1p*30/pi,'b',T,tetaXNM1p*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetap (RPM)')
        legend('tetaX1p','tetaXN1M1p','tetaXNM1p')
        grid
        %Aceleração angular:
        figure
        plot(T,tetaX1pp,'k')
        xlabel('Tempo (s)')
        ylabel('tetapp')
        legend('tetaX1pp (rad/s2)')
        grid
        figure
        plot(T,tetaXN1M1pp,'b',T,tetaXNM1pp,'r')
        xlabel('Tempo (s)')
        ylabel('tetapp (rad/s2)')
        legend('tetaXN1M1pp (rad/s2)','tetaXNM1pp (rad/s2)')
        grid
        
        teta2d = Zd(:,1);
        teta3d = Zd(:,2);
        figure
        plot(T,teta2d,'k',T,teta3d,'r')
        xlabel('Tempo (s)')
        ylabel('tetad (rad)')
        legend('teta2d','teta3d')
        grid
        teta2pd = Zd(:,3); 
        teta3pd = Zd(:,4);
        figure
        plot(T,tetaXN1M1p*30/pi,'k',T,teta2pd*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('tetaXN1M1p','teta2pd')
        grid
        figure
        plot(T,tetaXNM1p*30/pi,'k',T,teta3pd*30/pi,'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('tetaXNM1p','teta3pd')
        grid
        
        %Posição longitudinal:
        ubN1M1 = uTilb(:,6*(N1+1)-5);
        ubNM1 = uTilb(:,6*(N+1)-5);
        figure
        plot(T,u1b*10^6,'b',T,ubN1M1*10^6,'k',T,ubNM1*10^6,'r')
        xlabel('Tempo (s)')
        ylabel('ub')
        legend('u1b (10^{-6} m)','ubN1M1 (10^{-6} m)','ubNM1 (10^{-6} m)')
        grid
        %Reações relativas às restrições:
        %(1) Lineares:
        figure
        plot(T,lambdaV1,'b',T,lambdaW1,'k',T,lambdaVNM1,'r',T,lambdaWNM1,'g')
        xlabel('tempo (s)')
        ylabel('lambda (N)')
        legend('lambdaV1','lambdaW1','lambdaVNM1','lambdaWNM1')
        grid
        %(2) Angulares:
        figure
        plot(T,lambdaTetaZ1,'b',T,lambdaTetaY1,'k',T,lambdaTetaZNM1,'r',T,lambdaTetaYNM1,'g')
        xlabel('tempo (s)')
        ylabel('lambda (Nm)')
        legend('lambdaTetaZ1','lambdaTetaY1','lambdaTetaZNM1','lambdaTetaYNM1')
        grid
        
        %Outros gráficos:
        figure
        plot(tetaXNM1 - tetaX1,tetaXNM1p)
        xlabel('tetaXNM1-tetaX1 (rad)')
        ylabel('tetaXNM1p (rad/s)')
        grid
        Normalb = repositorio.F_Nb_NM1;
        Tatr3c = repositorio.TcxNM1vt;
        Nb = Normalb(uNM1b,uNM1bp);
        Tatr3 = Tatr3c(tetaXNM1p,uNM1b,uNM1bp);
        figure
        plot(T,Nb,'b')
        xlabel('Tempo (s)')
        ylabel('Normal (N)')
        grid
        figure
        plot(T,Tatr3,'b')
        xlabel('Tempo (s)')
        ylabel('Tatr (Nm)')
        grid
        
        %Rubbing:
        vN1M1 = uTilb(:,6*(N1+1)-4);
        yL = vN1M1;
        wN1M1 = uTilb(:,6*(N1+1)-2);
        zL = wN1M1;
        uL = sqrt(yL.^2 + zL.^2);
        figure
        plot(T,yL*10^3,'b',T,zL*10^3,'k',T,uL*10^3,'r')
        xlabel('Tempo (s)')
        ylabel('Deslocamento lateral (mm)')
        legend('yL','zL','uL')
        grid
        fg = argumentos.fg;
        teta_fg = 0:0.001:(2*pi);
        figure
        plot(yL*10^3,zL*10^3,'r',fg*10^3*cos(teta_fg),fg*10^3*sin(teta_fg),'k')
        xlabel('yL (mm)')
        ylabel('zL (mm)')
        axis equal
        legend('órbita','folga')
        grid
        
        figure
        plot(T,fiBd,'b',T,fipBd,'r')
        xlabel('Tempo (s)')
        ylabel('fipBd')        
        legend('fiBd','fipBd')
        grid
        figure
        plot(T,fippBd,'k')
        xlabel('Tempo (s)')
        ylabel('fippBd')
        grid
        
        s = sVt(T,X);
        figure
        plot(T,fiBd,'k',T,-fiBd,'k',T,s,'r')
        xlabel('Tempo (s)')
        ylabel('s, fiBd')
        grid
        figure
        plot(T,mid,'b',T,mipd,'k',T,mippd,'r')
        xlabel('Tempo (s)')
        ylabel('mid, mipd, mippd')
        legend('mid (rad)','mipd (rad/s)','mippd (rad/s2)')
        grid
        mi = X(:,4);
        figure
        plot(T,mid,'b',T,mi,'k')
        xlabel('Tempo (s)')
        ylabel('mid, mi')
        legend('mid (rad)','mi (rad)')
        grid
    case 12
        T = argumentos.T;
        ZZ = argumentos.ZZ;
        mid = ZZ(:,1);
        mipd = ZZ(:,2);
        mippd = ZZ(:,3);
        Zd = ZZ(:,4:7);
        figure
        plot(T,mid,'b',T,mipd,'k',T,mippd,'r')
        xlabel('Tempo (s)')
        ylabel('mid, mipd, mippd')
        legend('mid (rad)','mipd (rad/s)','mippd (rad/s2)')
        grid
        teta2d = Zd(:,1);
        teta3d = Zd(:,2);
        figure
        plot(T,teta2d,'k',T,teta3d,'r')
        xlabel('Tempo (s)')
        ylabel('tetad (rad)')
        legend('teta2d','teta3d')
        grid
        teta2pd = Zd(:,3);
        teta3pd = Zd(:,4);
        figure
        plot(T,teta2pd*(30/pi),'k',T,teta3pd*(30/pi),'r')
        xlabel('Tempo (s)')
        ylabel('tetapd (RPM)')
        legend('teta2pd','teta3pd')
        grid
        figure
        plot(teta3d - teta2d,teta2pd*(30/pi),'k',teta3d - teta2d,teta3pd*(30/pi),'r')
        xlabel('teta3d - teta2d (rad)')
        ylabel('tetapd (RPM)')
        legend('teta2pd','teta3pd')
        grid
    otherwise
        disp('other value')
end
end







