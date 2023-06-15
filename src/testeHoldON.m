clear
clc

%Teste 1:
%{
x = rand(50,1);
y = rand(50,1);
figure
hold on
for i = 1:1:length(x)
    i
    plot(x(i),y(i),'rx')
    xlim([-1 1])
    ylim([-1 1])
    grid
    pause(0.5)
    grid
end
grid
hold off
%}
%%
%Teste 2:
%{
t = 0;
tf = 10;
dt = 0.01;
figure
hold on
tic
while t <= tf
    t = toc;
    y = sin(pi*t);
    plot(t,y,'k.')
    xlim([0 tf])
    ylim([-1.5 1.5])
    grid
    pause(dt)
    grid
end
grid
hold off
%}
%%
%Teste 3:
%{
t = 0;
%Tamanho da janela de observação:
tJ = 2;
%Tempo de simulação:
tf = 10;
%Atraso:
dt = 0.01;

figure
hold on
tic
while t <= tf
    t = toc;
    y = sin(pi*t);
    if t<=tJ
        plot(t,y,'k.')
        xlim([0 tJ])
        ylim([-1.5 1.5])
        grid
        pause(dt)
        grid
    else
        plot(t,y,'k.')
        xlim([t-tJ t])
        ylim([-1.5 1.5])
        grid
        pause(dt)
        grid
    end
end
grid
hold off
%}
%%
%Teste 4:
%{
t = 0;
%Tamanho da janela de observação:
tJ = 5;
%Parcela da janela a partir da qual a mesma se movimenta (0<nJ<1):
nJ = 3/4;
%Tempo de simulação:
tf = 10;
%Atraso:
dt = 0.01;

figure
hold on
tic
while t <= tf
    t = toc;
    y1 = sin(pi*t); y2 = sin(pi*t - 2*pi/3); y3 = sin(pi*t - 4*pi/3);
    if t<=tJ*nJ
        plot(t,y1,'k.',t,y2,'r.',t,y3,'b.')
        xlim([0 tJ])
        ylim([-1.5 1.5])
        %legend('y1','y2','y3')
        grid
        pause(dt)
        grid
    else
        plot(t,y1,'k.',t,y2,'r.',t,y3,'b.')
        xlim([t-tJ+(tJ-tJ*nJ) t+(tJ-tJ*nJ)])
        ylim([-1.5 1.5])
        %legend('y1','y2','y3')
        grid
        pause(dt)
        grid
    end
end
grid
hold off
%}
%%
%Teste 5:
%
t = 0;
%Tamanho da janela de observação:
tJ = 5;
%Parcela da janela a partir da qual a mesma se movimenta (0<nJ<1):
nJ = 3/4;
%Tempo de simulação:
tf = 10;
%Atraso:
dt = 0.01;

figure
hold on
tic
while t <= tf
    t1 = t;
    y11 = sin(pi*t1); y21 = sin(pi*t1 - 2*pi/3); y31 = sin(pi*t1 - 4*pi/3);
    t = toc;
    y1 = sin(pi*t); y2 = sin(pi*t - 2*pi/3); y3 = sin(pi*t - 4*pi/3);
    if t<=tJ*nJ
        plot([t1 t],[y11 y1],'k-',[t1 t],[y21 y2],'r-',[t1 t],[y31 y3],'b-')
        xlim([0 tJ])
        ylim([-1.5 1.5])
        %legend('y1','y2','y3')
        grid
        pause(dt)
        grid
    else
        plot([t1 t],[y11 y1],'k-',[t1 t],[y21 y2],'r-',[t1 t],[y31 y3],'b-')
        xlim([t-tJ+(tJ-tJ*nJ) t+(tJ-tJ*nJ)])
        ylim([-1.5 1.5])
        %legend('y1','y2','y3')
        grid
        pause(dt)
        grid
    end
end
grid
hold off
%}









