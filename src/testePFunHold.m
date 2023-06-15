clear
clc

global iter
iter = 0;
x0 = 100;
[ y0 ] = testfCustHold( x0 );
%{
figure
plot(iter,y0,'ko-')
xlabel('iterações')
ylabel('custo')
grid 
%}
hold on
iter = iter + 1;
fun = @(x)(testfCustHold( x ));
[xOtm,yOtm] = fminsearch(fun,x0);
hold off

ad = 0:1:10;
figure
plot(ad,ad.^2)
grid


