function [ y ] = testfCustHold( x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global iter
y = (x - 3)^2 + 5;

plot(iter,y,'ko-')
xlabel('itera��es')
ylabel('custo')
grid
iter = iter + 1;
pause(1)
end

