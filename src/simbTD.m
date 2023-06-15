clear
clc

syms nu alfa1 alfa2 alfa3
fi = tanh(alfa1*nu) - alfa2*nu/((alfa3*nu)^2 + 1);

sfi = simplify(subs(fi,nu,0));
%sfi = 0 ok!!

%{
dfidnu = diff(fi,nu);
dfidnu = simplify(dfidnu);
sdfidnu = simplify(subs(dfidnu,nu,0));
%Para que dfidnu(nu=0) = 0, então alfa1 = alfa2
%}

syms alfa
alfa1 = alfa;
alfa2 = alfa;
fi = tanh(alfa*nu) - alfa*nu/((alfa3*nu)^2 + 1);
dfidnu = diff(fi,nu);
dfidnu = simplify(dfidnu);
sdfidnu = simplify(subs(dfidnu,nu,0));
%ok!

d2fidnu2 = diff(dfidnu,nu);
d2fidnu2 = simplify(d2fidnu2);
sd2fidnu2 = simplify(subs(d2fidnu2,nu,0));
%ok!

%{
d3fidnu3 = diff(d2fidnu2,nu);
d3fidnu3 = simplify(d3fidnu3);
sd3fidnu3 = simplify(subs(d3fidnu3,nu,0));
%Para que d3fidnu3(nu=0) = 0, então - 2*alfa^3 + 6*alfa*alfa3^2 = 0, ou
%melhor, alfa3^2 = alfa^2/3
%}

%{
alfa3 = alfa/sqrt(3);
fi = tanh(alfa*nu) - alfa*nu/((alfa3*nu)^2 + 1);
dfidnu = diff(fi,nu);
d2fidnu2 = diff(dfidnu,nu);
d3fidnu3 = diff(d2fidnu2,nu);
sd3fidnu3 = simplify(subs(d3fidnu3,nu,0));
%Ok!
%}

fi = tanh(alfa*nu) - alfa*nu/((1/3)*(alfa*nu)^2 + 1);
%{
dfidnu = diff(fi,nu);
d2fidnu2 = diff(dfidnu,nu);
d3fidnu3 = diff(d2fidnu2,nu);
sd3fidnu3 = simplify(subs(d3fidnu3,nu,0));
%ok!
%}

%%
syms t w Nr0b(w) dNr0bdw(w) d2Nr0bdw2(w) d3Nr0bdw3(w) d4Nr0bdw4(w) d5Nr0bdw5(w)
syms y1d(t) y1pd(t) y1ppd(t) y1pppd(t) y1ppppd(t) y1pppppd(t)
syms y2d(t) y2pd(t) y2ppd(t) y2pppd(t) y2ppppd(t) y2pppppd(t)

y2pppppd = d5Nr0bdw5(y1d(t))*(y1pd(t))^5 +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppd(t) +...
    15*d3Nr0bdw3(y1d(t))*y1pd(t)*(y1ppd(t))^2 +...
    10*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppd(t) +...
    10*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1d(t)*y1ppppd(t) +...
    dNr0bdw(y1d(t))*y1pppppd(t);

y2pppppdRef = d5Nr0bdw5(y1d(t))*(y1pd(t))^5 +...
    d4Nr0bdw4(y1d(t))*4*(y1pd(t))^3*y1ppd(t) +...
    6*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppd(t) +...
    6*d3Nr0bdw3(y1d(t))*2*(y1pd(t))*(y1ppd(t))^2 +...
    6*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppd(t) +...
    3*d3Nr0bdw3(y1d(t))*y1pd(t)*(y1ppd(t))^2 +...
    3*d2Nr0bdw2(y1d(t))*2*(y1ppd(t))*y1pppd(t) +...
    4*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1d(t)*y1ppppd(t) +...
    dNr0bdw(y1d(t))*y1pppppd(t);
erro = simplify(y2pppppdRef - y2pppppd);



