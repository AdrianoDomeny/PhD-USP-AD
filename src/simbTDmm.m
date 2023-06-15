clear
clc

syms alfa gama_TD beta_TD t t09 wRef
nu = t/t09;
fi = tanh(alfa*nu) - (alfa*nu)/((alfa^2*nu^2)/3 + 1);
y1d = wRef*fi;
w = y1d;
Nr = gama_TD*beta_TD*exp(w/beta_TD) - gama_TD*(w + beta_TD);
y2d = Nr;
%%
%(1) y2pd = @(t)(dNr0bdw(y1d(t))*y1pd(t))
dNr0bdw = gama_TD*(exp(w/beta_TD) - 1);
y1pd = diff(y1d,t);
y2pd = dNr0bdw*y1pd;
%{
y2pdRef = diff(y2d,t);
erro = simplify(y2pdRef - y2pd);
%}
%%
%(2) y2ppd = @(t)(d2Nr0bdw2(y1d(t))*(y1pd(t))^2 + dNr0bdw(y1d(t))*y1ppd(t))
d2Nr0bdw2 = (gama_TD*exp(w/beta_TD))/beta_TD;
y1ppd = diff(y1pd,t);
y2ppd = d2Nr0bdw2*(y1pd)^2 + dNr0bdw*y1ppd;
%{
y2ppdRef = diff(y2pd,t);
erro = simplify(y2ppdRef - y2ppd);
%}
%%
%(3) 
%{
y2pppd = @(t)(d3Nr0bdw3(y1d(t))*(y1pd(t))^3 +...
    3*d2Nr0bdw2(y1d(t))*(y1pd(t))*y1ppd(t) +...
    dNr0bdw(y1d(t))*y1pppd(t));
%}
d3Nr0bdw3 = (gama_TD*exp(w/beta_TD))/beta_TD^2;
y1pppd = diff(y1ppd,t);
y2pppd = d3Nr0bdw3*(y1pd)^3 +...
    3*d2Nr0bdw2*(y1pd)*y1ppd +...
    dNr0bdw*y1pppd;
%{
y2pppdRef = diff(y2ppd,t);
erro = simplify(y2pppdRef - y2pppd);
%}
%%
%(4) 
%{
y2ppppd = @(t)(d4Nr0bdw4(y1d(t))*(y1pd(t))^4 +...
    6*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1ppd(t) +...
    3*d2Nr0bdw2(y1d(t))*(y1ppd(t))^2 +...
    4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1pppd(t) +...
    dNr0bdw(y1d(t))*y1ppppd(t));
%}
d4Nr0bdw4 = (gama_TD*exp(w/beta_TD))/beta_TD^3;
y1ppppd = diff(y1pppd,t);
y2ppppd = d4Nr0bdw4*(y1pd)^4 +...
    6*d3Nr0bdw3*(y1pd)^2*y1ppd +...
    3*d2Nr0bdw2*(y1ppd)^2 +...
    4*d2Nr0bdw2*y1pd*y1pppd +...
    dNr0bdw*y1ppppd;
%{
y2ppppdRef = diff(y2pppd,t);
erro = simplify(y2ppppdRef - y2ppppd);
%}
%%
%(5) 
%{
y2pppppd = @(t)(d5Nr0bdw5(y1d(t))*(y1pd(t))^5 +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppd(t) +...
    15*d3Nr0bdw3(y1d(t))*y1pd(t)*(y1ppd(t))^2 +...
    10*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppd(t) +...
    10*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppd(t) +...
    dNr0bdw(y1d(t))*y1pppppd(t));
%}
d5Nr0bdw5 = (gama_TD*exp(w/beta_TD))/beta_TD^4;
y1pppppd = diff(y1ppppd,t);
y2pppppd = d5Nr0bdw5*(y1pd)^5 +...
    10*d4Nr0bdw4*(y1pd)^3*y1ppd +...
    15*d3Nr0bdw3*y1pd*(y1ppd)^2 +...
    10*d3Nr0bdw3*(y1pd)^2*y1pppd +...
    10*d2Nr0bdw2*y1ppd*y1pppd +...
    4*d2Nr0bdw2*y1pd*y1ppppd +...
    d2Nr0bdw2*y1pd*y1ppppd +...
    dNr0bdw*y1pppppd;
%{
y2pppppdRef = diff(y2ppppd,t);
erro = simplify(y2pppppdRef - y2pppppd);
%}
%%
%(6) 
%{
y2ppppppd = @(t)(d6Nr0bdw6(y1d(t))*(y1pd(t))^6 +...
    d5Nr0bdw5(y1d(t))*5*(y1pd(t))^4*y1ppd(t) +...
    10*d5Nr0bdw5(y1d(t))*(y1pd(t))^4*y1ppd(t) +...
    10*d4Nr0bdw4(y1d(t))*3*(y1pd(t))^2*(y1ppd(t))^2 +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1pppd(t) +...
    15*d4Nr0bdw4(y1d(t))*(y1pd(t))^2*(y1ppd(t))^2 +...
    15*d3Nr0bdw3(y1d(t))*(y1ppd(t))^3 +...
    15*d3Nr0bdw3(y1d(t))*y1pd(t)*2*y1ppd(t)*y1pppd(t) +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1pppd(t) +...
    10*d3Nr0bdw3(y1d(t))*2*y1pd(t)*y1ppd(t)*y1pppd(t) +...
    10*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1ppppd(t) +...
    10*d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1pppd(t) +...
    10*d2Nr0bdw2(y1d(t))*(y1pppd(t))^2 +...
    10*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1ppppd(t) +...
    4*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1ppppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1ppppd(t) +...
    4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1pppppd(t) +...
    d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1ppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1ppd(t)*y1ppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1pd(t)*y1pppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1pd(t)*y1pppppd(t) +...
    dNr0bdw(y1d(t))*y1ppppppd(t));
%}
d6Nr0bdw6 = (gama_TD*exp(w/beta_TD))/beta_TD^5;
y1ppppppd = diff(y1pppppd,t);
y2ppppppd = d6Nr0bdw6*(y1pd)^6 +...
    d5Nr0bdw5*5*(y1pd)^4*y1ppd +...
    10*d5Nr0bdw5*(y1pd)^4*y1ppd +...
    10*d4Nr0bdw4*3*(y1pd)^2*(y1ppd)^2 +...
    10*d4Nr0bdw4*(y1pd)^3*y1pppd +...
    15*d4Nr0bdw4*(y1pd)^2*(y1ppd)^2 +...
    15*d3Nr0bdw3*(y1ppd)^3 +...
    15*d3Nr0bdw3*y1pd*2*y1ppd*y1pppd +...
    10*d4Nr0bdw4*(y1pd)^3*y1pppd +...
    10*d3Nr0bdw3*2*y1pd*y1ppd*y1pppd +...
    10*d3Nr0bdw3*(y1pd)^2*y1ppppd +...
    10*d3Nr0bdw3*y1pd*y1ppd*y1pppd +...
    10*d2Nr0bdw2*(y1pppd)^2 +...
    10*d2Nr0bdw2*y1ppd*y1ppppd +...
    4*d3Nr0bdw3*(y1pd)^2*y1ppppd +...
    4*d2Nr0bdw2*y1ppd*y1ppppd +...
    4*d2Nr0bdw2*y1pd*y1pppppd +...
    d3Nr0bdw3*(y1pd)^2*y1ppppd +...
    d2Nr0bdw2*y1ppd*y1ppppd +...
    d2Nr0bdw2*y1pd*y1pppppd +...
    d2Nr0bdw2*y1pd*y1pppppd +...
    dNr0bdw*y1ppppppd;
%{
y2ppppppdRef = diff(y2pppppd,t);
erro = simplify(y2ppppppdRef - y2ppppppd);
%}
%%
%(7) 
d7Nr0bdw7 = (gama_TD*exp(w/beta_TD))/beta_TD^6;
y1pppppppd = diff(y1ppppppd,t);
%{
y2pppppppd = @(t)(d7Nr0bdw7(y1d(t))*(y1pd(t))^7 + d6Nr0bdw6(y1d(t))*6*(y1pd(t))^5*y1ppd(t) +...
    d6Nr0bdw6(y1d(t))*5*(y1pd(t))^5*y1ppd(t) + d5Nr0bdw5(y1d(t))*5*4*(y1pd(t))^3*(y1ppd(t))^2 + d5Nr0bdw5(y1d(t))*5*(y1pd(t))^4*y1pppd(t) +...
    10*d6Nr0bdw6(y1d(t))*(y1pd(t))^5*y1ppd(t) + 10*d5Nr0bdw5(y1d(t))*4*(y1pd(t))^3*(y1ppd(t))^2 + 10*d5Nr0bdw5(y1d(t))*(y1pd(t))^4*y1pppd(t) +...
    10*d5Nr0bdw5(y1d(t))*3*(y1pd(t))^3*(y1ppd(t))^2 + 10*d4Nr0bdw4(y1d(t))*3*2*(y1pd(t))*(y1ppd(t))^3 + 10*d4Nr0bdw4(y1d(t))*3*(y1pd(t))^2*2*(y1ppd(t))*y1pppd(t) +...
    10*d5Nr0bdw5(y1d(t))*(y1pd(t))^4*y1pppd(t) + 10*d4Nr0bdw4(y1d(t))*3*(y1pd(t))^2*y1ppd(t)*y1pppd(t) + 10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t) +...
    15*d5Nr0bdw5(y1d(t))*(y1pd(t))^3*(y1ppd(t))^2 + 15*d4Nr0bdw4(y1d(t))*2*(y1pd(t))*(y1ppd(t))^3 + 15*d4Nr0bdw4(y1d(t))*(y1pd(t))^2*2*(y1ppd(t))*y1pppd(t) +...
    15*d4Nr0bdw4(y1d(t))*y1d(t)*(y1ppd(t))^3 + 15*d3Nr0bdw3(y1d(t))*3*(y1ppd(t))^2*y1pppd(t) +...
    15*d4Nr0bdw4(y1d(t))*(y1pd(t))^2*2*y1ppd(t)*y1pppd(t) + 15*d3Nr0bdw3(y1d(t))*2*(y1ppd(t))^2*y1pppd(t) + 15*d3Nr0bdw3(y1d(t))*y1pd(t)*2*(y1pppd(t))^2 + 15*d3Nr0bdw3(y1d(t))*y1pd(t)*2*y1ppd(t)*y1ppppd(t) +...
    10*d5Nr0bdw5(y1d(t))*(y1pd(t))^4*y1pppd(t) + 10*d4Nr0bdw4(y1d(t))*3*(y1pd(t))^2*y1ppd(t)*y1pppd(t) + 10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t) +...
    10*d4Nr0bdw4(y1d(t))*2*(y1pd(t))^2*y1ppd(t)*y1pppd(t) + 10*d3Nr0bdw3(y1d(t))*2*(y1ppd(t))^2*y1pppd(t) + 10*d3Nr0bdw3(y1d(t))*2*y1pd(t)*(y1pppd(t))^2 + 10*d3Nr0bdw3(y1d(t))*2*y1pd(t)*y1ppd(t)*y1ppppd(t) +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t) + 10*d3Nr0bdw3(y1d(t))*2*(y1pd(t))*y1ppd(t)*y1ppppd(t) + 10*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t) +...
    10*d4Nr0bdw4(y1d(t))*(y1pd(t))^2*y1ppd(t)*y1pppd(t) + 10*d3Nr0bdw3(y1d(t))*(y1ppd(t))^2*y1pppd(t) + 10*d3Nr0bdw3(y1d(t))*y1pd(t)*(y1pppd(t))^2 + 10*d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1ppppd(t) +...
    10*d3Nr0bdw3(y1d(t))*y1pd(t)*(y1pppd(t))^2 + 10*d2Nr0bdw2(y1d(t))*2*(y1pppd(t))*y1ppppd(t) +...
    10*d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1ppppd(t) + 10*d2Nr0bdw2(y1d(t))*y1pppd(t)*y1ppppd(t) + 10*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) +...
    4*d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t) + 4*d3Nr0bdw3(y1d(t))*2*(y1pd(t))*y1ppd(t)*y1ppppd(t) + 4*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t) +...
    4*d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1ppppd(t) + 4*d2Nr0bdw2(y1d(t))*y1pppd(t)*y1ppppd(t) + 4*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) +...
    4*d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t) + 4*d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) + 4*d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppppd(t) +...
    d4Nr0bdw4(y1d(t))*(y1pd(t))^3*y1ppppd(t) + d3Nr0bdw3(y1d(t))*2*(y1pd(t))*y1ppd(t)*y1ppppd(t) + d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t) +...
    d3Nr0bdw3(y1d(t))*y1pd(t)*y1ppd(t)*y1ppppd(t) + d2Nr0bdw2(y1d(t))*y1pppd(t)*y1ppppd(t) + d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) +...
    d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t) + d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) + d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppppd(t) +...
    d3Nr0bdw3(y1d(t))*(y1pd(t))^2*y1pppppd(t) + d2Nr0bdw2(y1d(t))*y1ppd(t)*y1pppppd(t) + d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppppd(t) +...
    d2Nr0bdw2(y1d(t))*y1pd(t)*y1ppppppd(t) + dNr0bdw(y1d(t))*y1pppppppd(t));
%}
y2ppppppd01 = d6Nr0bdw6*(y1pd)^6;
y2ppppppd02 = d5Nr0bdw5*5*(y1pd)^4*y1ppd;
y2ppppppd03 = 10*d5Nr0bdw5*(y1pd)^4*y1ppd;
y2ppppppd04 = 10*d4Nr0bdw4*3*(y1pd)^2*(y1ppd)^2;
y2ppppppd05 = 10*d4Nr0bdw4*(y1pd)^3*y1pppd;
y2ppppppd06 = 15*d4Nr0bdw4*(y1pd)^2*(y1ppd)^2;
y2ppppppd07 = 15*d3Nr0bdw3*(y1ppd)^3;
y2ppppppd08 = 15*d3Nr0bdw3*y1pd*2*y1ppd*y1pppd;
y2ppppppd09 = 10*d4Nr0bdw4*(y1pd)^3*y1pppd;
y2ppppppd10 = 10*d3Nr0bdw3*2*y1pd*y1ppd*y1pppd;
y2ppppppd11 = 10*d3Nr0bdw3*(y1pd)^2*y1ppppd;
y2ppppppd12 = 10*d3Nr0bdw3*y1pd*y1ppd*y1pppd;
y2ppppppd13 = 10*d2Nr0bdw2*(y1pppd)^2;
y2ppppppd14 = 10*d2Nr0bdw2*y1ppd*y1ppppd;
y2ppppppd15 = 4*d3Nr0bdw3*(y1pd)^2*y1ppppd;
y2ppppppd16 = 4*d2Nr0bdw2*y1ppd*y1ppppd;
y2ppppppd17 = 4*d2Nr0bdw2*y1pd*y1pppppd;
y2ppppppd18 = d3Nr0bdw3*(y1pd)^2*y1ppppd;
y2ppppppd19 = d2Nr0bdw2*y1ppd*y1ppppd;
y2ppppppd20 = d2Nr0bdw2*y1pd*y1pppppd;
y2ppppppd21 = d2Nr0bdw2*y1pd*y1pppppd;
y2ppppppd22 = dNr0bdw*y1ppppppd;
%{
erro = simplify(y2ppppppd -...
    (y2ppppppd01 + y2ppppppd02 + y2ppppppd03 + y2ppppppd04 + y2ppppppd05 +...
    y2ppppppd06 + y2ppppppd07 + y2ppppppd08 + y2ppppppd09 + y2ppppppd10 +...
    y2ppppppd11 + y2ppppppd12 + y2ppppppd13 + y2ppppppd14 + y2ppppppd15 +...
    y2ppppppd16 + y2ppppppd17 + y2ppppppd18 + y2ppppppd19 + y2ppppppd20 +...
    y2ppppppd21 + y2ppppppd22));
%}
y2pppppppd01 = d7Nr0bdw7*(y1pd)^7 + d6Nr0bdw6*6*(y1pd)^5*y1ppd;
%{
y2pppppppd01Ref = diff(y2ppppppd01,t);
erro = simplify(y2pppppppd01Ref - y2pppppppd01);
%}
y2pppppppd02 = d6Nr0bdw6*5*(y1pd)^5*y1ppd + d5Nr0bdw5*5*4*(y1pd)^3*y1ppd^2 + d5Nr0bdw5*5*(y1pd)^4*y1pppd;
%{
y2pppppppd02Ref = diff(y2ppppppd02,t);
erro = simplify(y2pppppppd02Ref - y2pppppppd02);
%}
y2pppppppd03 = 10*d6Nr0bdw6*(y1pd)^5*y1ppd + 10*d5Nr0bdw5*4*(y1pd)^3*y1ppd^2 + 10*d5Nr0bdw5*(y1pd)^4*y1pppd;
%{
y2pppppppd03Ref = diff(y2ppppppd03,t);
erro = simplify(y2pppppppd03Ref - y2pppppppd03);
%}
y2pppppppd04 = 10*d5Nr0bdw5*3*(y1pd)^3*(y1ppd)^2 + 10*d4Nr0bdw4*3*2*(y1pd)*(y1ppd)^3 + 10*d4Nr0bdw4*3*(y1pd)^2*2*(y1ppd)*y1pppd;
%{
y2pppppppd04Ref = diff(y2ppppppd04,t);
erro = simplify(y2pppppppd04Ref - y2pppppppd04);
%}
y2pppppppd05 = 10*d5Nr0bdw5*(y1pd)^4*y1pppd + 10*d4Nr0bdw4*3*(y1pd)^2*y1ppd*y1pppd + 10*d4Nr0bdw4*(y1pd)^3*y1ppppd;
%{
y2pppppppd05Ref = diff(y2ppppppd05,t);
erro = simplify(y2pppppppd05Ref - y2pppppppd05);
%}
y2pppppppd06 = 15*d5Nr0bdw5*(y1pd)^3*(y1ppd)^2 + 15*d4Nr0bdw4*2*(y1pd)*(y1ppd)^3 + 15*d4Nr0bdw4*(y1pd)^2*2*(y1ppd)*y1pppd;
%{
y2pppppppd06Ref = diff(y2ppppppd06,t);
erro = simplify(y2pppppppd06Ref - y2pppppppd06);
%}
y2pppppppd07 = 15*d4Nr0bdw4*y1pd*(y1ppd)^3 + 15*d3Nr0bdw3*3*(y1ppd)^2*y1pppd;
%{
y2pppppppd07Ref = diff(y2ppppppd07,t);
erro = simplify(y2pppppppd07Ref - y2pppppppd07);
%}
y2pppppppd08 = 15*d4Nr0bdw4*y1pd^2*2*y1ppd*y1pppd + 15*d3Nr0bdw3*2*y1ppd^2*y1pppd + 15*d3Nr0bdw3*y1pd*2*y1pppd^2 + 15*d3Nr0bdw3*y1pd*2*y1ppd*y1ppppd;
%{
y2pppppppd08Ref = diff(y2ppppppd08,t);
erro = simplify(y2pppppppd08Ref - y2pppppppd08);
%}
y2pppppppd09 = 10*d5Nr0bdw5*(y1pd)^4*y1pppd + 10*d4Nr0bdw4*3*(y1pd)^2*y1ppd*y1pppd + 10*d4Nr0bdw4*(y1pd)^3*y1ppppd;
%{
y2pppppppd09Ref = diff(y2ppppppd09,t);
erro = simplify(y2pppppppd09Ref - y2pppppppd09);
%}
y2pppppppd10 = 10*d4Nr0bdw4*2*y1pd^2*y1ppd*y1pppd + 10*d3Nr0bdw3*2*y1ppd^2*y1pppd + 10*d3Nr0bdw3*2*y1pd*y1pppd^2 + 10*d3Nr0bdw3*2*y1pd*y1ppd*y1ppppd;
%{
y2pppppppd10Ref = diff(y2ppppppd10,t);
erro = simplify(y2pppppppd10Ref - y2pppppppd10);
%}
y2pppppppd11 = 10*d4Nr0bdw4*(y1pd)^3*y1ppppd + 10*d3Nr0bdw3*2*(y1pd)*y1ppd*y1ppppd + 10*d3Nr0bdw3*(y1pd)^2*y1pppppd;
%{
y2pppppppd11Ref = diff(y2ppppppd11,t);
erro = simplify(y2pppppppd11Ref - y2pppppppd11);
%}
y2pppppppd12 = 10*d4Nr0bdw4*y1pd^2*y1ppd*y1pppd + 10*d3Nr0bdw3*y1ppd^2*y1pppd + 10*d3Nr0bdw3*y1pd*y1pppd^2 + 10*d3Nr0bdw3*y1pd*y1ppd*y1ppppd;
%{
y2pppppppd12Ref = diff(y2ppppppd12,t);
erro = simplify(y2pppppppd12Ref - y2pppppppd12);
%}
y2pppppppd13 = 10*d3Nr0bdw3*y1pd*(y1pppd)^2 + 10*d2Nr0bdw2*2*(y1pppd)*y1ppppd;
%{
y2pppppppd13Ref = diff(y2ppppppd13,t);
erro = simplify(y2pppppppd13Ref - y2pppppppd13);
%}
y2pppppppd14 = 10*d3Nr0bdw3*y1pd*y1ppd*y1ppppd + 10*d2Nr0bdw2*y1pppd*y1ppppd + 10*d2Nr0bdw2*y1ppd*y1pppppd;
%{
y2pppppppd14Ref = diff(y2ppppppd14,t);
erro = simplify(y2pppppppd14Ref - y2pppppppd14);
%}
y2pppppppd15 = 4*d4Nr0bdw4*(y1pd)^3*y1ppppd + 4*d3Nr0bdw3*2*(y1pd)*y1ppd*y1ppppd + 4*d3Nr0bdw3*(y1pd)^2*y1pppppd;
%{
y2pppppppd15Ref = diff(y2ppppppd15,t);
erro = simplify(y2pppppppd15Ref - y2pppppppd15);
%}
y2pppppppd16 = 4*d3Nr0bdw3*y1pd*y1ppd*y1ppppd + 4*d2Nr0bdw2*y1pppd*y1ppppd + 4*d2Nr0bdw2*y1ppd*y1pppppd;
%{
y2pppppppd16Ref = diff(y2ppppppd16,t);
erro = simplify(y2pppppppd16Ref - y2pppppppd16);
%}
y2pppppppd17 = 4*d3Nr0bdw3*y1pd^2*y1pppppd + 4*d2Nr0bdw2*y1ppd*y1pppppd + 4*d2Nr0bdw2*y1pd*y1ppppppd;
%{
y2pppppppd17Ref = diff(y2ppppppd17,t);
erro = simplify(y2pppppppd17Ref - y2pppppppd17);
%}
y2pppppppd18 = d4Nr0bdw4*(y1pd)^3*y1ppppd + d3Nr0bdw3*2*(y1pd)*y1ppd*y1ppppd + d3Nr0bdw3*(y1pd)^2*y1pppppd;
%{
y2pppppppd18Ref = diff(y2ppppppd18,t);
erro = simplify(y2pppppppd18Ref - y2pppppppd18);
%}
y2pppppppd19 = d3Nr0bdw3*y1pd*y1ppd*y1ppppd + d2Nr0bdw2*y1pppd*y1ppppd + d2Nr0bdw2*y1ppd*y1pppppd;
%{
y2pppppppd19Ref = diff(y2ppppppd19,t);
erro = simplify(y2pppppppd19Ref - y2pppppppd19);
%}
y2pppppppd20 = d3Nr0bdw3*y1pd^2*y1pppppd + d2Nr0bdw2*y1ppd*y1pppppd + d2Nr0bdw2*y1pd*y1ppppppd;
%{
y2pppppppd20Ref = diff(y2ppppppd20,t);
erro = simplify(y2pppppppd20Ref - y2pppppppd20);
%}
y2pppppppd21 = d3Nr0bdw3*y1pd^2*y1pppppd + d2Nr0bdw2*y1ppd*y1pppppd + d2Nr0bdw2*y1pd*y1ppppppd;
%{
y2pppppppd21Ref = diff(y2ppppppd21,t);
erro = simplify(y2pppppppd21Ref - y2pppppppd21);
%}
y2pppppppd22 = d2Nr0bdw2*y1pd*y1ppppppd + dNr0bdw*y1pppppppd;
%{
y2pppppppd22Ref = diff(y2ppppppd22,t);
erro = simplify(y2pppppppd22Ref - y2pppppppd22);
%}
y2pppppppd = y2pppppppd01 + y2pppppppd02 + y2pppppppd03 + y2pppppppd04 +...
    y2pppppppd05 + y2pppppppd06 + y2pppppppd07 + y2pppppppd08 +...
    y2pppppppd09 + y2pppppppd10 + y2pppppppd11 + y2pppppppd12 +...
    y2pppppppd13 + y2pppppppd14 + y2pppppppd15 + y2pppppppd16 +...
    y2pppppppd17 + y2pppppppd18 + y2pppppppd19 + y2pppppppd20 +...
    y2pppppppd21 + y2pppppppd22;
%{
y2pppppppdRef = diff(y2ppppppd,t);
erro = simplify(y2pppppppdRef - y2pppppppd);
%}
%%
syms y2f0 wUf
y2df = y2f0*sin(wUf*t); 
y2pdf = y2f0*wUf*cos(wUf*t);
%{
y2pdfRef = diff(y2df,t);
erro = simplify(y2pdfRef - y2pdf);
%}
y2ppdf = (-1)*y2f0*wUf^2*sin(wUf*t);
%{
y2ppdfRef = diff(y2pdf,t);
erro = simplify(y2ppdfRef - y2ppdf);
%}
y2pppdf = (-1)*y2f0*wUf^3*cos(wUf*t);
%{
y2pppdfRef = diff(y2ppdf,t);
erro = simplify(y2pppdfRef - y2pppdf);
%}
y2ppppdf = (-1)^2*y2f0*wUf^4*sin(wUf*t);
%{
y2ppppdfRef = diff(y2pppdf,t);
erro = simplify(y2ppppdfRef - y2ppppdf);
%}
y2pppppdf = (-1)^2*y2f0*wUf^5*cos(wUf*t);
%{
y2pppppdfRef = diff(y2ppppdf,t);
erro = simplify(y2pppppdfRef - y2pppppdf);
%}
y2ppppppdf = (-1)^3*y2f0*wUf^6*sin(wUf*t);
%{
y2ppppppdfRef = diff(y2pppppdf,t);
erro = simplify(y2ppppppdfRef - y2ppppppdf);
%}
y2pppppppdf = (-1)^3*y2f0*wUf^7*cos(wUf*t);
%{
y2pppppppdfRef = diff(y2ppppppdf,t);
erro = simplify(y2pppppppdfRef - y2pppppppdf);
%}
ad = 1;















