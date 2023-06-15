function [ Faj_Hat_c ] = fEFFadcjc( nElem, le, uj, upj, argumentos, repositorio )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
indTatr = argumentos.indTatr;
Nu = repositorio.Nu;
NtetaX = repositorio.NtetaX;
ad_N = repositorio.ad_N;

uNM1 = ad_N(nElem)*Nu(1,le)*uj;
tetaXpNM1 = ad_N(nElem)*NtetaX(1,le)*upj;

FnNM1 = repositorio.FnbNM1;
TxNM1 = repositorio.TxNM1nElemDc;

Faj_Hat_c = ad_N(nElem)*...
    (FnNM1(uNM1)*Nu(1,le).' +...
    indTatr*TxNM1(nElem,tetaXpNM1,uNM1)*NtetaX(1,le).');
end