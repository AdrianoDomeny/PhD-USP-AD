function [ Kfn ] = fModKfn( KsnM1, KfnM1, A22nM1, A21nM1, B2nM1  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
KbnM1 = KsnM1 + KfnM1*(A22nM1\(A21nM1 + B2nM1*KsnM1));
Kfn = [KbnM1, KfnM1];
end

