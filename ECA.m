function [DSC] = ECA(Surv,Ref,K,d)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
T = length(Ref);
X = zeros(T,K);
for i = 1: K
    temp = zeros(1,d*(i-1));
    X(:,i) = [temp Ref(1:T-d*(i-1))].';
end
Surv = Surv.';
DSC = Surv - X*((X'*X)\X'*Surv);
return

