function [DSC] = SCA(Surv,Ref,K,d,M)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
T = length(Ref);
X = zeros(T,K) + j*zeros(T,K);
P = eye(T);
Q = zeros(T,K)+ j*zeros(T,K);
for i = 1: K
    temp = zeros(1,d*(i-1));
    X(:,i) = [temp Ref(1:T-d*(i-1))]';
end

for i = M:-1:1
    x = P * X(:,i);
    Q = (eye(T) - ((x*x')/(x'*x)));
    P = P * Q;
end


DSC = P * Surv';

return

