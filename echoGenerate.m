function [Echo,offset] = echoGenerate(direct,tx,rx,targetLocation,targetVelocity,fSample,fCarrier)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
            C = 299792458;
            dt = 1/fSample;

            a = targetVelocity;
            b = targetLocation - tx;
            c = targetLocation - rx;
            d = tx-rx;
            rtx = dot(a,b) / norm(b);   %Component of Tx/Targ motion
            rrx = dot(a,c) / norm(c);   %Component of Rx/Targ motion
            FShift = rtx + rrx;
            FShift = FShift*fCarrier/C;
            delay = norm(b)+norm(c)-norm(d);
            offset = round((delay/C)/dt);
            timevec = 0:dt:dt*(length(direct)-1);
            
            Echo = direct.*exp(1i*2*pi*FShift*timevec);
end

