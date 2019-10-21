%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Weiner Filter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yw, ew] = weinerf(order, reference, observed)
bw = firwiener(order,reference,observed); % Optimal FIR Wiener filter
yw = filter(bw,1,reference);   % Estimate of x using Wiener filter
ew = observed - yw;            % Estimate of actual sinusoid
return