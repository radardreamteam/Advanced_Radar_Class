%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Weiner Filter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yw, ew] = weinerf(order, observed, ref)
bw = firwiener(order,observed,ref); % Optimal FIR Wiener filter
yw = filter(bw,1,observed);   % Estimate of x using Wiener filter
ew = ref - yw;            % Estimate of actual sinusoid
return