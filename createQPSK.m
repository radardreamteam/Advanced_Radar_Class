%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developed by Rick Pooler <Rick.Pooler@jhuapl.edu> based on a script 
% written by Oscar Summerlock 
%       Generates a qpsk signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RFSignal, t] = createQPSK(cyclesPerSymbol, sampsPerCycle, fc)

j = 1i;

sampsPerSymbol = sampsPerCycle * cyclesPerSymbol;
qpskBBSigFreq = fc / cyclesPerSymbol;
qpskSym1 = repmat((.7+j*.7),1, sampsPerSymbol);
qpskSym2 = repmat((.7-j*.7),1, sampsPerSymbol);
qpskSym3 = repmat((-.7+j*.7),1, sampsPerSymbol);
qpskSym4 = repmat((-.7-j*.7),1, sampsPerSymbol);
qpskSymStream = cat(2, qpskSym1, qpskSym2, qpskSym3, qpskSym4);
Ls = 100*length(qpskSymStream);
qpskBBSig = repmat(qpskSymStream,1,Ls/length(qpskSymStream)); % qpsk baseband signal

fs = sampsPerCycle * fc;
t = (0:(Ls-1))/fs; % time axis (s)
RFSignal = qpskBBSig.*exp(2*pi*j*fc*t); % upconvert to RF

return