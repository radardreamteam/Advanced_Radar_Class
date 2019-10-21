%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Generates a QPSK Signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RFSignal, timeVector] = genrateQPSK_Signal(cyclesPerSymbol, sampsPerCycle, fc)

j = 1i;

sampsPerSymbol  = sampsPerCycle * cyclesPerSymbol;
qpskBBSigFreq   = fc / cyclesPerSymbol;
qpskSym1        = repmat((.7+j*.7),1, sampsPerSymbol);
qpskSym2        = repmat((.7-j*.7),1, sampsPerSymbol);
qpskSym3        = repmat((-.7+j*.7),1, sampsPerSymbol);
qpskSym4        = repmat((-.7-j*.7),1, sampsPerSymbol);
qpskSymStream   = cat(2, qpskSym1, qpskSym2, qpskSym3, qpskSym4);
Ls              = 100*length(qpskSymStream);

% QPSK baseband signal
qpskBBSig       = repmat(qpskSymStream,1,Ls/length(qpskSymStream)); 

fs = sampsPerCycle * fc;

% time vector axis (s)
timeVector = (0:(Ls-1))/fs; 
% Upconvert to RF (carier frequency)
RFSignal = qpskBBSig.*exp(2*pi*j*fc*timeVector); 

return