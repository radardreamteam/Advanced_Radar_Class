%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Generates a QPSK Signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RFSignal, timeVector] = ZgenrateQPSK_Signal(cyclesPerSymbol, sampsPerCycle, fc,sigNumber)

j = 1i;

sampsPerSymbol  = sampsPerCycle * cyclesPerSymbol;
qpskSym = zeros(4,1,sampsPerSymbol);
qpskBBSigFreq   = fc / cyclesPerSymbol;
qpskSym(1,:,:)  = repmat((.7+j*.7),1, sampsPerSymbol);
qpskSym(2,:,:)  = repmat((.7-j*.7),1, sampsPerSymbol);
qpskSym(3,:,:)  = repmat((-.7+j*.7),1, sampsPerSymbol);
qpskSym(4,:,:)  = repmat((-.7-j*.7),1, sampsPerSymbol);

qpskSymStream   = repmat((.7+j*.7),1, sampsPerSymbol);
for k = 1:1:sigNumber-1
    m   = randperm(4,1);
    add = reshape(qpskSym(m,:,:),[1,160]);
    qpskSymStream = cat (2,qpskSymStream,add);
end
    
%qpskSymStream   = cat(2, qpskSym1, qpskSym2, qpskSym3, qpskSym4);
Ls              = length(qpskSymStream);

% QPSK baseband signal
%qpskBBSig       = repmat(qpskSymStream,1,Ls/length(qpskSymStream)); 
qpskBBSig = qpskSymStream ;
fs        = sampsPerCycle * fc;

% time vector axis (s)
timeVector = (0:(Ls-1))/fs; 
% Upconvert to RF (carier frequency)
RFSignal = qpskBBSig.*exp(2*pi*j*fc*timeVector); 

return