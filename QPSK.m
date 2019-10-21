function [BaseBandSignal,SymbolFreq] = QPSK(symbols,fb,samplingF,cyclesPerSymbol)

sampsPerCycle = round(samplingF/fb);
j = 1i;
sampsPerSymbol = sampsPerCycle * cyclesPerSymbol;
qpskBBSig_real = zeros(1,length(symbols)*sampsPerSymbol);
qpskBBSig = complex(qpskBBSig_real,0);
for ind = 1:length(symbols)
    currentSymbol = symbols(ind);
    currentSig = repmat(currentSymbol,1,sampsPerSymbol);
    qpskBBSig(1+(ind-1)*sampsPerSymbol:ind*sampsPerSymbol) = currentSig;

end
qpskBBSigFreq = fb / cyclesPerSymbol;

%%
Ls = length(qpskBBSig);
fs = sampsPerCycle * fb;
t = (0:(Ls-1))/fs; % time axis (s)
BaseBandSignal = qpskBBSig.*exp(2*pi*j*fb*t);
SymbolFreq = fb / cyclesPerSymbol;
return