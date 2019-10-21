function [Signal] = generate_DVBS(integration_Time,fb,SamplingFreq,maxSymbol,numOfFrames,pilotOn)
data_num = integration_Time*SamplingFreq;
cyclesPerSymbol = 27;
% numOfFrames = 3;
% pilotOn     = 1;
symbols = DVBstream(numOfFrames,pilotOn,maxSymbol);
[BaseBandSignal,~] = QPSK(symbols,fb,SamplingFreq,cyclesPerSymbol);
Signal = BaseBandSignal(1:round(data_num));

return