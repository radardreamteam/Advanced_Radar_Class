%Generates the range doppler map for a given receiver over a specified
%range of inputs

function [rdmap, ranges, freqShift] = rdfft(tempEcho, samplingFreq, ranges, freqShift, dirPath)
%Inputs
%   indirPath is an Nx1 column vector of the samples received at receiver of
%   interest.
%   samplingFreq is the sampling freq in Hz
%   range is the maximum range of interest in metres. range>0
%   freqs is similar to ranges except for frequencies of interest. freqs is
%   in Hz

%Implemented using FFT's:
%Rxy(t) = x(t) (conv) y*(-t)
%fft(Rxy) = X(f) Y*(f)
%Rxy(t) = ifft( X(f)Y*(f) )

C = 299792458;
dt = 1/samplingFreq;
t = 0:dt:dt*(64000-1);
t = t';

p=0;
for f = freqShift
    p=p+1; %index
    i=0;
    
    %frequency correc the signal
    fc_sig = tempEcho.*exp(-1i*2*pi*f*t);
    
    for r = ranges
        i=i+1; %index stuff
        
        %Correct the signal range
        travel_time = r/C;
        %number of samples delayed by:
        num_delay = round(travel_time/dt);
        %range corrected signal
        sig = fc_sig(num_delay+1:end);
        comp = phi(1:numel(sig));
        %Calculate metric
        %       rdmap(i,j) = sig'*phi(1:numel(sig));
%         [~,~,R] = canoncorr(sig,comp);
%         rdmap(i,j) = R; %% corr or corr2
%         temp = corrcoef(sig,comp);
%         rdmap(i,j) = abs(temp(1,2));
        rdmap(i,p) = abs(corr(sig,comp));
        
    end
end