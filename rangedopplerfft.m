%Generates the range doppler map for a given receiver over a specified
%range of inputs

function [rdmap, ranges, freqs] = rangedopplerfft(tar, freq, range, freqs, rxdirect)
%Inputs
%   phi is an Nx1 column vector of the samples received at receiver of
%   interest.
%   freq is the sampling freq in Hz
%   range is the maximum range of interest in metres. range>0
%   freqs is similar to ranges except for frequencies of interest. freqs is
%   in Hz

%Implemented using FFT's:
%Rxy(t) = x(t) (conv) y*(-t)
%fft(Rxy) = X(f) Y*(f)
%Rxy(t) = ifft( X(f)Y*(f) )


%CONSTANTS
%C = Speed of light (in m/s)
C = 299792458;

%The period (ie time for 1 sample)
dt = 1/freq;

%Work out the maximum ranges
travel_time = range/C;
n = ceil(travel_time/dt); %The number of shifts we need to do
N = n + size(tar,1); %Size to pad to for fft

%Work out the non direct (ie target path) signal


%Pad out the inputs so we don't get overlap in frequency domain
Fdirect = conj(fft(rxdirect,N)); %Fourier transform of direct path signal

%The range doppler map
rdmap = zeros(n,numel(freqs));

%Time vector for frequency shifting
t = 0:dt:dt*(numel(tar)-1);
t = t';

%We will be taking conjugates in frequency domain, so
freqs = -freqs;

%j=0; %index
parfor f = 1:length(freqs)
    %j=j+1; %index
    
    %frequency correct the signal
    fc_sig = tar.*exp(-1i*2*pi*freqs(f)*t);
%     for k =1:1:length(t)
%         fc_sig(k) = tar(k).*exp(-1i*2*pi*f*t(k));
%     end
    %Compute conjugate of FFT
    Fsig = fft(fc_sig,N);
    
    %Find the correlation for all ranges
    cor = ifft(Fdirect.*Fsig);
    
    rdmap(:,f) = cor(1:n);
end

%Get real values
rdmap = abs(rdmap);
% rdmap = abs(rdmap) + real(rdmap);

%return the ranges
ranges = linspace(0,range,n);
