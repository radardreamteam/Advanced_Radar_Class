%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generate the DVB-s signal (QPSK modulation) 
%
%   All Eq. references can be found in the following paper:
%       'Use of DVB-S satellite TV signal as a source of opportunity for
%       passive coherent location.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all 

generatePlot = 1;
%% Define Constants
j = 1i;
propSpeed 	= 299792458; % m/s

%% Simulation parameters
% Initial System parameters from Balke Sheet
directSigPower          = -122.94+30;
echoSigPower            = -127; 
refGain                 = 20;
% The goal is to eliminate this attenaution via signal processing
dirPathAttenuation      = -50;%dB, Null
survAntGain             = 32; %40;

%% System Parameters
SystemParameters;
% Update System parameters
echoSigPower            = tgtPower_dBm(1);
directSigPower          = directPathPower_dBm;

%% Choose the Direct Signal Suppression Technique
% 1 - NLMS, 2 - Wiener, 3 - RLS, 4 - RLC, 4 - FBLMS, 5 - ECA, 6 - SCA
whichDSIsuppression = 'FBLMS';

%% Create input signal
sigNumber   = 400;

numOfFrames = 3;
pilotOn     = 1; % insert pilot module
maxSymbol   = 5000; % Max symbol number we generate
integrationTime = 0.0025; % second
BaseBandSignal  = generate_DVBS(integrationTime,3e6,samplingFreq,maxSymbol,numOfFrames,pilotOn);
sigLength   = length(BaseBandSignal);
dirPath     = zeros(1, sigLength + samp_offset);
indirPath   = zeros(1, sigLength + samp_offset);
tempEcho    = zeros(1, sigLength + samp_offset);
tempDir     =  zeros(1, sigLength + samp_offset );
sig         = zeros(1, sampsPerCycle*cyclesPerSymbol*400+samp_offset );
taxis       = (0:(length(dirPath)-1))/samplingFreq;

%% Time domain signal
dirPath(1:sigLength) = BaseBandSignal;
indirPath(samp_offset+1:end)       = dirPath(1:sigLength)*exp(j*phaseOffset).*exp(j*2*pi*FShift*taxis(1:sigLength));
tempEcho(samp_offset+1:end)        = indirPath(samp_offset+1:end); 

%The direct signal from the Tx should arrive at both Rx at the same time,
%there is no shift.
tempDir                            = dirPath;
%make the reference signal
refSignal                          = sqrt(10^((directSigPower+refGain)/10))*tempDir;
%Create surveillance channel
echoSigPower                       = echoSigPower + survAntGain;
echoSignal                         = sqrt(10^(echoSigPower/10))*tempEcho;
directSignal                       = 10^(dirPathAttenuation/10)*10^(directSigPower/10)*tempDir;
survChannel                        = echoSignal + directSignal ;

%% Add noise to the Surveillance and Reference Channels
% Create Noise 
noiseSigSurv     = sqrt(10^(noisePower_dBm/10))*(randn(1,length(tempDir))+j*randn(1,length(tempDir)));
survNoisyChannel = survChannel + noiseSigSurv;

noiseSigRef     = sqrt(10^(noisePower_dBm/10))*(randn(1,length(tempDir))+j*randn(1,length(tempDir)));
NoisyrefSignal  = refSignal + noiseSigRef;
% Generate Doppler Vector to simulate
freqVector      = -500:2:500;

% Frequency domain signal
dirPath_fft = fftshift(abs(fft(dirPath))/length(dirPath));
indirPath_fft = fftshift(abs(fft(indirPath))/length(indirPath));

%%% Csingle Channel Plots
% axes
taxis_sig = (0:(length(dirPath)-1))/samplingFreq;
faxis_sig = linspace(-samplingFreq/2,samplingFreq/2,length(dirPath));

% Time domain of signals
f1 = figure('visible','off'); 
hold on; plot(taxis_sig*1e6, real(dirPath));
hold on; plot(taxis_sig*1e6, real(indirPath)); axis tight; grid
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
title('Time Domain DVB-S Single Channel Signal')
legend('Direct Path','Indirect Path')

% Freq domain of signals
f2 = figure('visible','off'); 
indFreqSide = find((faxis_sig >= 0),1);

hold on; plot(faxis_sig(indFreqSide:end)*1e-6, 20*log10(dirPath_fft(indFreqSide:end)));

axis([0 6 -120 5]); grid
xlabel('Frequency (MHz)')
ylabel('Amplitude (dBm)')
title('DVB-S Single signal Spectrum')


%% Multi- Channel Implementation
% DOA estimation using MUSIC algorithm
% Build a 25*25 rectangular array, and use MUSIC algorithm to seperate these two signals 
array = phased.URA('Size',[25 25],'ElementSpacing',[lamda/2 lamda/2]);
array.Element.FrequencyRange = [500.0e6 5000.0e6];
% Target's angle
doa1 = [10;10];
% Given the direct signal in the surveillance channel is not helping, we
% want to place it as far away possible in the sidelobes.
doa2 = [80; -80]; %[30;-20];
survChannelArray  = collectPlaneWave(array,[echoSignal',directSignal'],[doa1,doa2],txSat(1).freq); %samplingFreq);
noiseArray        = sqrt(10^((noisePower_dBm)/10))*(randn(length(tempDir),625)+j*randn(length(tempDir),625));
survChannelArray  = survChannelArray + noiseArray ;

%% DIRECT SIGNAL SUPPRESSION
% NLMS Implementation on Multi-Channel
% We must assume we know the doa of the direct path signal
direct_distribute = collectPlaneWave(array,directSignal',doa2,samplingFreq);
noiseArray        = sqrt(10^(noisePower_dBm/10))*(randn(length(tempDir),625)+j*randn(length(tempDir),625));
direct_distribute = direct_distribute + noiseArray;
% NewsurvChannelArray = zeros(size(survChannelArray,1),size(survChannelArray,2));
tic
switch whichDSIsuppression
    case 'NLMS'
        parfor i = 1:size(survChannelArray,2)
            filterOrder = 32;
            nlms = dsp.LMSFilter('Length',filterOrder,'Method','Normalized LMS','StepSizeSource','Input port');
            nlms.reset();
            [outputNLMS,errNLMS,weightsNLMS] = nlms(direct_distribute(:,i),survChannelArray(:,i),0.001);
            NewsurvChannelArray(:,i) = errNLMS;
        end
    case 'Wiener'
        parfor i = 1:size(survChannelArray,2)
            filterOrder = 32;
            [outputW ,errW] = wienerf(filterOrder,direct_distribute(:,i),survChannelArray(:,i));
            NewsurvChannelArray(:,i) = errW;
        end
    case 'RLS'
        
    case 'RLC'
        
    case 'FBLMS'
        filterLength = 32;
        [numRow,numCol] = size(survChannelArray);
        remData = mod(numRow,filterLength);
        truncatedData = survChannelArray(1:(numRow - remData),:);
        directDistData = direct_distribute(1:(numRow - remData),:);
        parfor i = 1:size(truncatedData,2)            
            fdaf = dsp.FrequencyDomainAdaptiveFilter('Length',filterLength,'BlockLength',filterLength,'StepSize',0.001);
            [outputFBLMS,errFBLMS] = fdaf(directDistData(:,i),truncatedData(:,i));
            NewsurvChannelArray(:,i) = errFBLMS;
        end
    case 'ECA'
        K = 100;
        d =1;
        parfor i = 1:size(survChannelArray,2)
        
        NewsurvChannelArray(:,i) = ECA(survChannelArray(:,i)',direct_distribute(:,i)',K,d);
        end
    case 'SCA'
        
    otherwise
        disp('Not avalid DSI Algorithm')
end
toc
%% Range-Doppler MAP 
%compensate the phase shift to each element
reference_distribute = collectPlaneWave(array,NoisyrefSignal',doa2,samplingFreq);
% multi_channel_rdmap = zeros(size(rdmap,1),size(rdmap,2),size(reference_distribute,2));
for i = size(reference_distribute,2)
    [rdmap, ranges, freqs] = rangedopplerfft(NewsurvChannelArray(:,i),samplingFreq , 2*timeDelay*propSpeed , freqVector, reference_distribute(:,i));
    multi_channel_rdmap(:,:,i) = rdmap;
end
[X,Y] = meshgrid(ranges, freqs);
rdmap_compansated = sum(multi_channel_rdmap,3);

f3 = figure('Name',['Phased array DVB-S After ' whichDSIsuppression ' DSI'],'visible','off'); 
contourf(X*1e-3,Y,rdmap_compansated')
xlabel('Range (Km)')
ylabel('Doppler shift (Hz)')
title(['RDM After ' whichDSIsuppression ' DSI Suppression'])

f4 = figure('Name',['Phased array DVB-S After ' whichDSIsuppression ' DSI'],'visible','off'); 
mesh(X*1e-3,Y,(rdmap_compansated')); axis tight
xlabel('Range (Km)')
ylabel('Doppler shift (Hz)')
title(['RDM After ' whichDSIsuppression ' DSI Suppression'])

% Beamformed Array ouput without Direct Path Suppression
analogBeam = sum(survChannelArray,2);
[rdmap, ranges, freqs] = rangedopplerfft(analogBeam,samplingFreq , 2*timeDelay*propSpeed , freqVector, NoisyrefSignal');
f5 = figure('Name',['Phased array DVB-S random after ' whichDSIsuppression ' without compensation'],'visible','off'); 
contourf(X*1e-3,Y,rdmap')
xlabel('Range (Km)')
ylabel('Doppler shift (Hz)')
title('DVB-S Beamformed No DSI-Suppression')

%save the variables I need for 2D CFAR
save('myData/CFAR_files','X','Y','rdmap_compansated');

%% angle estimation
estimator = phased.MUSICEstimator2D('SensorArray',array,...
    'OperatingFrequency',txSat(1).freq,...
    'NumSignalsSource','Property',...
    'DOAOutputPort',true,'NumSignals',2,...
    'AzimuthScanAngles',-60:.5:60,...
    'ElevationScanAngles',-20:.5:80);
[~,doas] = estimator(NewsurvChannelArray);
f6 = figure('visible','off'); 
plotSpectrum(estimator); axis tight
title('DOA MUSIC Spectrum');
xlabel('Azimuth')
ylabel('Elevation')



%% Bistatic parameter estimation
 [tx2tgEst,rx2tgEst] = paraExtraction(Ddiff, norm(tx2rx), theta1);
 
 
%% Save Figures - Whether using Windows or Mac/Linux
% 'Windows' = Windows, 'Mac' = Mac/Linux
% Note: Need to create folder myPlot in project directory
os = 'Windows';

% For Windows
if (strcmp(os,'Windows'))
    
    saveas(f1,'.\myPlot\DVB_S_Time_Signals.fig', 'fig')
    saveas(f2,'.\myPlot\Signal_Spectrum','fig');
    saveas(f3,['.\myPlot\Array_After_' whichDSIsuppression '_DSI'],'fig')
    saveas(f4,['.\myPlot\Array_' whichDSIsuppression '_DSI_3D'],'fig')
    saveas(f5,'.\myPlot\DVB_S_Beamformed_No_DSI','fig')
    saveas(f6,'.\myPlot\Music_Spectrum','fig');

% For Mac/Linux
elseif (strcmp(os,'Mac'))

    saveas(f1,'./myPlot/DVB_S_Time_Signals.fig', 'fig')
    saveas(f2,'./myPlot/Signal_Spectrum','fig');
    saveas(f3,['./myPlot/Array_After_' whichDSIsuppression '_DSI'],'fig')
    saveas(f4,['./myPlot/Array_' whichDSIsuppression '_DSI_3D'],'fig')
    saveas(f5,'./myPlot/DVB_S_Beamformed_No_DSI','fig')
    saveas(f6,'./myPlot/Music_Spectrum','fig');
    
end



