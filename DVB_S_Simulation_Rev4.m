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
elementCancellation     = 0;
beamformedCancellation  = 1;
% Initial System parameters from Balke Sheet
directSigPower          = -122.94+30;
echoSigPower            = -127; 
refGain                 = 20;
% The goal is to eliminate this attenaution via signal processing
dirPathAttenuation      = 0;%dB, Null
%survAntGain             = 32; %40;

%% System Parameters
SystemParameters;
% Update System parameters
echoSigPower            = tgtPower_dBm;
directSigPower          = directPathPower_dBm;

%% Choose the Direct Signal Suppression Technique
% 1 - NLMS, 2 - Wiener, 3 - RLS, 4 - FBLMS, 5 - ECA, 6 - SCA
whichDSIsuppression = 'Wiener';

%% Create input signal
sigNumber   = 400;

numOfFrames = 3;
pilotOn     = 1; % insert pilot module
maxSymbol   = 5000*4; % Max symbol number we generate
integrationTime = 1e-3; % second
BaseBandSignal  = generate_DVBS(integrationTime,3e6,samplingFreq,maxSymbol,numOfFrames,pilotOn);
sigLength   = length(BaseBandSignal);
dirPath     = zeros(1, sigLength + max(samp_offset));
taxis       = (0:(length(dirPath)-1))/samplingFreq;
Echo        = cell(1,3);
for kk = 1:numTargets
    Echo{kk} = zeros(1,sigLength + max(samp_offset));
end

%% Time domain signal
dirPath(1:sigLength) = BaseBandSignal;
refSignal            = sqrt(10^((directSigPower+dirPathAttenuation)/10))*dirPath;
doas = uv2azel([uDirectPath;vDirectPath]);
Echos= refSignal';
for k3 = 1:numTargets
    Echo{k3}(samp_offset(k3)+1:samp_offset(k3)+sigLength) = ...
        sqrt(10^((tgtPower_dBm(k3))/10))*dirPath(1:sigLength)*exp(j*phaseOffset).*exp(j*2*pi*FShift(k3)*taxis(1:sigLength));
    doas = [doas uv2azel([uTarget(k3);vTarget(k3)])];
    Echos = [Echos Echo{k3}'];
end
array = phased.URA('Size',[rx.numRow  rx.numCol],'ElementSpacing',[lamda/2 lamda/2]);
elementNum = rx.numRow *rx.numCol;
array.Element.FrequencyRange = [1.5e9 3.0e9];
%% Generate survChannelArray, BeamformedSurv,direct_ref 
direct_ref      = sqrt(10^((directSigPower+refGain)/10))*dirPath;
noise           = sqrt(10^(noisePower_dBm/10))*(randn(1,length(dirPath))+j*randn(1,length(dirPath)));
direct_ref = direct_ref + noise;
survChannelArray  = collectPlaneWave(array,Echos,doas,txSat(1).freq,propSpeed); %samplingFreq);

noiseArray        = sqrt(10^((noisePower_dBm)/10))*((randn(length(dirPath),elementNum)+j*randn(length(dirPath),elementNum)));
survChannelArray  = survChannelArray + noiseArray ;
clear noiseArray;
beamformer = phased.PhaseShiftBeamformer('SensorArray',array,'PropagationSpeed',propSpeed,'OperatingFrequency',txSat(1).freq,'DirectionSource','Input port','WeightsOutputPort',true);
%% Define your beam direction
ANG = doas(:,2);
[~,W] = beamformer(survChannelArray, ANG);

[PAT,AZ_ANG,EL_ANG] = pattern(array,txSat(1).freq,[-80:1:80],[-20:1:80],'CoordinateSystem','rectangular','Type','directivity','PropagationSpeed',propSpeed,'Weight',W);
%tempEcho(samp_offset+1:end)        = indirPath(samp_offset+1:end); 
% BeamformedSurv = sum(Echos,2);
doas = round(doas);
BeamformedSurv = sqrt(10^(PAT(doas(2,1)+21,doas(1,1)+81)/10))*refSignal';
for k3 = 2:numTargets+1
    BeamformedSurv = BeamformedSurv+sqrt(10^(PAT(doas(2,k3)+21,doas(1,k3)+81)/10))*Echos(:,k3);
end
clear Echos Echo noise;
%The direct signal from the Tx should arrive at both Rx at the same time,
%there is no shift.
%% Do a spectrum
% Generate Doppler Vector to simulate
freqVector      = -500:2:500;

% Frequency domain signal
dirPath_fft = fftshift(abs(fft(dirPath))/length(dirPath));
%indirPath_fft = fftshift(abs(fft(indirPath))/length(indirPath));

%%% Csingle Channel Plots
% axes
taxis_sig = (0:(length(dirPath)-1))/samplingFreq;
faxis_sig = linspace(-samplingFreq/2,samplingFreq/2,length(dirPath));

% Time domain of signals
f1 = figure('visible','off'); 
%hold on; 
plot(taxis_sig*1e6, real(dirPath));
%hold on; plot(taxis_sig*1e6, real(indirPath)); axis tight; grid
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
title('Time Domain DVB-S Single Channel Signal')
legend('Direct Path')%,'Indirect Path')

% Freq domain of signals
f2 = figure('visible','off'); 
indFreqSide = find((faxis_sig >= 0),1);

hold on; plot(faxis_sig(indFreqSide:end)*1e-6, 20*log10(dirPath_fft(indFreqSide:end)));

axis([0 6 -120 5]); grid
xlabel('Frequency (MHz)')
ylabel('Amplitude (dBm)')
title('DVB-S Single signal Spectrum')

%% DIRECT SIGNAL SUPPRESSION
% NLMS Implementation on Multi-Channel
% We must assume we know the doa of the direct path signal
direct_distribute = collectPlaneWave(array,direct_ref',doas(:,1),samplingFreq);
% noiseArray        = sqrt(10^(noisePower_dBm/10))*(randn(length(tempDir),625)+j*randn(length(tempDir),625));
% direct_distribute = direct_distribute + noiseArray;
NewsurvChannelArray = zeros(size(survChannelArray,1),size(survChannelArray,2));
switch whichDSIsuppression
    case 'NLMS'
        if (elementCancellation == 1)
            parfor i = 1:size(survChannelArray,2)
                filterOrder = 32;
                nlms = dsp.LMSFilter(filterOrder,'Method','Normalized LMS','StepSizeSource','Input port');
                nlms.reset();
                [outputNLMS,errNLMS,weightsNLMS] = nlms(direct_distribute(:,i),survChannelArray(:,i),0.001);
                NewsurvChannelArray(:,i) = errNLMS;
            end
        end
        if (beamformedCancellation == 1)
          filterOrder = 32;
          nlms = dsp.LMSFilter(filterOrder,'Method','Normalized LMS','StepSizeSource','Input port');
          nlms.reset();
          [outputNLMS,errNLMS,weightsNLMS] = nlms(direct_ref',BeamformedSurv,0.001);
          DSIed = errNLMS;
        end
    case 'Wiener'
        if (elementCancellation == 1)
            parfor i = 1:size(survChannelArray,2)
                filterOrder = 32;
                [outputW ,errW] = wienerf(filterOrder-1,direct_distribute(:,i),survChannelArray(:,i));
                NewsurvChannelArray(:,i) = errW;
            end
        end
        if (beamformedCancellation == 1)
          filterOrder = 32;
          [outputW ,errW] = wienerf(filterOrder-1,direct_ref',BeamformedSurv);
          DSIed = errW;
        end
    case 'RLS'
        
    case 'FBLMS'
        %Change the 'BlockLength' parameter to interger divisible of the filter inputs
        if (elementCancellation == 1)
            parfor i = 1:size(survChannelArray,2)
                filterOrder = 32;
                fdaf = dsp.FrequencyDomainAdaptiveFilter(filterOrder,'BlockLength',31,'StepSize',0.001);
                [outputFBLMS,errFBLMS] = fdaf(direct_distribute(:,i),survChannelArray(:,i));
                NewsurvChannelArray(:,i) = errFBLMS;
            end
        end
        if (beamformedCancellation == 1)
            filterOrder = 32;
            fdaf = dsp.FrequencyDomainAdaptiveFilter(filterOrder,'BlockLength',31,'StepSize',0.001);
            [outputFBLMS,errFBLMS] = fdaf(direct_ref',BeamformedSurv);
            DSIed = errFBLMS;
        end
    case 'ECA'
        if (beamformedCancellation == 1)
            K = 25;
            d =1;
            DSIed = ECA(BeamformedSurv', direct_ref,25,1);
        end
        if (elementCancellation == 1)
            
        for i = 1:size(survChannelArray,2)       
            NewsurvChannelArray(:,i) = ECA(survChannelArray(:,i)',direct_distribute(:,i)',K,d);
        end
        end
    case 'SCA'
        
    otherwise
        disp('Not avalid DSI Algorithm')
end

%% Range-Doppler MAP Array elements
%compensate the phase shift to each element
%reference_distribute = collectPlaneWave(array,NoisyrefSignal',doa2,samplingFreq);
if (elementCancellation == 1)
    multi_channel_rdmap = zeros(size(rdmap,1),size(rdmap,2),size(reference_distribute,2));
    for i = size(direct_distribute,2)
        [rdmap, ranges, freqs] = rangedopplerfft(NewsurvChannelArray(:,i),samplingFreq ,...
            1.3*max(timeDelay)*propSpeed , freqVector, direct_distribute(:,i));
        multi_channel_rdmap(:,:,i) = rdmap;
    end
    DSIed = sum(multi_channel_rdmap,3);
end
%% Range-Doppler MAP Array Beamformed
[rdmap_compansated, ranges, freqs] = rangedopplerfft(DSIed,samplingFreq , 2*max(timeDelay)*propSpeed , freqVector, direct_ref');
[X,Y] = meshgrid(ranges, freqs);
%rdmap_compansated = sum(multi_channel_rdmap,3);

f3 = figure('Name',['Phased array DVB-S After ' whichDSIsuppression ' DSI'],'visible','on'); 
contourf(X*1e-3,Y,rdmap_compansated')
xlabel('Range (Km)')
ylabel('Doppler shift (Hz)')
title(['RDM After ' whichDSIsuppression ' DSI Suppression'])

f4 = figure('Name',['Phased array DVB-S After ' whichDSIsuppression ' DSI'],'visible','on'); 
rd_db = 10*log10((rdmap_compansated'));
%rd_db = rd_db-max(rd_db, [],'all');
mesh(X*1e-3,Y,rd_db); 
axis([0 1 0 1 -120 0]);
axis 'auto xy';

xlabel('Range (Km)')
ylabel('Doppler shift (Hz)')
title(['RDM After ' whichDSIsuppression ' DSI Suppression 3D'])

% Beamformed Array ouput without Direct Path 
[rdmap, ranges, freqs] = rangedopplerfft(BeamformedSurv,samplingFreq , 2*max(timeDelay)*propSpeed , freqVector, direct_ref');
f5 = figure('Name',['Phased array DVB-S random after ' whichDSIsuppression ' without compensation'],'visible','on'); 
contourf(X*1e-3,Y,rdmap')
xlabel('Range (Km)')
ylabel('Doppler shift (Hz)')
title('DVB-S Beamformed No DSI-Suppression')


%% Implement CFAR


%save the variables I need for 2D CFAR
save('myData/CFAR_files','X','Y','rdmap_compansated');

%% angle estimation
estimator = phased.MUSICEstimator2D('SensorArray',array,...
    'OperatingFrequency',txSat(1).freq,...
    'NumSignalsSource','Property',...
    'DOAOutputPort',true,'NumSignals',4,...
    'AzimuthScanAngles',-60:.5:60,...
    'ElevationScanAngles',-20:.5:80);
[~,doas] = estimator(survChannelArray);
f6 = figure('visible','on'); 
plotSpectrum(estimator); axis tight
title('DOA MUSIC Spectrum');
xlabel('Azimuth')
ylabel('Elevation')



%% Bistatic parameter estimation
% [tx2tgEst,rx2tgEst] = paraExtraction(Ddiff, norm(tx2rx), theta1);
 
 
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



