%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generate the DVB-s signal
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
% Set clutter
clutter_contamination = 0;
% Initial System parameters from Balke Sheet
directSigPower          = -122.94+30;
echoSigPower            = -127; 
refGain                 = 20;
% The goal is to eliminate this attenaution via signal processing
dirPathAttenuation      = -20;%dB, Null

%% System Parameters
SystemParameters;
% Update System parameters
echoSigPower            = tgtPower_dBm;
directSigPower          = directPathPower_dBm;

%% Clutter Parameters
clutterSigPower = [directSigPower - 10, directSigPower - 5,directSigPower - 20];
clutterSig_Delay = [200, 300, 400];
clutterSig_Doppler = [5,7,10];
clutterSig_Doa = [70 79 69 ;51 43 61];

%% Choose Methods
% Choose operating methods
% 1 - beamform, 2 - element
whichCancellation = 'element';
% Choose the Direct Signal Suppression Technique
% 1 - NLMS, 2 - Wiener, 3 - RLS, 4 - FBLMS, 5 - ECA, 6 - SCA, 7 - CLEAN
whichDSIsuppression = 'FBLMS';
% Choose a sub-cancellation algorithm for CLEAN
sub_algorithm = 'Wiener';

%% Create input signal
sigNumber   = 400;
numOfFrames = 3;
pilotOn     = 1; % insert pilot module
maxSymbol   = 5000*4; % Max symbol number we generate, optimized for beamform
integrationTime = 4e-3; % second
BaseBandSignal  = generate_DVBS(integrationTime,3e6,samplingFreq,maxSymbol,numOfFrames,pilotOn);
sigLength   = length(BaseBandSignal);
dirPath     = zeros(1, sigLength + max(samp_offset));
indirPath   = zeros(1, sigLength + max(samp_offset));
taxis       = (0:(length(dirPath)-1))/samplingFreq;
Echo        = cell(1,3);
for kk = 1:numTargets
    Echo{kk} = zeros(1,sigLength + max(samp_offset));
end

%% Create Input Signals
switch whichCancellation
    case 'element'
        elementCancellation = 1;
        beamformedCancellation = 0;
        survAntGain             = 32; %only needed in element-by-element
       %% Make Time Domain Signal
        dirPath(1:sigLength) = BaseBandSignal;
        for k3 = 1:numTargets
            Echo{k3}(samp_offset(k3)+1:samp_offset(k3)+sigLength) = dirPath(1:sigLength)*exp(j*phaseOffset).*exp(j*2*pi*FShift(k3)*taxis(1:sigLength));
        end
        echoSigPower                       = echoSigPower + survAntGain;
        for k0 = 1:numTargets
            indirPath(k0,:)  = sqrt(10^(echoSigPower(k0)/ 10))*Echo{k0};
        end
        %The direct signal from the Tx should arrive at both Rx at the same time,there is no shift.
        %make the reference signal
        refSignal                          = sqrt(10^((directSigPower+refGain)/10))*dirPath;
        %Create surveillance channel
        echoSignal                         = indirPath;
        directSignal                       = 10^(dirPathAttenuation/10)*10^(directSigPower/10)*dirPath;
        survChannel                        = sum(echoSignal) + directSignal ;
        %Create Noise 
        noiseSigSurv     = sqrt(10^(noisePower_dBm/10))*(randn(1,length(dirPath))+j*randn(1,length(dirPath)));
        survNoisyChannel = survChannel + noiseSigSurv;

        noiseSigRef     = sqrt(10^(noisePower_dBm/10))*(randn(1,length(dirPath))+j*randn(1,length(dirPath)));
        NoisyrefSignal  = refSignal + noiseSigRef;
       %% Estimate DOA using MUSIC
        % Build a 25*25 rectangular array, and use MUSIC algorithm to seperate these two signals 
        array = phased.URA('Size',[25 25],'ElementSpacing',[lamda/2 lamda/2]);
        array.Element.FrequencyRange = [500.0e6 5000.0e6];
        % Target's angle
        doaT1 = [10;10];
        doaT2 = [10;10];
        doaT3 = [10;10];
        % Given the direct signal in the surveillance channel is not helping, we
        % want to place it as far away possible in the sidelobes.
        doaTx = [80; -80]; %[30;-20];
        survChannelArray  = collectPlaneWave(array,[echoSignal(1,:)',echoSignal(2,:)',echoSignal(3,:)',directSignal'],[doaT1,doaT2,doaT3,doaTx],txSat(1).freq); %samplingFreq);
        noiseArray        = sqrt(10^((noisePower_dBm)/10))*(randn(length(dirPath),625)+j*randn(length(dirPath),625));
        survChannelArray  = survChannelArray + noiseArray ;
        % We must assume we know the doa of the direct path signal
        direct_distribute = collectPlaneWave(array,directSignal',doaTx,samplingFreq);
        NewsurvChannelArray = zeros(size(survChannelArray,1),size(survChannelArray,2));
    case 'beamform'
        beamformedCancellation = 1;
        elementCancellation = 0;
       %% Make Time Domain Signal
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
        %Generate Surveillence Channel Array, Beamfored Surveillence %Channel, and Direct Reference Channel
        direct_ref      = sqrt(10^((directSigPower+refGain)/10))*dirPath;
        noise           = sqrt(10^(noisePower_dBm/10))*(randn(1,length(dirPath))+j*randn(1,length(dirPath)));
        direct_ref = direct_ref + noise;
        direct_reft = transpose(direct_ref);
        survChannelArray  = collectPlaneWave(array,Echos,doas,txSat(1).freq,propSpeed); %samplingFreq);
        noiseArray        = sqrt(10^((noisePower_dBm)/10))*((randn(length(dirPath),elementNum)+j*randn(length(dirPath),elementNum)));
        survChannelArray  = survChannelArray + noiseArray;
        clear noiseArray;
       %% Determine Beam Direction
        beamformer = phased.PhaseShiftBeamformer('SensorArray',array,'PropagationSpeed',propSpeed,'OperatingFrequency',txSat(1).freq,'DirectionSource','Input port','WeightsOutputPort',true);
        ANG = doas(:,4);
        [~,W] = beamformer(survChannelArray, ANG);
        [PAT,AZ_ANG,EL_ANG] = pattern(array,txSat(1).freq,[-80:1:80],[-80:1:80],'CoordinateSystem','rectangular','Type','directivity','PropagationSpeed',propSpeed,'Weight',W);
        doas = round(doas);
        BeamformedSurv = sqrt(10^(PAT(doas(2,1)+21,doas(1,1)+81)/10))*refSignal';
        for k3 = 2:numTargets+1
            BeamformedSurv = BeamformedSurv+sqrt(10^(PAT(doas(2,k3)+21,doas(1,k3)+81)/10))*Echos(:,k3);
        end
    otherwise
    disp('Wrong Cancellation')
end

%% Add Clutter
if(clutter_contamination == 1)
    for k4 = 1:length(clutterSigPower)
        clutter_temp = zeros(1,sigLength + max(samp_offset));
        clutter_temp(clutterSig_Delay(k4)+1:clutterSig_Delay(k4)+sigLength) = ...
            sqrt(10^(((clutterSigPower(k4)+PAT(clutterSig_Doa(2,k4)+21,clutterSig_Doa(1,k4)+81)))/10))*dirPath(1:sigLength)*exp(j*phaseOffset).*exp(j*2*pi*clutterSig_Doppler(k4)*taxis(1:sigLength));
        switch whichCancellation
            case 'beamform'
                %WIP
                BeamformedSurv = BeamformedSurv + clutter_temp';
            case 'element'
                survChannelArray = survChannelArray + clutter_temp';
        end
    end
    clear clutter_temp;
end
clear Echos Echo noise;
%% Spectrum Plots
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

%% DIRECT SIGNAL SUPPRESSION

switch whichDSIsuppression
    case 'NLMS'
        filterOrder = 32;
        if (elementCancellation == 1)
            parfor i = 1:size(survChannelArray,2)
                nlms = dsp.LMSFilter(filterOrder,'Method','Normalized LMS','StepSizeSource','Input port');
                nlms.reset();
                [outputNLMS,errNLMS,weightsNLMS] = nlms(direct_distribute(:,i),survChannelArray(:,i),0.001);
                NewsurvChannelArray(:,i) = errNLMS;
            end
        end
        if (beamformedCancellation == 1)
          nlms = dsp.LMSFilter(filterOrder,'Method','Normalized LMS','StepSizeSource','Input port');
          nlms.reset();
          [outputNLMS,errNLMS,weightsNLMS] = nlms(direct_ref',BeamformedSurv,0.001);
          DSIed = errNLMS;
        end
    case 'Wiener'
        filterOrder = 32;
        if (elementCancellation == 1)
            parfor i = 1:size(survChannelArray,2)
                [outputW ,errW] = wienerf(filterOrder-1,direct_distribute(:,i),survChannelArray(:,i));
                NewsurvChannelArray(:,i) = errW;
            end
        end
        if (beamformedCancellation == 1)
          [outputW ,errW] = wienerf(filterOrder-1,conj(direct_reft),BeamformedSurv);
          DSIed = errW;
        end
    case 'RLS'
        filterOrder = 32;
        FF = 0.98; %Forgetting Factor
        IP = 1; %Initialization Parameter
        if (elementCancellation == 1)
            parfor i = 1:size(survChannelArray,2)
                [outputRLS,errRLS] = RLSFilter(filterOrder,direct_distribute(:,i),survChannelArray(:,i),FF,IP);
                NewsurvChannelArray(:,i) = errRLS;
            end
        end
        if (beamformedCancellation == 1)
          [outputRLS ,errRLS] = RLSFilter(filterOrder-1,direct_ref',BeamformedSurv,FF,IP);
          DSIed = errRLS;
        end
    case 'FBLMS'
        filterLength = 768;
        if (elementCancellation == 1)
            %truncate some data to make blocklength a factor of number of data
            [numRow,numCol] = size(survChannelArray);
            remData = mod(numRow,filterLength);
            truncatedData = survChannelArray(1:(numRow - remData),:);
            directDistData = direct_distribute(1:(numRow - remData),:);
            %accommodate NewsurvChannelArray pre-allocation to above
            NewsurvChannelArray = zeros(size(truncatedData,1),size(truncatedData,2));
            parfor i = 1:size(truncatedData,2)            
                fdaf = dsp.FrequencyDomainAdaptiveFilter('Length',filterLength,'BlockLength',filterLength,'StepSize',1);
                fdaf.reset()
                [outputFBLMS,errFBLMS] = fdaf(directDistData(:,i),truncatedData(:,i));
                NewsurvChannelArray(:,i) = errFBLMS;
            end
        end
        if (beamformedCancellation == 1)
            %truncate some data to make blocklength a factor of number of data
            [numRow,numCol] = size(BeamformedSurv);
            remData = mod(numRow,filterLength);
            truncatedData = BeamformedSurv(1:(numRow - remData));            
            directDistData = direct_ref(1:(numRow - remData));
            %accommodate NewsurvChannelArray pre-allocation to above
            fdaf = dsp.FrequencyDomainAdaptiveFilter('Length',filterLength,'BlockLength',filterLength,'StepSize',0.001);
            fdaf.reset()
            [outputFBLMS,errFBLMS] = fdaf(directDistData(:),truncatedData);
            DSIed = errFBLMS;
        end
    case 'ECA'
        K = 25;
        d = 1;
        if (beamformedCancellation == 1)
            DSIed = ECA(BeamformedSurv', direct_ref,K,d);
        end
        if (elementCancellation == 1)
                b = 10; %number of batches
                [numRow,numCol] = size(survChannelArray);
                remData = mod(numRow,b);
                truncatedData = survChannelArray(1:(numRow - remData),:);
                directDistData = direct_distribute(1:(numRow - remData),:);
                Ns = size(truncatedData,1);
                Nr = size(directDistData,1);
            for p=1:b
                bIdx_s=(p-1)*Ns/b+1:p*Ns/b; %batch indices for survillence
                bIdx_r=(p-1)*Nr/b+1:p*Nr/b; %batch indices for reference
                survChannelArrayBatch= truncatedData(bIdx_s,:);
                refChannelArrayBatch= directDistData(bIdx_r,:);
                NewsurvChannelArrayBatch = zeros(size(survChannelArrayBatch,1),size(survChannelArrayBatch,2));
                parfor i = 1:size(survChannelArrayBatch,2)        
                    NewsurvChannelArrayBatch(:,i) = ECA(survChannelArrayBatch(:,i)',refChannelArrayBatch(:,i)',K,d);
                end
                NewsurvChannelArray(bIdx_s,:) = NewsurvChannelArrayBatch;
            end
        end
    case 'SCA'
        % WIP
        K = 25;
        d = 1;
        M = 10;
        if (beamformedCancellation == 1)
            DSIed = SCA(BeamformedSurv', direct_ref,K,d,M);
        end
        if (elementCancellation == 1)
            b = 10; %batch number
            [numRow,numCol] = size(survChannelArray);
            remData = mod(numRow,b);
            truncatedData = survChannelArray(1:(numRow - remData),:);
            directDistData = direct_distribute(1:(numRow - remData),:);
            Ns = size(truncatedData,1);
            Nr = size(directDistData,1);
            for p=1:b
                bIdx_s=(p-1)*Ns/b+1:p*Ns/b;
                bIdx_r=(p-1)*Nr/b+1:p*Nr/b;
                survChannelArrayBatch= truncatedData(bIdx_s,:);
                refChannelArrayBatch= directDistData(bIdx_r,:);
                NewsurvChannelArrayBatch = zeros(size(survChannelArrayBatch,1),size(survChannelArrayBatch,2));
                parfor i = 1:size(survChannelArrayBatch,2)       
                    NewsurvChannelArrayBatch(:,i) = SCA(survChannelArrayBatch(:,i)',refChannelArrayBatch(:,i)',K,d, M);
                end
                NewsurvChannelArray(bIdx_s,:) = NewsurvChannelArrayBatch;
            end
        end
    case 'CLEAN'
            K = 3;
            d = 1;
            shiftMax = 3;
            maxIter = 8;
            % Initial Value
            iter = 0;
            [rdmap, ranges, freqs] = rangedopplerfft(BeamformedSurv,samplingFreq , 2*max(timeDelay)*propSpeed , freqVector, direct_ref');
            [X,Y] = meshgrid(ranges, freqs);
            [DSIed,maxscore] = CLEAN(BeamformedSurv,direct_ref,K,d,shiftMax,rdmap,ranges,freqs,propSpeed,samplingFreq,sub_algorithm);
            while(maxscore > 5e-06 && iter <= maxIter)
                [rdmap, ranges, freqs] = rangedopplerfft(DSIed,samplingFreq , 2*max(timeDelay)*propSpeed , freqVector, direct_ref');
                iter = iter + 1;
                [DSIed, maxscore] = CLEAN(DSIed,direct_ref,K,d,shiftMax,rdmap,ranges,freqs,propSpeed,samplingFreq,sub_algorithm);
            end
    otherwise
        disp('Not avalid DSI Algorithm')
end

%% Range-Doppler MAP 
switch whichCancellation
    case 'element'        
        %compensate the phase shift to each element
        reference_distribute = collectPlaneWave(array,NoisyrefSignal',[0; 0],samplingFreq);
        % multi_channel_rdmap = zeros(size(rdmap,1),size(rdmap,2),size(reference_distribute,2));
        for i = size(reference_distribute,2)
            [rdmap, ranges, freqs] = rangedopplerfft(NewsurvChannelArray(:,i),samplingFreq,1.3*max(timeDelay)*propSpeed , freqVector, reference_distribute(:,i));
            multi_channel_rdmap(:,:,i) = rdmap;
        end
        rdmap_compansated = sum(multi_channel_rdmap,3);
    case 'beamform'
        [rdmap_compansated, ranges, freqs] = rangedopplerfft(DSIed,samplingFreq, 2*max(timeDelay)*propSpeed , freqVector, direct_ref');
end
[X,Y] = meshgrid(ranges, freqs);


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
switch whichCancellation
    case 'element'  
        analogBeam = sum(survChannelArray,2);
        [rdmap, ranges, freqs] = rangedopplerfft(analogBeam,samplingFreq , 1.3*max(timeDelay)*propSpeed , freqVector, NoisyrefSignal');
    case 'beamform'
        [rdmap, ranges, freqs] = rangedopplerfft(BeamformedSurv,samplingFreq , 2*max(timeDelay)*propSpeed , freqVector, direct_ref');
end
f5 = figure('Name',['Phased array DVB-S random after ' whichDSIsuppression ' without compensation'],'visible','off'); 
contourf(X*1e-3,Y,rdmap')
xlabel('Range (Km)')
ylabel('Doppler shift (Hz)')
title('DVB-S Beamformed No DSI-Suppression')

%save the variables I need for 2D CFAR
save('myData/CFAR_files','X','Y','rdmap_compansated');

% %% angle estimation
% estimator = phased.MUSICEstimator2D('SensorArray',array,...
%     'OperatingFrequency',txSat(1).freq,...
%     'NumSignalsSource','Property',...
%     'DOAOutputPort',true,'NumSignals',2,...
%     'AzimuthScanAngles',-60:.5:60,...
%     'ElevationScanAngles',-20:.5:80);
% [~,doas] = estimator(NewsurvChannelArray);
% f6 = figure('visible','off'); 
% plotSpectrum(estimator); axis tight
% title('DOA MUSIC Spectrum');
% xlabel('Azimuth')
% ylabel('Elevation')



%% Bistatic parameter estimation
% need fix SystemParameters about theta1
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
    %saveas(f6,'.\myPlot\Music_Spectrum','fig');

% For Mac/Linux
elseif (strcmp(os,'Mac'))

    saveas(f1,'./myPlot/DVB_S_Time_Signals.fig', 'fig')
    saveas(f2,'./myPlot/Signal_Spectrum','fig');
    saveas(f3,['./myPlot/Array_After_' whichDSIsuppression '_DSI'],'fig')
    saveas(f4,['./myPlot/Array_' whichDSIsuppression '_DSI_3D'],'fig')
    saveas(f5,'./myPlot/DVB_S_Beamformed_No_DSI','fig')
    %saveas(f6,'./myPlot/Music_Spectrum','fig');
    
end



