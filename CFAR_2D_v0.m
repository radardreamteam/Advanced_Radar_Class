%% 2D CFAR
tic;
close all; clear all;
load('myData/CFAR_files');

plotFigures();
%create the cfar object
cfar2D = phased.CFARDetector2D('Method','CA','GuardBandSize',10,'TrainingBandSize',80,'ProbabilityFalseAlarm',8e-3);


%% My works but doesn't detect anything
%sub in my own values do this later
resp=rdmap_compansated;
rngGrid=X*1e-3;
rngGrid=rngGrid(1,:)';
dopGrid=flipud(Y);
dopGrid=dopGrid(:,1);

%add in next step
[~,rangeIndx] = min(abs(rngGrid-[3 15]));
[~,dopplerIndx] = min(abs(dopGrid-[-320 300]));
[columnInds,rowInds] = meshgrid(dopplerIndx(1):dopplerIndx(2),rangeIndx(1):rangeIndx(2));
CUTIdx = [rowInds(:) columnInds(:)]';


%plots or whatever
%this is same as calling step function
%detections = cfar2D(resp,CUTIdx);
[detections]=step(cfar2D,resp,CUTIdx);
if_zero_nothing_detected=max(detections)
helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections,X,Y)

%time it
elapsed_time=toc()