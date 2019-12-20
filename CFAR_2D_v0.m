%% 2D CFAR
tic;
close all; clear all;
%load('myData/new/CFAR_files');

%plotFigures();
%Dimitri plotter function dimension reduces the cfar data
dimitri_plotter();
%create the cfar object
%if you would like to test against Matlabs version
%cfar2D = phased.CFARDetector2D('Method','GOCA','GuardBandSize',4,'TrainingBandSize',...
 %   12,'ProbabilityFalseAlarm',1e-5,'ThresholdOutputPort',true);
 
 % Always Specify GOCA even though it will be VI-CFAR
cfar2D = CFARDetector2D_dimitri('Method','GOCA','GuardBandSize',4,'TrainingBandSize',...
    12,'ProbabilityFalseAlarm',1e-6,'ThresholdOutputPort',true);

%% prepare grids for use with matlab adjusted CFAR algorithms
%sub in my own values do this later
reducing=true;
if reducing% reduction
    resp=rdmap_reduced;
    rngGrid=X_reduced*1e-3;%scale range grid
    rngGrid=rngGrid(1,:)';%grab only one dimension of the range grid
    dopGrid=flipud(Y_reduced);%invert doppler grid to fix scaling??
    dopGrid=dopGrid(:,1);%grab only one dimension
    
    [~,rangeIndx] = min(abs(rngGrid-[10 45]));
    [~,dopplerIndx] = min(abs(dopGrid-[-300 300]));
else
    resp=rdmap_compansated;
    rngGrid=X*1e-3;%scale range grid
    rngGrid=rngGrid(1,:)';%grab only one dimension of the range grid
    dopGrid=flipud(Y);%invert doppler grid to fix scaling??
    dopGrid=dopGrid(:,1);%grab only one dimension
    
    [~,rangeIndx] = min(abs(rngGrid-[15 40]));
    [~,dopplerIndx] = min(abs(dopGrid-[-400 400]));
end


%add in next step

[columnInds,rowInds] = meshgrid(dopplerIndx(1):dopplerIndx(2),rangeIndx(1):rangeIndx(2));
CUTIdx = [rowInds(:) columnInds(:)]';


%plots or whatever
%this is same as calling step function
%detections = cfar2D(resp,CUTIdx);
[detections,th]=step(cfar2D,resp,CUTIdx);%take a look at the computer threshold
%if_zero_nothing_detected=max(detections)

%add in optional masking
mask=true;
if mask
    %get sub sections of both thresholds and the RD map
    
end

if reducing
    helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections,th,X_reduced,Y_reduced,CUTIdx)
else
    helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections,th,X,Y)
end
%display max of threshold and dopplermap
doppler_max=max(max(rdmap_compansated))
thresh_max=max(max(th))

%time it
elapsed_time=toc()