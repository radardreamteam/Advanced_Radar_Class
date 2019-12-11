function [DSC,value] = CLEAN(Surv,Ref,K,d,shiftMax,rdmap,ranges,freqs,C,samplingFreq,sub_algorithm)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%DSC = (eye(T) - X*inv(X'*X)*X')*Surv';
% A = X'*X;
% A = inv(A);
% A = X*A;
% A = A*X';
% DSC = (eye(T) - A)*Surv';
j = 1i;
index = [0 0];
Bandwidth = 12;
clutter_band = rdmap(:,(size(rdmap,2)-1)/2-Bandwidth+1:(size(rdmap,2)-1)/2+Bandwidth+1);
[value1, ind1] = max(clutter_band,[],1);
[value, index(2)] = max(value1);
index(1) = ind1(index(2));
peakDelay = round(ranges(index(1))/C*samplingFreq);
peakDoppler = freqs(index(2)+(size(rdmap,2)-1)/2-Bandwidth);

T = length(Ref);
X = zeros(T,1);

taxis       = (0:(T-1))/samplingFreq;
col = 1;
for k = -shiftMax : 1 : shiftMax
    temp_Ref = Ref.*(exp(j*2*pi*(peakDoppler+k)*taxis(1:T)));
    for i = -K: 1: K
        if(peakDelay ~= 0 || i>=0)
            temp = zeros(1,peakDelay+d*i);
            X(:,col) = [temp temp_Ref(1:T-length(temp))]';
            col = col + 1;
        end
    end
end

switch sub_algorithm
%ECA for cancellation
    case 'ECA'
          DSC = Surv - X*((X'*X)\X'*Surv);
    case 'NLMS'
          filterOrder = 32;
          DSIed = Surv;
          nlms = dsp.LMSFilter(filterOrder,'Method','Normalized LMS','StepSizeSource','Input port');
          for i = 1:size(X,2)
              nlms.reset();
              [~,errNLMS,~] = nlms(X(:,i),DSIed,0.001);
              DSIed = errNLMS;
          end
          DSC = DSIed;
    case 'Wiener'
          filterOrder = 32;
          DSIed = Surv;
          for i = 1:size(X,2)
              [~ ,errW] = wienerf(filterOrder-1,X(:,i),DSIed);
              DSIed = errW;
          end
          DSC = DSIed;
    case 'FBLMS'
          filterLength = 512;
          DSIed = Surv;
          [numRow,numCol] = size(X);
          remData = mod(numRow,filterLength);
          truncatedData = DSIed(1:(numRow - remData),:);
          directDistData = X(1:(numRow - remData),:);
          fdaf = dsp.FrequencyDomainAdaptiveFilter('Length',filterLength,'BlockLength',filterLength,'StepSize',0.001);
          fdaf.reset()
          for i = 1:size(X,2)
            [~,errFBLMS] = fdaf(directDistData(:,i),truncatedData);
            DSIed = errFBLMS;
          end
          DSC = DSIed;
    case 'RLS'
        filterOrder = 32;
        DSIed = Surv;
        FF = 0.98; %Forgetting Factor
        IP = 1; %Initialization Parameter
        for i = 1:size(X,2)
            [~,errRLS] = RLSFilter(filterOrder,X(:,i),DSIed,FF,IP);
            DSIed = errRLS;
        end
        DSC = DSIed;
    otherwise
        disp('Not a valid algorithm for clean');
end
return