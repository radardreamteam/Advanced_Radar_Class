function helperDetectionsMap(resp,rng_grid,dop_grid,rangeIndx,dopplerIndx,detections,X,Y)
% This function is only in support of CFARDetectionExample. It may be
% removed in a future release.

%   Copyright 2016 The MathWorks, Inc.

%figure
%subplot(2,1,1);
detectionMap = zeros(size(resp));
detectionMap(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) = ...
  reshape(double(detections),rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);
%h = imagesc(dop_grid,rng_grid,detectionMap);
%xlabel('Doppler (Hz)'); ylabel('Range (m)'); title('Range Doppler CFAR Detections');
%h.Parent.YDir = 'normal';

% subplot(2,1,2);
figure
contourf(X'*1e-3,Y',resp)
hold on;
h2 = imagesc(rng_grid,flipud(dop_grid),detectionMap');
ylabel('Doppler (Hz)'); xlabel('Range (m)'); title('Range Doppler CFAR Detections axis flipped');

h2.Parent.YDir = 'normal';
hold off;