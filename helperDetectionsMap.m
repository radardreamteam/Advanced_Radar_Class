function helperDetectionsMap(resp,rng_grid,dop_grid,rangeIndx,dopplerIndx,detections,th,X,Y,CUTIdx)
% This function is only in support of CFARDetectionExample. It may be
% removed in a future release.

%   Copyright 2016 The MathWorks, Inc.

%it has been adapted for VI-CFAR by Dimtri Herr

%detections map
%make detection map
detectionMap = zeros(size(resp));
 detectionMap(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) = ...
   reshape(double(detections),rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);


%plot the thresholds here
%reshape thresholds like we do with detections in helperDetectionsMap
threshMap=zeros(size(resp));
threshMap(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) =...
    reshape(th,rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);



%active normalize adaptive threshold masking
%take sub section that is threshmap of actual area
%take subsection of thresholds
sub_thresh=reshape(th,rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);%this is the thresholds
%take subsection of rdmap
sub_map=resp(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2));
%now I need subsection of the X and Y axises ugh
sub_X=X(dopplerIndx(1):dopplerIndx(2),rangeIndx(1):rangeIndx(2));
sub_Y=Y(dopplerIndx(1):dopplerIndx(2),rangeIndx(1):rangeIndx(2));

%compute normalization mask
new_sub_map=sub_map.*normalize(sub_thresh);

%recompute detections
new_det=zeros(size(new_sub_map));
new_det(find(new_sub_map>sub_thresh))=1;
if max(max(new_det))~=1
    disp("no detections");
else
    disp("object detected");
end

%plot subsections to confirm viability of simulation parameters
local=false;
if ~local
    f7=figure('visible','off');
else
    f7=figure();
end

m2=mesh(sub_X*1e-3,sub_Y,(new_sub_map'), 'edgecolor', 'k'); axis tight
hold on;
m1=mesh(sub_X*1e-3,sub_Y,(sub_thresh'));
alpha(m1,0.1);

title("3D mesh of sub thresholds colored(thresholds) vs. data(black)");
xlabel('Range (Km)')
ylabel('Doppler shift (Hz)')

%plot the detections of new thing
if ~local
    f8=figure('visible','off');
else
    f8=figure();
end
contourf(sub_X'*1e-3,sub_Y',new_sub_map)
hold on;
h2 = imagesc(sub_X(1,:)*1e-3,sub_Y(:,1),new_det', 'AlphaData', .2);
ylabel('Doppler (Hz)'); xlabel('Range (Km)'); title('CFAR Detections overlayed on RDM');

h2.Parent.YDir = 'normal';
hold off;

% save the plots 
saveas(f7,'./myPlot/CFAR_Mesh.fig', 'fig')
saveas(f8,'./myPlot/CFAR_Detection.fig','fig');

