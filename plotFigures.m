% new script for plotting stuff in myplot
close all;

operating_system='win';

if operating_system=='mac'
    openfig('myPlot/Array_After_NLMS_DSI.fig','visible');
%     openfig('myPlot/Array_NLMS_DSI_3D.fig','visible');
%     openfig('myPlot/DVB_S_Beamformed_No_DSI.fig','visible');
%     openfig('myPlot/DVB_S_Time_Signals.fig','visible');
%     openfig('myPlot/Music_Spectrum.fig','visible');
%     openfig('myPlot/Signal_Spectrum.fig','visible');
else %PC
    openfig('myPlot\Array_After_FBLMS_DSI.fig','visible');
    openfig('myPlot\Array_FBLMS_DSI_3D.fig','visible');
    openfig('myPlot\DVB_S_Beamformed_No_DSI.fig','visible');
    openfig('myPlot\DVB_S_Time_Signals.fig','visible');
    openfig('myPlot\Music_Spectrum.fig','visible');
    openfig('myPlot\Signal_Spectrum.fig','visible');
end