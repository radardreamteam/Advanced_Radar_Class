%% System Parameters
N = 2; % number of independent sources
subspaceplot = 1021; % subspace to plot
% Transmitters Information
txSat(1).pos_lla   = [0, -85, 35686e3];
txSat(1).pos_ecef  = lla2ecef(txSat(1).pos_lla);
txSat(2).pos_lla   = [0, -85, 53686];
txSat(2).pos_ecef  = lla2ecef(txSat(2).pos_lla);
%-----------------------------------------------------
txSat(1).erp_dBm   = 58.45 + 30; %dBW to dBm
txSat(1).freq      = 2335.305e6;
lamda              = propSpeed/txSat(1).freq;
txSat(1).lambda    = propSpeed/txSat(1).freq;
txSat(2).erp_dBm   = 58.45 + 30; %dBW to dBm
txSat(2).freq      = 2333.405e6;
txSat(2).lambda    = propSpeed/txSat(2).freq;
% Receiver Information
rx.pos_lla      = [39.329956, -76.620485, 1552];
rx.pos_ecef     = lla2ecef(rx.pos_lla);
% Antenna Parameters
% Antenna is a rectangular array
rx.numRow       = 25;
rx.numCol       = 25;
% Set up antenna parameters
rx.yaw          =  135*pi/180; % positive = counterclockwise from east
rx.pitch        = -20*pi/180;  % positive = down from e-n plane
% X-axis is normal to array face, z-axis is "up" on array face
rx.dcm_enu2ant  = euler_to_matrix(rx.yaw,rx.pitch,0)';
% Setup Refernce antenna
rx_ref.pos_lla  = [39.329956, -76.620485, 1552];
rx_ref.pos_ecef = lla2ecef(rx_ref.pos_lla);
rx_ref.gain     = 20; %dBi
Loss            = 6; %loss budget
% Define Target information
tgt(1).pos_lla     = [39.366666, -76.743338, 15000];
tgt(1).pos_ecef    = lla2ecef(tgt(1).pos_lla);
tgt(1).vel_ecef     = [0,50 ,0];
% Target Radar Cross Section
tgt(1).rcs_dBsm    = 25;

%% Calculate direct-path antenna U/V coordinates
% Get transform from ecef to ENU at antenna position
numTargets = 1;
for k0 = 1:numTargets
    targetENU                              = zeros(3,1);
    [targetENU(1), targetENU(2), targetENU(3)] = ecef2enu(tgt(k0).pos_ecef(1), ...
        tgt(k0).pos_ecef(2),tgt(k0).pos_ecef(3),rx.pos_lla(1),rx.pos_lla(2), ...
        rx.pos_lla(3),wgs84Ellipsoid());
    
    % Get direct path in antenna coordinates
    targetAnt     = rx.dcm_enu2ant*targetENU;
    uTargetAnt    = targetAnt / norm(targetAnt);
    % Convert to U/V coordinates
    uTarget(k0)   = uTargetAnt(2);
    vTarget(k0)   = uTargetAnt(3);
end

%% Convert the RX array position to ENU
rx_enu      = zeros(3,1);
[rx_enu(1), rx_enu(2), rx_enu(3)] = ecef2enu(rx.pos_ecef(1),rx.pos_ecef(2), ...
    rx.pos_ecef(3),rx.pos_lla(1),rx.pos_lla(2),rx.pos_lla(3),wgs84Ellipsoid());

%% Calculate Free-space path loss
baselineRange         = norm(txSat(1).pos_ecef - rx.pos_ecef);
path_loss_baseline_dB = 20*log10(4*pi*baselineRange/txSat(1).lambda);
directPathPower_dBm   = txSat(1).erp_dBm - path_loss_baseline_dB;

%% Calculate direct-path antenna U/V coordinates
% Get transform from ecef to ENU at antenna position
directPathENU      = zeros(3,1);
[directPathENU(1), directPathENU(2), directPathENU(3)] = ecef2enu(txSat(1).pos_ecef(1),txSat(1).pos_ecef(2), ...
    txSat(1).pos_ecef(3),rx.pos_lla(1),rx.pos_lla(2),rx.pos_lla(3),wgs84Ellipsoid());

%% Get direct path in antenna coordinates
directPathAnt      = rx.dcm_enu2ant*directPathENU;
uDirectPathAnt     = directPathAnt / norm(directPathAnt);
% Convert to U/V coordinates
uDirectPath        = uDirectPathAnt(2);
vDirectPath        = uDirectPathAnt(3);


%% Signal and Noise Power
%Target signal power
for k1 = 1:numTargets
    tx2tgtRange(k1) = norm(txSat(k1).pos_ecef - tgt(k1).pos_ecef);
    rx2tgtRange(k1) = norm(rx.pos_ecef - tgt(k1).pos_ecef);
    rxGain_dB       = 32; %Rx beamformed gain
    %recieved Signal power
    tgtPower_dBm(k1)    = txSat(k1).erp_dBm + 20*log10(txSat(k1).lambda^2) + ...
        tgt(k1).rcs_dBsm - 30*log10(4*pi) - 20*log10(tx2tgtRange(k1)) - ...
        20*log10(rx2tgtRange(k1)) + rxGain_dB - Loss;
end

% Receiver noise power
k               = 1.38e-23; %Boltzman constant (J/K)
TsysTemp        = 300; %77; %temperature (Kelvin)
chanBW          = 6.00e6; %approximate signal bandwidth
noiseFigure_dB  = 0.4; 
noisePower_dBm  = 10*log10(k*TsysTemp*chanBW)+noiseFigure_dB+30; %+30 to go to dBm.

%% Configure input signal parameters
samplingFreq    = 12e6; % sampling frequency (Hz)
fcSignal        = txSat(1).freq; %2e6; % signal frequency (Hz)
numSamps        = 2^10; % number of samples
t = (0:numSamps-1)/samplingFreq; % time (s)
% phase offset between echo and direct signal (rad)
phaseOffset     = 0; 
% Time delay in micro seconds of the reflected path (tx-tgt-rx)
%timeDelay       = (tx2tgtRange + rx2tgtRange)/propSpeed; 
Ddiff           = tx2tgtRange + rx2tgtRange - baselineRange;
timeDelay       = (tx2tgtRange + rx2tgtRange - baselineRange)/propSpeed;
%samp_offset    = 40000; %round( timeDelay * samplingFreq ); % TDOA [samples] (integer) %
samp_offset     = round( timeDelay * samplingFreq );
sampsPerCycle   = 10;
% modulation rate must be less than 10% RF freq to be narrowband
cyclesPerSymbol = 16; % was 16
%% Calculate doppler shift
tgtVel = tgt(1).vel_ecef;
tgt2tx = tgt(1).pos_ecef-txSat(1).pos_ecef  ;
tgt2rx = tgt(1).pos_ecef-rx_ref.pos_ecef  ;
tx2rx  = txSat(1).pos_ecef - rx_ref.pos_ecef; %known parameter
dopplerTx = dot(tgtVel,tgt2tx)/norm(tgt2tx);
dopplerRx = dot(tgtVel,tgt2rx)/norm(tgt2rx);
FShift = dopplerTx + dopplerRx;
FShift = FShift*txSat(1).freq  /propSpeed;

theta1 = acos(dot(tx2rx,tgt2rx)/(norm(tx2rx)*norm(tgt2rx)));% the angle used in parameter extraction