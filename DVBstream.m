% According to https://www.mathworks.com/help/comm/examples/dvb-s-2-link-including-ldpc-coding.html
% The DVB signal:symbol rate = 2.2e7/s  fc = 10714e6 Hz 
% Typically 1 frame has 33900 symbols, which last about 33900/110k =
% 300ms
% This function construct a given frame number stream of using DVB modulation. If the stream is longer than 
% the maxSymbol, then use the maxSymbol as the final length.
% e.g  stream = DVBstream(1,1,40000) or  stream = DVBstream(1,1,30000)
% This function output a N*1 array.
function stream = DVBstream(numFrames,pilotOn,maxSymbol)

%% Initialization
% The <matlab:edit('configureDVBS2Demo.m') configureDVBS2Demo.m> script
% initializes some simulation parameters and generates a structure, dvb.
% The fields of this structure are the
%parameters of the DVB-S.2 system at
% hand. It also creates the System objects making up the DVB-S.2 system.

subsystemType = 'QPSK 5/6';   % Constellation and LDPC code rate
EsNodB        = 9;              % Energy per symbol to noise PSD ratio in dB
%numFrames     = 20;             % Number of frames to simulate
%pilotOn       = 1;              %% Define whether add pilot block into it

% Initialize
configureDVBS2Demo

% Display system parameters
dvb;

%%
% The following is a list of objects this example uses:
%
% *Simulation objects:*
%
%  enc              - BCH encoder
%  dec              - BCH decoder
%  LDPCEnc          - LDPC encoder
%  LDPCDec          - LDPC decoder
%  intrlvr          - Block interleaver
%  deintrlvr        - Block deinterleaver
%  pskModulator     - PSK modulator
%  pskDemodulator   - PSK demodulator
%  chan             - AWGN channel
%
% *Performance measurement objects:*
%
%  PER          - Packet error rate calculator
%  BERLDPC      - LDPC decoder output error rate calculator
%  BERMod       - Demodulator output error rate calculator
%  constDiag    - Scatter plot of channel output
%  varCalc      - Variance of the noise on a frame
%  meanCalc     - Average of the noise variance
%
% The following is a list of functions this example uses:
%
% *Simulation functions:*
%
%  dvbsapskmod      - DVBSAPSK modulator
%  dvbsapskdemod    - DVBSAPSK demodulator

%% LDPC Encoder and Decoder
% Create LDPC encoder and decoder System objects and set the parity check
% matrix according to Section 5.3.1 of the DVB-S.2 standard [ <#12 1> ]. You set
% the IterationTerminationCondition property to 'Parity check satisfied' to
% stop the decoder iterations when all the parity checks are satisfied,
% which reduces the decoding time. Set the MaximumIterationCount property to
% 50, to limit the number of simulation iterations. Set the
% NumIterationsOutputPort to true to output the number of iterations
% performed for each codeword.

encldpc = comm.LDPCEncoder(dvb.LDPCParityCheckMatrix);

decldpc = comm.LDPCDecoder(dvb.LDPCParityCheckMatrix, ...
    'IterationTerminationCondition', 'Parity check satisfied', ...
    'MaximumIterationCount',         dvb.LDPCNumIterations, ...
    'NumIterationsOutputPort',       true);

%% Stream Processing Loop
% This section of the code calls the processing loop for a DVB-S.2 system.
% The main loop processes the data frame-by-frame, where the system
% parameter dvb.NumPacketsPerBBFrame determines the number of data packets
% per BB frame. The first part of the for-loop simulates the system. The
% simulator encodes each frame using BCH and LDPC encoders as inner and
% outer codes, respectively. The encoded bits pass through an interleaver.
% The modulator maps the interleaved bits to symbols from the predefined
% constellation. The modulated symbols pass through an AWGN channel. The
% demodulator employs an approximate log-likelihood algorithm to obtain soft
% bit estimates. The LDPC decoder decodes the deinterleaved soft bit values and
% generates hard decisions. The BCH decoder works on these hard decisions to
% create the final estimate of the received frame.
% 
% The second part of the for-loop collects performance measurements such as the
% bit error rate and a scatter plot. It also estimates the received SNR value.

bbFrameTx  = false(encbch.MessageLength,1);
numIterVec = zeros(numFrames, 1);
falseVec   = false(dvb.NumPacketsPerBBFrame, 1);

finalOut = [];
for frameCnt=1:numFrames
    
    % Transmitter, channel, and receiver
    bbFrameTx(1:dvb.NumInfoBitsPerCodeword) = ...
          logical(randi([0 1], dvb.NumInfoBitsPerCodeword, 1));
    
    bchEncOut = encbch(bbFrameTx);    
    ldpcEncOut = encldpc(bchEncOut);
    %intrlvrOut = intrlvr(ldpcEncOut); %2019a
    intrlvrOut = intrlv(ldpcEncOut,dvb.InterleaveOrder); %2019b
    
    if dvb.ModulationOrder == 4 || dvb.ModulationOrder == 8
        modOut = pskModulator(intrlvrOut);
    else
        modOut = dvbsapskmod(intrlvrOut, dvb.ModulationOrder, 's2', ...
            dvb.CodeRate, 'InputType', 'bit', 'UnitAveragePower', true);
    end
    % Define pilot block
    pilot = sqrt(2)*(0.5+0.5*1i)*ones(64,1);
    SOF = hexToBinaryVector('18D2E82');
    % This part is not the same as Pg31. The latter part should not be random but calculated by some formula. To be revised later  
    PLHeader = cat(1,SOF',(randi([0,1],65,1)));
    FrameOut = [];
    
    if pilotOn == 0
        FrameOut = cat(PLHeader,modOut);
    else
        for i = 0:16:(length(modOut)/90)-16
            FrameOut = cat(1,FrameOut,modOut(90*i+1:90*(i+16),1),pilot);
        end
        FrameOut = cat(1,FrameOut,modOut(90*(i+16)+1:end)); % add tail
        FrameOut = cat(1,PLHeader,FrameOut);                % add head
    end
end
    
    finalOut = cat(1,finalOut,FrameOut);

if maxSymbol < length(finalOut)
    stream = zeros(maxSymbol,1);
    stream = finalOut(1:maxSymbol);
else
    stream = finalOut;
end
end


