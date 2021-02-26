% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

addpath('functions');
addpath('classes');

clear variables
rng('default');
% rng('shuffle');
warning ('on','all');  % Turn on when testing

%% %%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%% %%

% System Parameters
nSubcarriersPerBlock = 12;  % Number of SCs / PRB
nSymbolsPerBlock = 14;  % Number of symbols per frame
fSubcarrierSpacing = 15e3;  % SC spacing in Hz
nFft = 256;  % FFT size
nCyclicPrefix = floor(nFft * 0.067);  % CP length in samples (6.7% is the calibration assumption defined in 3GPP R1-165989)
samplingRate = nFft * fSubcarrierSpacing;  % Sampling rate

modOrder = 16; % 4 for QPSK, 16 for 16QAM, 64 for 64QAM, 256 for 256QAM
nTxAntennas = 1;
nRxAntennas = 1;

nSubcarriers = nFft;
nBitsPerBlock = nSymbolsPerBlock * nSubcarriers * log2(modOrder);

% Parameter for LTE Channel Coder
turboDecodingNumIterations = 8;
codeRate = 0.5;  % Coding rate of source data bits over transmition redundant data bits
codedTransportBlockSize = nBitsPerBlock * codeRate;
crcBitLength = 24;  % defined in 3GPP TS36.2xx
maxCodeBlockSize = 6144;  % defined in 3GPP TS36.2xx
numCodeBlockSegments = ceil((codedTransportBlockSize - crcBitLength) / maxCodeBlockSize);
if numCodeBlockSegments > 1
    tbs = codedTransportBlockSize-crcBitLength*(numCodeBlockSegments+1);  % CRC24A + CRC24B*numSegments
else
    tbs = codedTransportBlockSize-crcBitLength;  % CRC24A only
end

% Parameters for Scrambler
rnti = hex2dec('003D');
cellId = 0;
frameNo  = 0; % radio frame
codeword = 0;

% Parameter for OFDM Modulator
nGuardBands = [(nFft - nSubcarriers)/2; (nFft - nSubcarriers)/2];  % Number of guard band subcarriers in both sides

% Parameter for ChannelEstimator
activeSubcarrierIndices = 1+nGuardBands(1):nFft-nGuardBands(2);

% Parameters for Fading Channel
velocity_kmph = 500;                   % user velocity [km/h]
delaySpread_ns = 300;                  % rms delay spred for TDL channel model only, in nano second
fCenter = 0.8e9;                       % center carrier frequency [Hz]
velocity = velocity_kmph*1000/(60*60); % convert km/h -> m/s
waveLength = physconst('lightspeed')/fCenter;
fDoppler = velocity / waveLength;
txCorrMat = 1;
rxCorrMat = 1;


%% %%%%%%%%%%%%%%%%%%%%%%% ALGORITHM SELECTION %%%%%%%%%%%%%%%%%%%%%%% %%

% ----- Path Estimator -----
channelEstimationAlgorithm = 'Pilot-based iterative path estimation';
% channelEstimationAlgorithm = 'PN-based iterative estimation';
% channelEstimationAlgorithm = 'Ideal';
fprintf("Channel estimation: %s\n", channelEstimationAlgorithm);

% ----- Equalizer -----
% equalizerAlgorithm = 'Vectorized equalizer';
% eqAlgorithm = 'blockedMMSE';
equalizerAlgorithm = 'Deconvolutional equalizer';
eqAlgorithm = 'Wiener';
fprintf("Channel equalization: %s\n", equalizerAlgorithm);
    
% ----- Channel -----
channelModel = 'EVA';
fprintf("Channel profile: %s\n", channelModel);

%---- Channel Coding ----
channelCoding = 'LTE';
fprintf("Channel coding: %s\n", channelCoding);

%% %%%%%%%%%%%%%%%%%%%%%%% OBJECT GENERATION %%%%%%%%%%%%%%%%%%%%%%% %%

switch channelCoding
    case 'None'
        % Bit Generator
        hBitGenerator = RandomBitGenerator('DataLength', nBitsPerBlock);
        % Symbol Mapper
        hSymbolMapper = QamMapperV2( ...
            'ModOrder', modOrder, ...
            'InputType', 'bit', ...
            'OutputType', 'bit', ...  % 'bit' for uncoded, 'llr' for coded
            'UnitAveragePower', true, ...
            'LLROverflowPrevention', true);
    case 'LTE'
        % Bit Generator
        hBitGenerator = RandomBitGenerator('DataLength', tbs);
        % Symbol Mapper
        hSymbolMapper = QamMapperV2( ...
            'ModOrder', modOrder, ...
            'InputType', 'bit', ...
            'OutputType', 'llr', ...  % 'bit' for uncoded, 'llr' for coded
            'UnitAveragePower', true, ...
            'LLROverflowPrevention', true);
        % Symbol Mapper for uncoded
        hSymbolMapperUncoded = QamMapperV2( ...
            'ModOrder', modOrder, ...
            'InputType', 'bit', ...
            'OutputType', 'bit', ...  % 'bit' for uncoded, 'llr' for coded
            'UnitAveragePower', true, ...
            'LLROverflowPrevention', true);
        % LTE Channel Coder
        hChannelCoder = LteChannelCoder( ...
            'TurboDecodingNumIterations', turboDecodingNumIterations, ...
            'ModOrder', modOrder, ...
            'NumLayers', 1, ...
            'OutputLength', nBitsPerBlock, ...
            'LinkDirection', 'Downlink', ...
            'RedundancyVersion', 0, ...
            'TBS', tbs);

        % Scrambler
        hScrambler = LteScrambler(rnti, cellId, frameNo, codeword);
end

% OFDM Modulator
hOfdmModulator = CpOfdmModulator(...
    'FFTLength',            nFft, ...
    'NumGuardBandCarriers', nGuardBands, ...
    'NumSymbols',           nSymbolsPerBlock, ...
    'CyclicPrefixLength',   nCyclicPrefix, ...
    'InsertDCNull',         false, ...
    'PilotInputPort',       false, ...
    'PilotOutputPort',      false, ...
    'NumTransmitAntennas',  nTxAntennas, ...
    'NumReceiveAntennas',   nRxAntennas, ...
    'SubcarrierIndexOrder', 'ZeroToFFTSize');

% OTFS Precoder
hPrecoder = OtfsPrecoder( ...
    'NumDelayBins', nSubcarriers, ...
    'NumDopplerBins', nSymbolsPerBlock, ...
    'NumTransmitAntennas', nTxAntennas, ...
    'NumReceiveAntennas', nRxAntennas);

% Path Estimator
switch channelEstimationAlgorithm
    case 'Pilot-based iterative path estimation'
        hEstimator = OtfsPilotResponseBasedPathParameterEstimator( ...
            'DividingNumber', 10, ...
            'CyclicPrefixLength', nCyclicPrefix, ...  % Ncp
            'NumDopplerBins', nSymbolsPerBlock, ...   % N
            'NumDelayBins', nSubcarriers, ...         % M
            'SamplingRate', samplingRate ...
            );
        threshAlpha = 1/50;
        threshBeta = 1/10;
        localUpdateThreshold = 0.01;
    case 'PN-based iterative estimation'
%         threshold = 1/150; % This is the best performance regardless the computation time
        threshold = 1/25;  % This is the same number of paths to be estimated with the proposed CE
        hEstimator = OtfsPNSeqBasedPathParameterEstimator( ...
            'DividingNumber', 10, ...
            'CyclicPrefixLength', nCyclicPrefix, ...  % Ncp
            'NumDopplerBins', nSymbolsPerBlock, ...   % N
            'NumDelayBins', nSubcarriers, ...         % M
            'SamplingRate', samplingRate, ...
            'Threshold', threshold);
    case 'Ideal'
end

% OTFS Equalizer
switch equalizerAlgorithm
    case 'Deconvolutional equalizer'
        hEqualizer = OtfsDeconvolutionalEqualizer( ...
            'NumSymbols', nSymbolsPerBlock, ...       % N
            'CyclicPrefixLength', nCyclicPrefix, ...  % Ncp
            'NumSubcarriers', nSubcarriers, ...       % M
            'SamplingRate', samplingRate, ...
            'DividingNumber', 10, ... 
            'EqualizationAlgorithm', eqAlgorithm);   
    case 'Vectorized equalizer'
        hEqualizer = OtfsVectorizedEqualizer( ...
            'EqualizationAlgorithm', 'ZF', ...
            'NumSymbols', nSymbolsPerBlock, ...       % N
            'CyclicPrefixLength', nCyclicPrefix, ...  % Ncp
            'NumSubcarriers', nSubcarriers, ...       % M
            'SamplingRate', samplingRate, ...
            'OutputConditionNumber', false, ...
            'EqualizationAlgorithm', eqAlgorithm);   
end

% Fading Channel
hFading = SoSBasedChannel( ...  % based on in-house implementation
    'ChannelModel', channelModel, ...
    'RMSDelaySpread', delaySpread_ns, ...
    'SamplingRate', samplingRate, ...
    'DopplerFrequency', fDoppler, ...
    'NumTxAntennas', nTxAntennas, ...
    'NumRxAntennas', nRxAntennas, ...
    'TxCorrMatrix', txCorrMat, ...
    'RxCorrMatrix', rxCorrMat, ...
    'OutputCIR', true, ...
    'CyclicPrefix', nCyclicPrefix, ...
    'FFTSize', nFft, ...
    'ImpulseSource', 'Input', ...
    'SequentialOrRandomChannel', 'Sequential', ...
    'FDFMethod', 'ApplyFDFToPathParameters');

% AWGN
hAwgn = AwgnChannel('N0');


%% %%%%%%%%%%%%%%%%%%%%%%% SIMULATION SETTINGS %%%%%%%%%%%%%%%%%%%%%%% %%
% Create an empty storage to store simulation status
snrRange = 10:2:30; % in dB
simIterations = 20000 * ones(1,length(snrRange));

simstatus = table;  % use table-type variable since it's good to see in the workspace
simstatus.SNR = snrRange';
simstatus.BER = nan * zeros(length(snrRange),1);
simstatus.UncodedBER = nan * zeros(length(snrRange),1);
simstatus.BLER = nan * zeros(length(snrRange),1);
simstatus.TotalBitErrors = zeros(length(snrRange),1); 
simstatus.TotalUncodedBitErrors = zeros(length(snrRange),1); 
simstatus.TotalBlockErrors = zeros(length(snrRange),1); 
simstatus.Iteration = zeros(length(snrRange),1);
simstatus.PathCounts = zeros(length(snrRange),1);
simstatus.ComputeTime = zeros(length(snrRange),1);

time = strrep(strrep(strrep(string(datetime),'/',''),' ',''),':','');
filename = strcat('Resume-', mfilename, '-', time);

%% Load data when finding files that are created when simulation was interrupted
resumefiles = ls(sprintf('Resume-%s-*.mat', mfilename));
if not(isempty(resumefiles))
    fprintf("Found file(s) to resume a simulation:\n")
    disp(resumefiles)
    for file = resumefiles'
        prompt = strcat('Do you want to resume the simulation using "', file', '"? [Y/n/F]:');
        answer = input(prompt, 's');
        if answer == 'Y'
            disp("Now Loading...")
            load(file, 'simstatus')  % Load only the simulation status
            snrRange = simstatus.SNR';
            filename = file;
            break;
        elseif answer == 'F'
            disp("Now Loading...")
            load(file)  % Load all parameters
            snrRange = simstatus.SNR';
            filename = file;
            break;
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%% START SIMULATION %%%%%%%%%%%%%%%%%%%%%%% %%
% For ideal channel estimation
delayDopplerImpulse = zeros(nSubcarriers, nSymbolsPerBlock); delayDopplerImpulse(1,1)=sqrt(nSymbolsPerBlock*nSubcarriers);
ddImpulseInTFDomain = hPrecoder.encode(delayDopplerImpulse);
ddImpulseInTimeDomain = hOfdmModulator.modulate(ddImpulseInTFDomain);

% For PN-sequence based channel estimation
switch channelEstimationAlgorithm
    case 'PN-based estimation'
        txPNSeq = zeros((nCyclicPrefix+nFft)*nSymbolsPerBlock*10,1);  % for long PN sequence
%         txPNSeq = zeros((nCyclicPrefix+nFft)*nSymbolsPerBlock*1,1);  % for short PN sequence
        nPNSeq = 1023;
        tmpSeq = repmat(genltegoldseq(nPNSeq, de2bi(31,31)), 1, ceil(length(txPNSeq)/nPNSeq));  % 1023-length with an initial value of 31
        tmpSeq(tmpSeq==0) = -1;  % map 0 or 1 bit into -1 or 1 bit
        txPNSeq = tmpSeq(1:length(txPNSeq))';
    case 'PN-based iterative estimation'
        txPNSeq = zeros((nCyclicPrefix+nFft)*nSymbolsPerBlock,1);
        nPNSeq = 1023;
        tmpSeq = repmat(genltegoldseq(nPNSeq, de2bi(31,31)), 1, ceil(length(txPNSeq)/nPNSeq));  % 1023-length with an initial value of 31
        tmpSeq(tmpSeq==0) = -1;  % map 0 or 1 bit into -1 or 1 bit
        txPNSeq = tmpSeq(1:length(txPNSeq))';
end

charCount = 49+17;
for snrdb = snrRange
    rng('default');
    rng(106);

    fprintf('\nSNR = %2d dB \n', snrdb);
    snrIndex = find(snrRange==snrdb);
    snr = 10^(snrdb/10);
    noiseVar = 1/snr;  % Total power
    N0 = noiseVar/nFft;  % N0 is the power spectral density of noise per unit of bandwidth, which is band-limited.

    totalBitErrors = 0;
    fprintf(repmat(' ', 1, charCount));  % Print spaces in advance to avoid deleting the previously displayed characters

    tic
    for count = simstatus.Iteration(snrIndex)+1:simIterations(snrIndex)
        %% Transmitter
        % Bit Generation
        txBits = hBitGenerator.generate();
%         txBits = zeros(size(txBits));  % FOR TEST
        % Channel Encoding
        txCodedBits = hChannelCoder.encode(txBits);
        % Scrambler
        txScrampledBits = hScrambler.scramble(txCodedBits);
        % Symbol Mapping (QAM modulation)
        txSymbols = hSymbolMapper.map(txScrampledBits);

        % OTFS Modulation
        txBlocks = reshape(txSymbols, nSubcarriers, nSymbolsPerBlock);
        txPrecodedBlocks = hPrecoder.encode(txBlocks);
        % OFDM Modulation
        txSignals = hOfdmModulator.modulate(txPrecodedBlocks);
        
        %% Channel

        % Fading Channel
        switch channelEstimationAlgorithm
            case {'Pilot-based iterative path estimation', 'Ideal'} 
                hFading.initRayleighFading();
                [distortedSignals, cir] = hFading.apply(txSignals, ddImpulseInTimeDomain);
            case 'PN-based iterative estimation'
                [distortedSignals, distortedPNSeq] = hFading.apply(txSignals, txPNSeq);
        end

        % AWGN Channel
        [noisySignals, noise] = hAwgn.add(distortedSignals, N0);
        
        %% Receiver
        % OFDM Demodulator
        rxPrecodedBlocks = hOfdmModulator.demodulate(noisySignals(1:(nFft+nCyclicPrefix)*nSymbolsPerBlock));
        % OTFS Demodulator
        rxBlocks = hPrecoder.decode(rxPrecodedBlocks);
        % Equalization (Ideal)
        switch channelEstimationAlgorithm
            case 'Pilot-based iterative path estimation'
                noisyCir = hAwgn.add(cir, N0);
                chanEst = hPrecoder.decode(hOfdmModulator.demodulate(noisyCir(1:(nFft+nCyclicPrefix)*nSymbolsPerBlock)));
                [estGains, estDopplers, estDelays, estOffsets] = hEstimator.estimate(chanEst, threshAlpha, threshBeta, noiseVar);
                rxEqBlocks = hEqualizer.equalize(rxBlocks, estGains, estDopplers, estDelays, estOffsets, noiseVar, localUpdateThreshold);
            case 'PN-based iterative estimation'
                noisyPNSeq = hAwgn.add(distortedPNSeq, noiseVar*sqrt(nSubcarriers));
                [estGains, estDopplers, estDelays, estOffsets] = hEstimator.estimateIteratively(noisyPNSeq, txPNSeq, fDoppler*2);
                rxEqBlocks = hEqualizer.equalize(rxBlocks, estGains, estDopplers, estDelays, estOffsets, noiseVar, 0.01);
            case 'Ideal'
                chanEst = hPrecoder.decode(hOfdmModulator.demodulate(cir(1:(nFft+nCyclicPrefix)*nSymbolsPerBlock)));
                idealGains = hFading.pathRayleighGains;
                idealDopplers = hFading.pathDopplers;
                idealDelays = hFading.pathDelays;
                idealOffsets = hFading.pathOffsets;
                rxEqBlocks = hEqualizer.equalize(rxBlocks, idealGains, idealDopplers, idealDelays, idealOffsets, noiseVar);
        end
%         figure(1); clf; plot(rxEqBlocks, '.'); grid on; grid minor; xlim([-1.5 1.5]); ylim([-1.5 1.5]); % hold on; plot(rxIdealEqBlocks, '.'); hold off;
        rxSymbols = rxEqBlocks(:);
        % Symbol Demapping
        rxSoftBits = hSymbolMapper.demap(rxSymbols, noiseVar);
        rxUncodedHardBits = hSymbolMapperUncoded.demap(rxSymbols);
        % Descrambling
        rxDescrambledBits = hScrambler.descramble(rxSoftBits);
        % Channel Decoding
        [rxHardBits, blockError] = hChannelCoder.decode(rxDescrambledBits);
        
        %% Result
        % Bit Error Rate
        numBitErrors = sum(xor(txBits, rxHardBits));
        simstatus.TotalBitErrors(snrIndex) = simstatus.TotalBitErrors(snrIndex) + numBitErrors;
        simstatus.BER(snrIndex) = simstatus.TotalBitErrors(snrIndex)/(count*nBitsPerBlock);
        simstatus.TotalBlockErrors(snrIndex) = simstatus.TotalBlockErrors(snrIndex) + blockError;
        simstatus.BLER(snrIndex) = simstatus.TotalBlockErrors(snrIndex)/count;
        numBitErrors = sum(xor(txScrampledBits, rxUncodedHardBits));
        simstatus.TotalUncodedBitErrors(snrIndex) = simstatus.TotalUncodedBitErrors(snrIndex) + numBitErrors;
        simstatus.UncodedBER(snrIndex) = simstatus.TotalUncodedBitErrors(snrIndex)/(count*length(txScrampledBits));
        fprintf(repmat('\b',1, charCount));  % Delete the prevously displayed characters, and then update the display
        fprintf('[%5d /%5d] BLER: %1.4f, BER: %1.7f (%9d/ %9d)', count, simIterations(snrIndex), simstatus.BLER(snrIndex), simstatus.UncodedBER(snrIndex), simstatus.TotalUncodedBitErrors(snrIndex), count*length(txScrampledBits));
        
        if mod(count,100)==0
            simstatus.Iteration(snrIndex) = count;
            save(filename);  % save all parameters to resume
        end
        pause(0.01)
    end
    computationTime = toc;
    computationTimePerIteration = computationTime/simIterations(snrIndex);
    fprintf('  Computation Time : %f sec. (as per iteration)\n', computationTimePerIteration);
    simstatus.ComputeTime(snrIndex) = computationTimePerIteration;
end

%% Save results
save(filename)

%% Show results
if usejava('jvm')  % if GUI is available
    figure
    semilogy(simstatus.SNR, simstatus.BER);
    ylabel('BER')
    xlabel('SNR (dB)')
    
    figure
    semilogy(simstatus.SNR, simstatus.BLER);
    ylabel('BLER')
    xlabel('SNR (dB)')     
end
