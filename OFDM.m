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
crcBitLength = 24;  % defined in 3GPP TS36.213
maxCodeBlockSize = 6144;  % defined in 3GPP TS36.213
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
velocity_kmph = 500;                     % user velocity [km/h]
delaySpread_ns = 300;                  % rms delay spred for TDL channel model only, in nano second
fCenter = 0.8e9;                         % center carrier frequency [Hz]
velocity = velocity_kmph*1000/(60*60); % convert km/h -> m/s
waveLength = physconst('lightspeed')/fCenter;
fDoppler = velocity / waveLength;
txCorrMat = 1;
rxCorrMat = 1;


%% %%%%%%%%%%%%%%%%%%%%%%% ALGORITHM SELECTION %%%%%%%%%%%%%%%%%%%%%%% %%
%---- Channel ----
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

% OFDM Channel Estimator
hChannelEstimator = FreqDomainChannelEstimator( ...
    'EstimationMethod', 'Ideal', ...
    'InputType', 'TimeDomainChannelImpulseResponse', ...
    'CPOFDMModulator', hOfdmModulator, ...
    'FFTSize', nFft, ...
    'CyclicPrefix', nCyclicPrefix, ...
    'NumGuardBandCarriers', nGuardBands);

% OFDM Equalizer
hEqualizer = FreqDomainEqualizer( ...
    'EqualizationAlgorithm', 'MMSE');

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
    'ImpulseSource', 'TimeDomainImpulse', ...
    'SequentialOrRandomChannel', 'Random', ...
    'FDFMethod', 'ApplyFDFToPathParameters');

% AWGN
hAwgn = AwgnChannel('N0');


%% %%%%%%%%%%%%%%%%%%%%%%% SIMULATION SETTINGS %%%%%%%%%%%%%%%%%%%%%%% %%
% Create an empty storage to store simulation status
snrRange = 10:2:30; % in dB
simIterations = 20000;

simstatus = table;  % use table-type variable since it's good to see in the workspace
simstatus.SNR = snrRange';
simstatus.BER = nan * zeros(length(snrRange),1);
simstatus.UncodedBER = nan * zeros(length(snrRange),1);
simstatus.BLER = nan * zeros(length(snrRange),1);
simstatus.TotalBitErrors = zeros(length(snrRange),1); 
simstatus.TotalUncodedBitErrors = zeros(length(snrRange),1); 
simstatus.TotalBlockErrors = zeros(length(snrRange),1); 
simstatus.Iteration = ones(length(snrRange),1);
time = strrep(strrep(strrep(string(datetime),'/',''),' ',''),':','');
filename = strcat('Resume-', mfilename, '-', time);

%% Load data when finding files that are created when simulation was interrupted
resumefiles = ls(sprintf('Resume-%s-*.mat', mfilename));
if not(isempty(resumefiles))
    fprintf("Found file(s) to resume a simulation:\n")
    disp(resumefiles)
    for file = resumefiles'
        prompt = strcat('Do you want to resume the simulation using "', file', '"? [Y/n]:');
        answer = input(prompt, 's');
        if answer == 'Y'
            disp("Now Loading...")
            load(file, 'simstatus')  % Load only the simulation status
            snrRange = simstatus.SNR';
            filename = file;
            simstatus.Iteration = simstatus.Iteration + 1;
            break;
        end
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%% START SIMULATION %%%%%%%%%%%%%%%%%%%%%%% %%
charCount = 49+17;
for snrdb = snrRange
    rng('default');
    rng(100);
    
    fprintf('\nSNR = %2d dB \n', snrdb);
    snrIndex = find(snrRange==snrdb);
    snr = 10^(snrdb/10);
    noiseVar = 1/snr;  % Total power power
    N0 = noiseVar/nFft;  % N0 is the power spectral density of noise per unit of bandwidth, which is band-limited.

    totalBitErrors = 0;
    fprintf(repmat(' ', 1, charCount));  % Print spaces in advance to avoid deleting the previously displayed characters
    aveChannelGain = zeros(simIterations,1);
    
    tic
    for count = simstatus.Iteration(snrIndex):simIterations
        %% Transmitter
        % Bit Generation
        txBits = hBitGenerator.generate();
        % Channel Encoding
        txCodedBits = hChannelCoder.encode(txBits);
        % Scrambler
        txScrampledBits = hScrambler.scramble(txCodedBits);
        % Symbol Mapping (QAM modulation)
        txSymbols = hSymbolMapper.map(txScrampledBits);
        % OFDM Modulation
        txBlocks = reshape(txSymbols, nSubcarriers, nSymbolsPerBlock);
        txSignals = hOfdmModulator.modulate(txBlocks);
        
        %% Channel
        % Fading Channel
        hFading.initRayleighFading();
        [distortedSignals, cir] = hFading.apply(txSignals);  % for SoSBasedChannel
        
        % AWGN Channel
        noisySignals = hAwgn.add(distortedSignals, N0);
        noisycir = hAwgn.add(cir, N0);
        
        %% Receiver
        % OFDM Demodulator
        rxBlocks = hOfdmModulator.demodulate(noisySignals(1:(nFft+nCyclicPrefix)*nSymbolsPerBlock));
        % Channel Estimation (Ideal, not used in this simulation)
        chanEst = hChannelEstimator.estimate(noisycir(1:(nFft+nCyclicPrefix)*nSymbolsPerBlock));
        % Equalization (Ideal)
        rxSymbols = rxBlocks(:);
        [rxEqSymbols, noiseVarEq] = hEqualizer.equalize(rxSymbols, chanEst, noiseVar);  
%         figure(1); clf; plot(rxEqSymbols, '.'); xlim([-1.5 1.5]); ylim([-1.5 1.5]);% hold on; plot(rxIdealEqBlocks, '.'); hold off;

        % Symbol Demapping
        rxSoftBits = hSymbolMapper.demap(rxEqSymbols, noiseVar);
        rxUncodedHardBits = hSymbolMapperUncoded.demap(rxEqSymbols);
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
        fprintf('[%5d /%5d] BLER: %1.4f, BER: %1.7f (%9d/ %9d)', count, simIterations, simstatus.BLER(snrIndex), simstatus.UncodedBER(snrIndex), simstatus.TotalUncodedBitErrors(snrIndex), count*length(txScrampledBits));
        
        if mod(count,10)==0
            simstatus.Iteration(snrIndex) = count;
            save(filename);  % save all parameters to resume
        end
    end
    computationTime = toc;
    computationTimePerIteration = computationTime/simIterations;
    fprintf('  Computation Time : %f sec. (as per iteration)\n', computationTimePerIteration);
end

%% Save all parameters
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
