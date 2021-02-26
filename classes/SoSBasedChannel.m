% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef SoSBasedChannel < matlab.System
    properties
        % Parameters defined when creating an object
        ChannelModel;  % ETU, TDL-A, TDL-B, TDL-C
        RMSDelaySpread;  % Root mean square values of the delay spread for TDL channel models in nano second
        SequentialOrRandomChannel;  % This indicates whether the channel is continous in time
        NormalizePathGains;  % This indicates whether the average path gians are normalized
        ImpulseSource;  % Input source of the impulse signal
        OutputCIR;  % This indicates whether the channel impulse response should be outputted
        FDFMethod;  % This indicates whether the FDF is applied to input signals or path parameters
        
        SamplingRate;
        FFTSize;
        CyclicPrefix;
        DopplerFrequency;
        NumTxAntennas;  % To Do
        NumRxAntennas;  % To Do
        TxCorrMatrix;  % To Do
        RxCorrMatrix;  % To Do

        % Parameters defined in this class
        delayProfile;  % determined according to the channel model
        pathAvgGainsdB;  % determined according to the channel model

        pathAvgGains;  % Average gain for each path
        pathRayleighGains;  % Gain for each path after Rayleigh faiding
        pathDelays;
        pathOffsets;
        pathDopplers;
        
        pathDelaysInSample;
        
        pathOrigDopplers;
        pathOrigRayleighGains;

        initialTimeIndex;  % time index of SoS-based channel
        FDF;  % Fractional delay filter object

    end
        
    methods
        function this = SoSBasedChannel(varargin) % Constructor
            if nargin > 0
                setProperties(this, nargin, varargin{:}, ...
                    'ChannelModel', ...
                    'RMSDelaySpread', ...
                    'SamplingRate', ...
                    'DopplerFrequency', ...
                    'NumBSAntennas', ...
                    'NumUEAntennas', ...
                    'TxCorrMatrix', ...
                    'RxCorrMatrix', ...
                    'OutputCIR', ...
                    'FFTSize', ...
                    'CyclicPrefix', ...
                    'ImpulseSource', ...
                    'SequentialOrRandomChannel', ...
                    'NormalizePathGains');

                % Set channel model
                switch this.ChannelModel
                    case 'EPA' % defined in 3GPP TS 36.104 Annex B
                        this.delayProfile   = [  0   30   70   90  110   190   410] * 1e-9; % in second
                        this.pathAvgGainsdB = [0.0 -1.0 -2.0 -3.0 -8.0 -17.2 -20.8]; % in dB

                    case 'EVA' % defined in 3GPP TS 36.104 Annex B
                        this.delayProfile   = [  0   30  150  310  370  710 1090  1730  2510] * 1e-9; % in second
                        this.pathAvgGainsdB = [0.0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9]; % in dB

                    case 'ETU' % defined in 3GPP TS 36.104 Annex B
                        this.delayProfile   = [   0   50  120 200 230 500 1600 2300 5000] * 1e-9; % in second
                        this.pathAvgGainsdB = [-1.0 -1.0 -1.0 0.0 0.0 0.0 -3.0 -5.0 -7.0]; % in dB

                    case 'TDLA' % defined in 3GPP TS 38.900 v1.0.0 (2016-06), Table 7.7.2-1 
                        this.delayProfile = this.RMSDelaySpread*[0,0.3819,0.4025,0.5868,0.461,0.5375,0.6708,0.575,0.7618,1.5375,1.8978,2.2242,2.1718,2.4942,2.5119,3.0582,4.081,4.4579,4.5695,4.7966,5.0066,5.3043,9.6586]' * 1e-9; % in second
                        this.pathAvgGainsdB = [-13.4,0,-2.2,-4,-6,-8.2,-9.9,-10.5,-7.5,-15.9,-6.6,-16.7,-12.4,-15.2,-10.8,-11.3,-12.7,-16.2,-18.3,-18.9,-16.6,-19.9, -29.7]'; % in dB

                    case 'TDLB' % defined in 3GPP TS 38.900 v1.0.0 (2016-06), Table 7.7.2-2 
                        this.delayProfile = this.RMSDelaySpread*[0,0.1072,0.2155,0.2095,0.287,0.2986,0.3752,0.5055,0.3681,0.3697,0.57,0.5283,1.1021,1.2756,1.5474,1.7842,2.0169,2.8294,3.0219,3.6187,4.1067,4.279,4.7834]' * 1e-9; % in second
                        this.pathAvgGainsdB = [0,-2.2,-4,-3.2,-9.8,-1.2,-3.4,-5.2,-7.6,-3,-8.9,-9,-4.8,-5.7,-7.5,-1.9,-7.6,-12.2,-9.8,-11.4,-14.9,-9.2,-11.3]'; % in dB

                    case 'TDLC' % defined in 3GPP TS 38.900 v1.0.0 (2016-06), Table 7.7.2-3 
                        this.delayProfile = this.RMSDelaySpread* [0, 0.2099, 0.2219, 0.2329, 0.2176, 0.6366, 0.6448, 0.656, 0.6584, 0.7935, 0.8213, 0.9336, 1.2285, 1.3083, 2.1704, 2.7105, 4.2589, 4.6003, 5.4902, 5.6077, 6.3065, 6.6374, 7.0427, 8.6523]' * 1e-9 + 1/this.SamplingRate; % in second
                        this.pathAvgGainsdB = [-4.4, -1.2, -3.5, -5.2, -2.5, 0, -2.2, -3.9, -7.4, -7.1, -10.7, -11.1, -5.1, -6.8, -8.7, -13.2, -13.9, -13.9, -15.8, -17.1, -16, -15.7, -21.6, -22.8]'; % in dB        

                    case 'TEST'
                        this.delayProfile = [0]';
                        this.pathAvgGainsdB = [0]';

                    case 'TEST2'
                        this.delayProfile = [1]'/this.SamplingRate;
                        this.pathAvgGainsdB = [0]';
                end
                
                switch this.SequentialOrRandomChannel
                    case 'Sequential'
                        this.initialTimeIndex = 0.0;
                    case 'Random'
                        this.initialTimeIndex = double(randi([1 intmax]))/this.SamplingRate;
                end
            end
            
            if isempty(this.NormalizePathGains) || this.NormalizePathGains == true
                this.NormalizePathGains = true;
                this.pathAvgGains = 10.^(this.pathAvgGainsdB/10) / sum(10.^(this.pathAvgGainsdB/10));
                this.pathAvgGainsdB = 10*log10(this.pathAvgGains);
            end
            
            this.initRayleighFading();
            this.initialTimeIndex = 0;
        end

        
        %% Rayleigh Fading
        % This implementation is based on the following paper in order to support MIMO channel (multiple uncorrelated paths)
        % [2] Patzold, Matthias, Cheng-Xiang Wang, and Bjorn Olav Hogstad. “Two New Sum-of-Sinusoids-Based Methods for the Efficient Generation of Multiple Uncorrelated Rayleigh Fading Waveforms.? IEEE Transactions on Wireless Communications 8, no. 6 (June 2009): 3122?31. https://doi.org/10.1109/TWC.2009.080769.
        function initRayleighFading(this)
           
            numSinusoids = length(this.pathAvgGains);
            % Path Dopplers
            if this.NumRxAntennas == 1  % for K=1, Jakes model (Monte Carlo Method)
                aoa = 2*pi*rand(1,numSinusoids);  % angle of arrival
                this.pathDopplers = this.DopplerFrequency * cos(aoa);
            else
                % Doppler frequency of each path
%                 % --- for the GMEDS1 for K=1 ---
%                 q = 1;
%                 alphai = -pi/(4*numSinusoids) * 1/(1+3);
%                 alphaq = pi/(4*numSinusoids) * 1/(1+3);
%                 alphai = q*pi/(2*numSinusoids) * (n - 1/2) + alphai;
%                 alphaq = q*pi/(2*numSinusoids) * (n - 1/2) + alphaq;
%                 this.pathDopplers(n,1) = this.DopplerFrequency * cos(alphai);
%                 this.pathDopplers(n,2) = this.DopplerFrequency * cos(alphaq);
                % --- for the GMEDS2 for K=1 ---
                q = 2;
%                 alphaper = pi/(2*numSinusoids);
%                 alpha0 = alphaper/2*0;
%                 alpha = q*pi/(2*numSinusoids) * (n - 1/2) + alpha0;
%                 this.pathDopplers(n) = this.DopplerFrequency * cos(alpha);
            end                
            
            % Path offset
            this.pathOffsets = 2*pi*rand(1,numSinusoids);
%             this.pathOffsets = zeros(numSinusoids,1);
            
            % Rayleigh amplitude
            this.pathRayleighGains = sqrt(this.pathAvgGains) .* abs(randn(1,numSinusoids)+1i*randn(1,numSinusoids))/sqrt(2);
            

            if strcmp(this.FDFMethod, 'ApplyFDFToPathParameters')
                % The FDF is not a causal filter, especially if the delay is less than 1.
                % In this case, the FDF cannot keep the energy of an input signal such as a discontinuous signal (e.g., delta function).
                % To keep the same energy between discontinuous signals and continous signals, operations below convert parameters to express fractional delays.
                this.pathOrigRayleighGains = this.pathRayleighGains;
                this.pathOrigDopplers = this.pathDopplers;
                [this.pathRayleighGains, this.pathDopplers, this.pathDelays, this.pathOffsets] = fdfpathpars(this.pathRayleighGains, this.pathDopplers, this.delayProfile, this.pathOffsets, this.SamplingRate);
                this.pathDelaysInSample = this.pathDelays*this.SamplingRate;
            elseif strcmp(this.FDFMethod, 'Floor')
                this.pathDelaysInSample = ceil(this.delayProfile*this.SamplingRate);
                this.pathDelays = this.pathDelaysInSample/this.SamplingRate;
            elseif strcmp(this.FDFMethod, 'Round')
                this.pathDelaysInSample = round(this.delayProfile*this.SamplingRate);
                this.pathDelays = this.pathDelaysInSample/this.SamplingRate;
            else
                this.pathDelaysInSample = this.delayProfile*this.SamplingRate;
                this.pathDelays = this.delayProfile;
            end            
            this.FDF = FractionalDelayFilter;
            this.FDF.init(this.pathDelaysInSample);

            switch this.SequentialOrRandomChannel
                case 'Sequential'
                    this.initialTimeIndex = 0;
                case 'Random'
                    this.initialTimeIndex = double(randi([1 intmax]))/this.SamplingRate;
            end            
        end
        
        function [out, cir, channelGain] = apply(this, in, varargin) % Channel Impulse Response: CIR
            %% ----- Input parameters -----
            %   - varargin{1}: impulse source when `this.ImpulseSource == 'Input'`
            
            %% Apply for impulse response
            if this.OutputCIR
                switch this.ImpulseSource
                    case 'Input'
                        impulse = varargin{1};
                    case 'TimeDomainImpulse'
                        impulse = zeros(size(in)); impulse(this.CyclicPrefix+1:this.FFTSize+this.CyclicPrefix:end) = 1.0; % After cyclic prefix samples
                end
                cir = this.applySoSChannel(impulse);
                
%                 % Update time index
%                 switch this.SequentialOrRandomChannel
%                     case 'Sequential'
%                         this.initialTimeIndex = this.initialTimeIndex + length(in)*1/this.SamplingRate;
%                     case 'Random'
%                         this.initialTimeIndex = double(randi([1 intmax]))/this.SamplingRate;
%                 end
            else
                cir = nan;
            end

            %% Apply for input signal
            out = this.applySoSChannel(in);
%             out = this.applySoSChannel(impulse);
           
            % Calc channel gain
            channelGain = (rms(out)/rms(in))^2;  % in power
            
            % Update time index
            switch this.SequentialOrRandomChannel
                case 'Sequential'
                    this.initialTimeIndex = this.initialTimeIndex + length(in)*1/this.SamplingRate;
                case 'Random'
                    this.initialTimeIndex = double(randi([1 intmax]))/this.SamplingRate;
            end
            
        end

        %% Sum of sinusoid based channel model
        function out = applySoSChannel(this, in)
            % Apply the fractional delay filter
            delayedInput = this.FDF.apply(in);
            
            % Sum of sinusoid
            out = zeros(size(in));
            numSinusoids = length(this.pathRayleighGains);  % # of sinusoids
            for n=1:numSinusoids
                t = ((0:length(in)-1)-this.pathDelaysInSample(n)) / this.SamplingRate + this.initialTimeIndex;
                phaseShift = exp( 1i * (2*pi * this.pathDopplers(n) * t' + this.pathOffsets(n)) );
                
                out = out + this.pathRayleighGains(n) * phaseShift .* delayedInput(1:length(in),n);
            end
        end
        
        %% Compatibility for the old interface
        function [out, cir, channelGain] = add(this, in, varargin)
            channelGain = nan;
            [out, cir] = this.apply(in, varargin{1});
        end

    end
end
        