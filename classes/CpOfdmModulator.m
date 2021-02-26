% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef CpOfdmModulator < matlab.System
    
    properties
        FFTLength;
        NumGuardBandCarriers;
        NumSymbols;
        CyclicPrefixLength;
        InsertDCNull;
        PilotInputPort;
        NumTransmitAntennas;
        PilotOutputPort;
        NumReceiveAntennas;
        SubcarrierIndexOrder;
        
        OfdmModObj;
        OfdmDemodObj;
        
        OfdmModMat;
        OfdmDemMat;
        CpAddMat;
        CpRemMat;
    end
    
    methods
        function this = CpOfdmModulator(varargin) % Constructor
            if nargin > 0
                setProperties(this, nargin, varargin{:}, ...
                    'FFTLength', ...
                    'NumGuardBandCarriers', ...
                    'NumSymbols', ...
                    'CyclicPrefixLength', ...
                    'InsertDCNull', ...
                    'PilotInputPort', ...
                    'NumTransmitAntennas', ...
                    'PilotOutputPort', ...
                    'NumReceiveAntennas', ...
                    'SubcarrierIndexOrder' ...
                );

                if isempty(this.SubcarrierIndexOrder)
                    this.SubcarrierIndexOrder = 'LowToHigh';  % default setting of comm.OFDMModulator/Demotulator
                end
                
                % OFDM modulator object
                this.OfdmModObj = comm.OFDMModulator( ...
                    'FFTLength',            this.FFTLength, ...
                    'NumGuardBandCarriers', this.NumGuardBandCarriers, ...
                    'InsertDCNull',         this.InsertDCNull, ...
                    'PilotInputPort',       this.PilotInputPort, ...
                    'CyclicPrefixLength',   this.CyclicPrefixLength, ...
                    'NumSymbols',           this.NumSymbols, ...
                    'NumTransmitAntennas',  this.NumTransmitAntennas);

                % OFDM demodulator object
                this.OfdmDemodObj = comm.OFDMDemodulator( ...
                    'FFTLength',            this.FFTLength, ...
                    'NumGuardBandCarriers', this.NumGuardBandCarriers, ...
                    'PilotOutputPort',      this.PilotOutputPort, ...
                    'CyclicPrefixLength',   this.CyclicPrefixLength, ...
                    'NumSymbols',           this.NumSymbols, ...
                    'NumReceiveAntennas',   this.NumReceiveAntennas);

%                 % Same operation with the order 'ZeroToFFTSize' as OFDMModulator/Demodulator
%                 modAllocMask = this.NumGuardBandCarriers(1)+1 : this.FFTLength-this.NumGuardBandCarriers(2);
%                 OfdmModulationFullMatrix = exp(1i*2*pi*(0:this.FFTLength-1).'*(0:this.FFTLength-1)/this.FFTLength);
%                 this.OfdmModMat = 1/this.FFTLength * OfdmModulationFullMatrix(:,modAllocMask);
%                 this.OfdmDemMat = this.OfdmModMat';
%                 this.CpAddMat = [[zeros(this.CyclicPrefixLength,this.FFTLength-this.CyclicPrefixLength) eye(this.CyclicPrefixLength)]; eye(this.FFTLength)];
%                 this.CpRemMat = [zeros(this.FFTLength, this.CyclicPrefixLength) eye(this.FFTLength)];

            end
        end

        function txSignal = modulate(this, txGridFreqTime)
            switch this.SubcarrierIndexOrder
                case 'LowToHigh'
                    txSignal = step(this.OfdmModObj, txGridFreqTime);
                case 'ZeroToFFTSize'
                    txSignal = step(this.OfdmModObj, fftshift(txGridFreqTime,1));
%                     txSignal = this.CpAddMat * this.OfdmModMat * txGridFreqTime;
%                     txSignal = txSignal(:);
            end
        end
        
        function rxGridFreqTime = demodulate(this, rxSignal)
            switch this.SubcarrierIndexOrder
                case 'LowToHigh'
                    rxGridFreqTime = step(this.OfdmDemodObj, rxSignal);
                case 'ZeroToFFTSize'
                    rxGridFreqTime = step(this.OfdmDemodObj, rxSignal);
                    rxGridFreqTime = ifftshift(rxGridFreqTime,1);
%                     rxGridFreqTime = this.FFTLength * this.OfdmDemMat * this.CpRemMat * reshape(rxSignal, this.CyclicPrefixLength+this.FFTLength, this.NumSymbols);
            end
        end
    end
end