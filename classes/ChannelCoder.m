% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef ChannelCoder < matlab.System
    
    properties
        % Input properties
        TurboDecodingNumIterations;
        TurboCoderInterleaverRandomSeed;
        OutputLength;
        TBS;        
        ChannelInterleaverRandomSeed;
        RateMatchingInterleaverRandomSeed;

        % Internal properties
        turboCodedBitLength;

        % Handlers
        hCrc24Coder;
        hTurboCoder;
        hRateMatchInterleaver;
        hChannelInterleaver;
        
        bitSeqA;
        bitSeqB;
        bitSeqC;
        bitSeqD;
        bitSeqE;
        bitSeqF;

    end
    
    methods
        function this = ChannelCoder(varargin)  % Constructor
            setProperties(this, nargin, varargin{:}, ...
                'TurboDecodingNumIterations', ...
                'TurboCoderInterleaverRandomSeed', ...
                'OutputLength', ...
                'TBS', ...
                'ChannelInterleaverRandomSeed', ...
                'RateMatchingInterleaverRandomSeed');

            % CRC 24 bit
            this.hCrc24Coder = Crc24Coder('A');
            crcLength = 24;
            % Turbo coder
            interleaverIndices = randintrlv(1:this.TBS+crcLength, this.TurboCoderInterleaverRandomSeed);
            this.hTurboCoder = TurboCoder( ...
                'InterleaverIndicesSource', 'Property', ...
                'InterleaverIndices', interleaverIndices, ...
                'NumIterations', this.TurboDecodingNumIterations);
            tailBitLength = log2(this.hTurboCoder.trellis.numStates) * log2(this.hTurboCoder.trellis.numOutputSymbols);
            this.turboCodedBitLength = 3*(this.TBS+crcLength) + 2*tailBitLength;            
            % Rate matching interleaver
            this.hRateMatchInterleaver = Interleaver( ...
                'InterleaverIndicesSource', 'Random', ...
                'InputLength', this.turboCodedBitLength, ...
                'InterleaverRandomSeed', this.RateMatchingInterleaverRandomSeed);
            % Channel interleaver
            this.hChannelInterleaver = Interleaver( ...
                'InterleaverIndicesSource', 'Random', ...
                'InputLength', this.OutputLength, ...
                'InterleaverRandomSeed', this.ChannelInterleaverRandomSeed);
           
        end
        
        function bitSeqF = encode(this, in)
            bitSeqA = in;
            bitSeqB = this.hCrc24Coder.generate(bitSeqA);
            bitSeqC = this.hTurboCoder.encode(bitSeqB);
            bitSeqD = this.hRateMatchInterleaver.interleave(bitSeqC);
            bitSeqE = bitSeqD(1:this.OutputLength);
            bitSeqF = this.hChannelInterleaver.interleave(bitSeqE);
            this.bitSeqA=bitSeqA;this.bitSeqB=bitSeqB;this.bitSeqC=bitSeqC;this.bitSeqD=bitSeqD;this.bitSeqE=bitSeqE;this.bitSeqF=bitSeqF;
        end
        
        function [bitSeqA, blkErr] = decode(this, in)
            bitSeqE = this.hChannelInterleaver.deinterleave(in);
            bitSeqD = zeros(this.turboCodedBitLength,1);
            bitSeqD(1:this.OutputLength) = bitSeqE;
            bitSeqC = this.hRateMatchInterleaver.deinterleave(bitSeqD);
            bitSeqB = this.hTurboCoder.decode(-bitSeqC); 
            [bitSeqA, blkErr] = this.hCrc24Coder.detect(bitSeqB);
        end

        function [bitSeqA, error, bitSeqFEst] = decodeWithUpdatedLlrOutput(this, in)
            bitSeqE = this.hChannelInterleaver.deinterleave(in);
            bitSeqD = zeros(this.turboCodedBitLength,1);
            bitSeqD(1:this.OutputLength) = bitSeqE;
            bitSeqC = this.hRateMatchInterleaver.deinterleave(bitSeqD);
            [bitSeqB, bitSeqCEst, bitSeqA, error] = this.hTurboCoder.decode(bitSeqC, this.hCrc24Coder); % TODO: implement CRC check in the Turbo decoding iteration loop
            bitSeqDEst = this.hRateMatchInterleaver.interleave(bitSeqCEst);
            bitSeqEEst = bitSeqDEst(1:this.puncturedBitLength);
            bitSeqFEst = this.hChannelInterleaver.interleave(bitSeqEEst);
        end
        
    end
end