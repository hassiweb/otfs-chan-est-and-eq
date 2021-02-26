% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef TurboCoder < matlab.System
    
    properties
        NumIterations;
        InterleaverIndices;
        InterleaverIndicesSource;

        trellis;
        InternalInterleaver;
        CoderObj;
        DecoderObj;
        AppDecObj1;
        AppDecObj2;
    end
    
    methods
        function this = TurboCoder(varargin)
            setProperties(this, nargin, varargin{:}, 'InterleaverIndicesSource', 'InterleaverIndices', 'NumIterations');
            
            this.trellis = poly2trellis(4, [13 15], 13);

            if strcmp(this.InterleaverIndicesSource, 'Property')
                this.CoderObj = comm.TurboEncoder( ...
                    'TrellisStructure', this.trellis, ...
                    'InterleaverIndicesSource', this.InterleaverIndicesSource, ...
                    'InterleaverIndices', this.InterleaverIndices);

                this.DecoderObj = comm.TurboDecoder( ...
                    'TrellisStructure', this.trellis, ...
                    'InterleaverIndicesSource', this.InterleaverIndicesSource, ...
                    'InterleaverIndices', this.InterleaverIndices, ...
                    'Algorithm', 'Max', ...
                    'NumIterations', this.NumIterations);  
            elseif strcmp(this.InterleaverIndicesSource, 'Input port')
                this.CoderObj = comm.TurboEncoder( ...
                    'TrellisStructure', this.trellis, ...
                    'InterleaverIndicesSource', this.InterleaverIndicesSource);

                this.DecoderObj = comm.TurboDecoder( ...
                    'TrellisStructure', this.trellis, ...
                    'InterleaverIndicesSource', this.InterleaverIndicesSource, ...
                    'Algorithm', 'Max', ...
                    'NumIterations', this.NumIterations);  
            end
            
            
            this.AppDecObj1 = comm.APPDecoder(...
                'TrellisStructure',this.trellis, ...
                'TerminationMethod', 'Terminated',...
                'Algorithm','Max', ...
                'CodedBitLLROutputPort',true);
            this.AppDecObj2 = comm.APPDecoder(...
                'TrellisStructure',this.trellis, ...
                'TerminationMethod', 'Terminated',...
                'Algorithm','Max', ...
                'CodedBitLLROutputPort',true);
        end
        
        function out = encode(this, in, varargin)
            if strcmp(this.InterleaverIndicesSource, 'Property')
                out = this.CoderObj.step(in);
            elseif strcmp(this.InterleaverIndicesSource, 'Input port')
                out = this.CoderObj.step(in, varargin{1});
            end
        end
        
        function out = decode(this, in, varargin)
            % Avoid NaN due to calculation of infinity numbers
            in(in==Inf) = 700;
            in(in==-Inf) = -700;            
            if strcmp(this.InterleaverIndicesSource, 'Property')
                out = this.DecoderObj.step(in);
            elseif strcmp(this.InterleaverIndicesSource, 'Input port')
                out = this.DecoderObj.step(in, varargin{1});
            end
        end
        
        function [outLapr, inLapr, outHard, error] = decodeWithUpdatedLlrOutput(this, in, varargin)  % vargin is used when the early termination is enabled
            if strcmp(this.InterleaverIndicesSource, 'Property')
                interleaverIndices = this.InterleaverIndices;
                if nargin > 2
                    Crc24Coder = varargin{1};
                end
            elseif strcmp(this.InterleaverIndicesSource, 'Input port')
                interleaverIndices = varargin{1};
                if nargin  > 3
                    Crc24Coder = varargin{2};
                end
            end

            % Avoid NaN due to calculation of infinity numbers
            in(in==Inf) = 700;
            in(in==-Inf) = -700;
            
            tailBetLength = log2(this.trellis.numStates) * log2(this.trellis.numOutputSymbols);
            K = (length(in) - 2*tailBetLength) / 3; % K = transport block size
        
            systematic = in(1:3:3*K);
            parity1 = [in(2:3:3*K); in(3*K+1:3*K+tailBetLength)];
            parity2 = [in(3:3:3*K); in(3*K+tailBetLength+1:end)];

            in1 = [systematic parity1(1:K)].';
            in1 = [in1(:); parity1(K+1:end)];
            in2 = [intrlv(systematic, interleaverIndices) parity2(1:K)].';
            in2 = [in2(:); parity2(K+1:end)];

            outLapr = zeros(K+tailBetLength/2,1); % Initial value of a priori information 
            for iter = 1:this.NumIterations
                [outLapr1, inLapr1] = step(this.AppDecObj1, outLapr, in1);
                outLapr(1:K) = intrlv(outLapr1(1:K), interleaverIndices);
                [outLapr2, inLapr2] = step(this.AppDecObj2, outLapr, in2);
                outLapr(1:K) = deintrlv(outLapr2(1:K), interleaverIndices);
                
                if nargin > 2  % Early termination
                    outLaprTmp = outLapr(1:K);
                    outHard = zeros(size(outLaprTmp));
                    outHard(outLaprTmp>0) = 1; 
                    outHard(outLaprTmp<=0) = 0; 
                    [~, error] = Crc24Coder.detect(outHard);
                    if ~error
                        break;
                    end
                end
                
                outLapr(outLapr==Inf) = 700;
                outLapr(outLapr==-Inf) = -700;
            end
            
            % Truncate tails
            outLapr = outLapr(1:K);
            % change order according to [systematic parity1 parity2]
            systematic2 = deintrlv(inLapr2(1:2:2*K), interleaverIndices);
            systematic = inLapr1(1:2:2*K) + systematic2;
%             parity1 = inLapr1(2:2:2*K);
%             parity2 = inLapr2(2:2:2*K);
            tmp = [systematic inLapr1(2:2:2*K) inLapr2(2:2:2*K)].';
            inLapr = [tmp(:); inLapr1(2*K+1:end); inLapr2(2*K+1:end)]; % Add tails
        end
    end
end
