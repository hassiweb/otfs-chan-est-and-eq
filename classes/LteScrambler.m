% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef LteScrambler < matlab.System
    
    properties
        hSeqGen;
        hInt2Bit;
        initialStates1;
        initialStates2;
    end
        
    methods
        function this = LteScrambler(rnti, cellId, frameNo, codeword) % Constructor
            if nargin > 0
                % Create Gold sequence
                maxG=43200;
                this.hSeqGen = comm.GoldSequence( ...
                    'FirstPolynomial',[1 zeros(1, 27) 1 0 0 1],...
                    'FirstInitialConditions', [zeros(1, 30) 1], ...
                    'SecondPolynomial', [1 zeros(1, 27) 1 1 1 1],...
                    'SecondInitialConditionsSource', 'Input port',... 
                    'Shift', 1600,...
                    'VariableSizeOutput', true,...
                    'MaximumOutputSize', [maxG 1]);
                this.hInt2Bit = comm.IntegerToBit('BitsPerInteger', 31);

                % Initial conditions to create Gold sequence
                slotNo = [frameNo*2, frameNo*2+1];
                c_init = rnti*(2^14) + codeword*(2^13) + floor(slotNo/2)*(2^9) + cellId; % For first slot

                % Convert initial condition to binary vector
                this.initialStates1 = step(this.hInt2Bit, c_init(1));
                this.initialStates2 = step(this.hInt2Bit, c_init(2));
            end
        end
        
        function scrambledSeq = scramble(this, bitSeq)
            % Generate the scrambling sequence
            nSamples    = length(bitSeq)/2;
            seq_tmp1    = zeros(nSamples,1);
            seq_tmp2    = zeros(nSamples,1);
            seq_tmp1(:) = step(this.hSeqGen, this.initialStates1, nSamples);
            seq_tmp2(:) = step(this.hSeqGen, this.initialStates2, nSamples);
            seq         = [seq_tmp1; seq_tmp2];

            % Scramble input with the scrambling sequence
            scrambledSeq = double(xor(bitSeq, seq));
            scrambledSeq = logical(scrambledSeq);
        end
        
        function descrambledSeq = descramble(this, bitSeq)
            % Generate the scrambling sequence
            nSamples    = length(bitSeq)/2;
            seq_tmp1    = zeros(nSamples,1);
            seq_tmp2    = zeros(nSamples,1);
            seq_tmp1(:) = step(this.hSeqGen, this.initialStates1, nSamples);
            seq_tmp2(:) = step(this.hSeqGen, this.initialStates2, nSamples);
            seq         = [seq_tmp1; seq_tmp2];

            seq2    = zeros(size(bitSeq));
            seq2(:) = seq(1:numel(bitSeq), 1);

            % If descrambler inputs are log-likelihood ratios (LLRs) then convert
            % sequence to a bipola format
            seq2 = 1-2.*seq2;

            % Descramble
            descrambledSeq = bitSeq.*seq2;
            %     descrambledSeq = xor(bitSeq(:,1), seq(:,1)); % hard detection
            
        end
        
    end
end
        