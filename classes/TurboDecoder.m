% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef TurboDecoder
    
    properties
        algorithm;
        trellis;
        interleaver;
        maxIter;
        AppDecObj1;
        AppDecObj2;
    end
    
    methods
        function this = TurboDecoder(trellis, interleaver, maxIter)
            this.trellis = trellis;
            this.interleaver = interleaver;
            this.maxIter = maxIter;
            
            this.AppDecObj1 = comm.APPDecoder(...
                'TrellisStructure',trellis, ...
                'TerminationMethod', 'Terminated',...
                'Algorithm','Max', ...
                'CodedBitLLROutputPort',true);
            this.AppDecObj2 = comm.APPDecoder(...
                'TrellisStructure',trellis, ...
                'TerminationMethod', 'Terminated',...
                'Algorithm','Max', ...
                'CodedBitLLROutputPort',true);
        end
        
        function [outLapr, inLapr] = decode(this, in)
            
            numTails = log2(this.trellis.numStates) * log2(this.trellis.numOutputSymbols);
            K = (length(in) - 2*numTails) / 3; % K = transport block size
        
            systematic = in(1:3:3*K);
            parity1 = [in(2:3:3*K); in(3*K+1:3*K+numTails)];
            parity2 = [in(3:3:3*K); in(3*K+numTails+1:end)];

            in1 = [systematic parity1(1:K)].';
            in1 = [in1(:); parity1(K+1:end)];
            in2 = [systematic(this.interleaver) parity2(1:K)].';
            in2 = [in2(:); parity2(K+1:end)];

            outLapr = zeros(K+numTails/2,1); % Initial value of a priori information 
            for iter = 1:this.maxIter
                [outLapr1, inLapr1] = step(this.AppDecObj1, outLapr, in1);
%                 disp(["a" sum(abs((-aEstLlrTq1(1:K)>=0) - a))])

                outLapr(1:K) = outLapr1(this.interleaver);

                [outLapr2, inLapr2] = step(this.AppDecObj2, outLapr, in2);
%                 disp(["a" sum(abs((-aEstLlrTq2(1:K)>=0) - a(this.interleaver)))])

                outLapr(this.interleaver) = outLapr2(1:K);

                %Clip
%                 LaprInner(LaprInner> maxLLR) = maxLLR;
%                 LaprInner(LaprInner<-maxLLR) =-maxLLR;      
            end
            
            % Truncate tails
            outLapr = outLapr(1:K);

            % change order according to [sys par1 par2]
            systematic2 = zeros(K,1);
            systematic2(this.interleaver) = inLapr2(1:2:2*K);
            systematic = inLapr1(1:2:2*K) + systematic2;
    %         parity1 = [bEstLlrTq1(2:2:2*K); bEstLlrTq1(2*K+1:end)];
    %         parity2 = [bEstLlrTq2(2:2:2*K); bEstLlrTq2(2*K+1:end)];
            tmp = [systematic inLapr1(2:2:2*K) inLapr2(2:2:2*K)].';
            inLapr = [tmp(:); inLapr1(2*K+1:end); inLapr2(2*K+1:end)];
%             inLapr = inLapr(outerInterleaver);

            %Clip
%             LaprOuter(LaprOuter> maxLLR) = maxLLR;
%             LaprOuter(LaprOuter<-maxLLR) =-maxLLR;  
            
            % Avoid overflow
            inLapr(inLapr==Inf) = 300;
            inLapr(inLapr==-Inf) = -300;  
            outLapr(outLapr==Inf) = 300;
            outLapr(outLapr==-Inf) = -300;  
        end
    end
end