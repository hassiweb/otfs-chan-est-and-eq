% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef OtfsPNSeqBasedPathParameterEstimator < matlab.System
    properties
        DividingNumber;
        CyclicPrefixLength;
        NumDopplerBins;
        NumDelayBins;
        SamplingRate;
        Threshold;
    end
    
    methods
        function this = OtfsPNSeqBasedPathParameterEstimator(varargin) % Constructor
            setProperties(this, nargin, varargin{:}, ...
                'DividingNumber', ...
                'CyclicPrefixLength', ...  % Ncp
                'NumDopplerBins', ...      % N
                'NumDelayBins', ...        % M
                'SamplingRate', ...
                'Threshold' ...
                );
        end
        

        function [estGains, estDopplers, estDelays, estPhaseOffsets] = estimate(this, rxPNSeq, txPNSeq, varargin)
            % ----- DESCRIPTION -----
            % This function estimates path parameters composed of the channel from the received PN sequence.
            % The path can be detected when the magnitude of the inner product b/w the received and transmitted sequences is higher than the threshold.
            % However, this estimation can work well when the length of the PN sequence is relatively long to the Doppler frequency.
            % ----- INPUT PARAMETERS -----
            %   - rxPNSeq: the received PN sequence
            %   - txPNSeq: the transmitted PN sequence
            %   - varargin:
            %     - expectedMaxDoppler: expected max Doppler to reduce computations by limiting the range of Doppler to calculate the inner products

            matchedFilterMatrix = zeros(this.NumDopplerBins*this.DividingNumber, this.CyclicPrefixLength);
            filterLength = length(rxPNSeq)-this.NumDelayBins;
            t = (0:filterLength-1)'/this.SamplingRate;

            % Identify a Doppler range to search
            if isempty(varargin)  % if varargin is empty, 
                dopplerSearchRange = 0:this.NumDopplerBins*this.DividingNumber-1; % search Doppler shifts within all available values given by the system
            else
                expectedMaxDoppler = varargin{1};
                maxDopplerBinIndex = ceil(expectedMaxDoppler * this.DividingNumber / this.SamplingRate * ((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins));
                dopplerSearchRange = (-maxDopplerBinIndex:maxDopplerBinIndex-1)+(this.NumDopplerBins*this.DividingNumber/2); % search Doppler shifts within given max Doppler shift
            end
            
            % Mached filter (inner product b/w the received and transmitted PN sequence)
%             for delay = 0:this.NumDelayBins-1  % search delays within all available values given by the system
            for delay = 0:this.CyclicPrefixLength-1  % search delays within the given cyclic prefix
                for omegaIdx = dopplerSearchRange  
                    fDoppler = (omegaIdx-this.NumDopplerBins*this.DividingNumber/2)/this.DividingNumber*this.SamplingRate / ((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins);
                    matchedFilterMatrix(omegaIdx+1, delay+1) = dot(exp(1i*2*pi*fDoppler*t).*txPNSeq(1:filterLength), rxPNSeq(1+delay:filterLength+delay))/filterLength;
                end
            end
            
            % Pick indices that have higher values than the threshold, i.e., picked indices are estimated paths
            estpathindices = find(abs(matchedFilterMatrix) > this.Threshold);

            estGains = zeros(length(estpathindices),1);
            estDopplers = zeros(length(estpathindices),1);
            estDelays = zeros(length(estpathindices),1);
            estPhaseOffsets = zeros(length(estpathindices),1);
            for p = 1:length(estpathindices)
                estGains(p) = abs(matchedFilterMatrix(estpathindices(p)));
                estDopplers(p) = (mod(estpathindices(p), this.NumDopplerBins*this.DividingNumber) - this.NumDopplerBins*this.DividingNumber/2 - 1)/this.DividingNumber;
                estDelays(p) = floor(estpathindices(p)/(this.NumDopplerBins*this.DividingNumber));
                estPhaseOffsets(p) = angle(matchedFilterMatrix(estpathindices(p)));
            end
            estDopplers = estDopplers * this.SamplingRate / ((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins);
            estDelays = estDelays / this.SamplingRate;
        end

        function [estGains, estDopplers, estDelays, estPhaseOffsets] = estimateIteratively(this, rxPNSeq, txPNSeq, varargin)
            % ----- DESCRIPTION -----
            % This function estimates path parameters composed of the channel from the received PN sequence.
            % The highest magnitude of inner product b/w the transmitted and received PN sequence with a delay and Doppler parameter is an estimated path.
            % After estimating a path, the path with parameters will be removed from the received PN sequence.
            % Then, this operation will be repeated until the magnitude is smaller than the threshold.
            % Note: This algorithm is just an idea, so this may NOT be formulated.
            % This algorithm has better performance than the threshold-based algorithm, but this algorithm requires high number of computations compared to other algorithms at the least.
            % ----- INPUT PARAMETERS -----
            %   - rxPNSeq: the received PN sequence
            %   - txPNSeq: the transmitted PN sequence
            %   - varargin:
            %     - expectedMaxDoppler: expected max Doppler to reduce computations by limiting the range of Doppler to calculate the inner products

            matchedFilterMatrix = zeros(this.NumDopplerBins*this.DividingNumber, this.NumDelayBins);
            filterOutputLength = length(rxPNSeq)-this.NumDelayBins;
            t = (0:filterOutputLength-1)'/this.SamplingRate;

            % Identify a Doppler range to search
            if isempty(varargin)  % if varargin is empty, 
                dopplerSearchRange = 0:this.NumDopplerBins*this.DividingNumber-1; % search Doppler shifts within all available values given by the system
            else
                expectedMaxDoppler = varargin{1};
                maxDopplerBinIndex = ceil(expectedMaxDoppler * this.DividingNumber / this.SamplingRate * ((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins));
                dopplerSearchRange = (-maxDopplerBinIndex:maxDopplerBinIndex-1)+(this.NumDopplerBins*this.DividingNumber/2); % search Doppler shifts within given max Doppler shift
            end
            
            % Storages for estimated parameters
            estGains = [];
            estDopplers = [];
            estDelays = [];
            estPhaseOffsets = [];

            
            residualRxPNSeq = rxPNSeq;
            while (1)
                % Matched filter (inner products between the received and transmitted sequences)
                for delay = 0:this.CyclicPrefixLength-1  % search delays within the given cyclic prefix
                    for omegaIdx = dopplerSearchRange
                        fDoppler = (omegaIdx-this.NumDopplerBins*this.DividingNumber/2)/this.DividingNumber*this.SamplingRate / ((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins);
                        matchedFilterMatrix(omegaIdx+1, delay+1) = dot(exp(1i*2*pi*fDoppler*t).*txPNSeq(1:filterOutputLength), residualRxPNSeq(1+delay:filterOutputLength+delay))/filterOutputLength;
                    end
                end

                % Find the highest magnitude of the inner product while the magnitude is higher than the threshold
                maxpathvalue = max(abs(matchedFilterMatrix(:)));
                if maxpathvalue < this.Threshold
                    break;
                end

                % Pick the highest magnitude of the inner product as an estimated path
                estpathindex = find(abs(matchedFilterMatrix)==maxpathvalue);
                estGain = abs(matchedFilterMatrix(estpathindex));
                estDoppler = (mod(estpathindex, this.NumDopplerBins*this.DividingNumber) - this.NumDopplerBins*this.DividingNumber/2 - 1)/this.DividingNumber;
                estDelay = floor(estpathindex/(this.NumDopplerBins*this.DividingNumber));
                estPhase = angle(matchedFilterMatrix(estpathindex));
                
                % Remove the estimated path from the received PN sequence
                fDoppler = estDoppler*this.SamplingRate/((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins);
                delayedtime = (0:length(rxPNSeq)-1-estDelay)'/this.SamplingRate;
                replica = txPNSeq(1:end-estDelay) .* (estGain*exp(1i*(2*pi*fDoppler*delayedtime + estPhase)));
                residualRxPNSeq(1+estDelay:end) = residualRxPNSeq(1+estDelay:end) - replica;
                
                % Store the estimated path parameters
                estGains = [estGains estGain];
                estDopplers = [estDopplers estDoppler];
                estDelays = [estDelays estDelay];
                estPhaseOffsets = [estPhaseOffsets estPhase];
                
            end
            estDopplers = estDopplers * this.SamplingRate / ((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins);
            estDelays = estDelays / this.SamplingRate;
            
        end
        
    end
end