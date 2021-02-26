% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef OtfsPilotResponseBasedPathParameterEstimator < matlab.System
    properties
        DividingNumber;
        CyclicPrefixLength;
        NumDopplerBins;
        NumDelayBins;
        SamplingRate;
        
        upsilon;
    end
    
    methods
        function this = OtfsPilotResponseBasedPathParameterEstimator(varargin) % Constructor
            setProperties(this, nargin, varargin{:}, ...
                'DividingNumber', ...
                'CyclicPrefixLength', ...  % Ncp
                'NumDopplerBins', ...      % N
                'NumDelayBins', ...        % M
                'SamplingRate' ...
                );
            
            kappa = 0:1/this.DividingNumber:1-1/this.DividingNumber;
            
            this.upsilon = zeros(this.NumDopplerBins,length(kappa));
            for k = 1:length(kappa)
                x = kappa(k) - (0:this.NumDopplerBins-1);
                for i = 1:this.NumDopplerBins
                    this.upsilon(:,k) = this.upsilon(:,k) + 1/this.NumDopplerBins * exp(1i*2*pi*(i-1)*x/this.NumDopplerBins).';
                end
            end
        end
        
        
        function [estGains, estDopplers, estDelays, estPhaseOffsets] = estimate(this, Hdd, alpha, beta, noiseVar)
            % ----- DESCRIPTION -----
            % This function estimates path parameters from the pilot response in the delay Doppler domain based on [1].
            % ----- INPUT PARAMETERS -----
            %   - Hdd: pilot response in the delay-Doppler domain (size: M x N)
            %   - alpha: tuning parameter for the summation of magnitude
            %   - beta: tuning parameter for noise variance
            %   - noiseVar: noise variance
            % ----- REFERENCE -----
            % [1] https://arxiv.org/abs/2010.15396

            M = this.NumDelayBins;
            N = this.NumDopplerBins;
            Ncp = this.CyclicPrefixLength;
            
            estGains = [];
            estDopplers = [];
            estDelays = [];
            estPhaseOffsets = [];
            Hdd = Hdd/sqrt(M*N);

            kappa = 0:1/this.DividingNumber:1-1/this.DividingNumber;
            pathindex = 1;
            for delay = 0:M-1
                Hprime = Hdd(delay+1,:);
                sumMagnitude = abs(sum(Hdd(delay+1,:)));  % Sum of amplitudes in this delay bin
                estGain = Inf;
                while(1)
                    % Take cross-correlation b/w Hdd and Upsilon
                    crosscorr = zeros(length(kappa), N);
                    for idx_kappa = 1:length(kappa)
                        crosscorr(idx_kappa,:) = ifft( fft(Hprime) .* conj(fft(this.upsilon(:,idx_kappa).')) );
                    end
                    crosscorrvec = circshift(crosscorr(:), N*this.DividingNumber/2);

                    % Condition 1: find paths until the cross-correlation becomes small
                    [corrval, largestindices] = sort(abs(crosscorrvec), 'descend');
                    largestindices(corrval < alpha * sumMagnitude) = [];  % remove indices that are less then the threshold

                    % Condition 2: the cross-correlation is not negligibly small (equally or less than the noise power)
                    largestindices(abs(crosscorrvec(largestindices)) < beta * sqrt(noiseVar)) = [];  % remove indices that are less then the noise power

                    if isempty(largestindices)
                        break;
                    end

                    largestidx = largestindices(1);

                    if estGain < abs(crosscorrvec(largestidx))
                        break;
                    end
                    estGain = abs(crosscorrvec(largestidx));
                    estDoppler = (largestidx-1)/this.DividingNumber-N/2;  % index -> value of Doppler
                    estDopplerShift = exp(1i*2*pi*estDoppler*(Ncp-delay)/((M+Ncp)*N));
                    estInitialPhase = crosscorrvec(largestidx)/abs(crosscorrvec(largestidx))*estDopplerShift^-1;

                    estGains = [estGains estGain];
                    estDopplers = [estDopplers estDoppler];
                    estDelays = [estDelays delay];
                    estPhaseOffsets = [estPhaseOffsets angle(estInitialPhase)];

%                     % Check estimated paths
%                     figure; clf;
%                     plot((0:1/this.DividingNumber:N-1/this.DividingNumber)-N/2, abs(crosscorrvec)); hold on; 
%                     stem(estDopplers(pathindex), estGains(pathindex));
%                     xlim([-N/2 N/2]); ylim([0 1]);
%                     xlabel('k+\kappa'); ylabel('|R_H_,_\Upsilon(k+\kappa)|');

                    intDoppler = floor(estDopplers(pathindex));
                    fracDopplerIdx = floor(mod(estDopplers(pathindex),1)*this.DividingNumber + 1 + 1e-9);
                    Hprime = Hprime - estGains(pathindex) * exp(1i*estPhaseOffsets(pathindex))*circshift(this.upsilon(:,fracDopplerIdx),intDoppler).';

                    pathindex = pathindex + 1;
                end
            end
            
            estDopplers = estDopplers / N / ((M+Ncp)/this.SamplingRate);
            estDelays = estDelays / this.SamplingRate;
%             estPhaseOffsets = estPhaseOffsets - angle(exp(1i*2*pi*estDopplers*(M+Ncp)*N/this.SamplingRate));
        end
        
       
        function estError = evaluateEstimationError(this, HddIdeal, estGains, estDopplers, estDelays, estOffsets)
            HddIdeal = HddIdeal/sqrt(this.NumDopplerBins*this.NumDelayBins);
            numPaths = length(estGains);
            estDelaysInSample = estDelays*this.SamplingRate;
            estDopplersInSample = estDopplers/(this.SamplingRate/((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins)); 
            
            HddEst = zeros(this.NumDelayBins, this.NumDopplerBins);
            for p = 1:numPaths
                for l = 0:this.NumDelayBins-1
                    if estDelaysInSample(p) == l
                        Upsilon_N = zeros(1,this.NumDopplerBins);
                        for i = 1:this.NumDopplerBins
                            x = (estDopplers(p)/(this.SamplingRate/((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins))-[0:this.NumDopplerBins-1]);
                            Upsilon_N = Upsilon_N + 1/this.NumDopplerBins * exp(1i*2*pi*(i-1)*x/this.NumDopplerBins);
                        end                        
                        HddEst(l+1,:) = HddEst(l+1,:) + estGains(p) * exp(1i*estOffsets(p)) * exp(1i*2*pi*estDopplersInSample(p)*(this.CyclicPrefixLength-estDelaysInSample(p)+l)/((this.NumDelayBins+this.CyclicPrefixLength)*this.NumDopplerBins)) * Upsilon_N;
                    end
                end
            end
            
            estError = 1/(this.NumDopplerBins*this.NumDelayBins) * sum( abs( HddIdeal(:) - HddEst(:) ).^2 );
        end        
    end
end

