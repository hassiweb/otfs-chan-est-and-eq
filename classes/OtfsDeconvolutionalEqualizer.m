% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef OtfsDeconvolutionalEqualizer < matlab.System
    properties
        EqualizationAlgorithm;
        NumSymbols;
        CyclicPrefixLength;
        NumSubcarriers;
        SamplingRate;
        DividingNumber;
        
        upsilon;
    end
    
    methods
        function this = OtfsDeconvolutionalEqualizer(varargin) % Constructor
            setProperties(this, nargin, varargin{:}, ...
                'NumSymbols', ...          % N
                'CyclicPrefixLength', ...  % Ncp
                'NumSubcarriers', ...      % M
                'SamplingRate', ...
                'DividingNumber', ...      % this determines kappa
                'EqualizationAlgorithm' ...  % ZF or Wiener
                );
            
            % Pre-calculation of kappa
            kappa = 0:1/this.DividingNumber:1-1/this.DividingNumber;

            % Pre-calculation of Upsilon
            this.upsilon = zeros(this.NumSymbols,length(kappa));
            for k = 1:length(kappa)
                x = kappa(k) - (0:this.NumSymbols-1);
                for i = 1:this.NumSymbols
                    this.upsilon(:,k) = this.upsilon(:,k) + 1/this.NumSymbols * exp(1i*2*pi*(i-1)*x/this.NumSymbols).';
                end
            end
        end
        
        function Xdd = equalize(this, Ydd, pathGains, pathDopplers, pathDelays, pathOffsets, noiseVar, varargin)
            % ----- INPUT PARAMETERS -----
            %   - Ydd: received symbols in the delay Doppler domain
            %   - pathGains, pathDopplers, pathDelays, pathOffsets: path parameters
            %   - noiseVar: noise variance
            %   - (EXPERIMENTAL) varargin: threshold to avoid zero division or too small value to divide if this argument is inputted.  (This has not been well tested)

            if isempty(varargin)  % if varargin is empty, 
                threshold = 0;
            else
                threshold = varargin{1};
            end
            
            P = length(pathGains);
            M = this.NumSubcarriers;
            N = this.NumSymbols;
            Ncp = this.CyclicPrefixLength;
            
            % Convert delay and Doppler into grid indices
            pathDelayInSample = pathDelays * this.SamplingRate;
            pathDopplerInSample = pathDopplers * N * ((M+Ncp)/this.SamplingRate);

            % Calc Lambda (an effective channel matrix)
            lambda = zeros(P,M,N);
            for p = 1:P
                for l = 0:M-1
                    intDoppler = mod(floor(round(mod(pathDopplerInSample(p),N)+1e-9, 1)),N);
                    fracDopplerIndex = mod(round((mod(pathDopplerInSample(p),N)-intDoppler)*this.DividingNumber+1e-9), this.DividingNumber) + 1;
                    hp = pathGains(p) * (pathDelayInSample(p)==l);
                    lambda(p,l+1,:) = hp * circshift(this.upsilon(:,fracDopplerIndex), intDoppler);
                end
            end
            
            Ydd_dft = fft2(Ydd);
            Xdd = zeros(M,N);
            for l = 0:M-1
                psi = exp(1i*2*pi*pathDopplerInSample.*(Ncp-pathDelayInSample+l)/((M+Ncp)*N));
                phi = exp(1i*pathOffsets);
                Hdd = sum(psi .* phi .* lambda);
                Hdd = permute(Hdd, [2,3,1]);
                
                Hdd_dft = fft2(Hdd);

                switch this.EqualizationAlgorithm
                    case 'ZF'
                        Xdd_tmp = ifft2(Ydd_dft ./ Hdd_dft);
                    case 'Wiener'
%                         eq = conj(Hdd_dft) ./ (abs(Hdd_dft).^2 + noiseVar);  % https://www.researchgate.net/publication/223015617_The_point_spread_function_revisited_Image_restoration_using_2-D_deconvolution

                        denominator = (abs(Hdd_dft).^2 + noiseVar);

                        % %%%%% EXPERIMENTAL 1: add small value to the denominator %%%%%%
%                         eq = conj(Hdd_dft) ./ (denominator + threshold);  % https://www.researchgate.net/publication/223015617_The_point_spread_function_revisited_Image_restoration_using_2-D_deconvolution

                        % %%%%% EXPERIMENTAL 2: local modification to avoid zero (and small value) division %%%%%%
%                         if not(isempty(denominator(denominator<threshold)))
%                             denominator = denominator + threshold;
%                         end
%                         eq = conj(Hdd_dft) ./ denominator;

                        % %%%%% EXPERIMENTAL 3: global modification (add estimated variance of the estimation error to the denominator) to avoid zero (or small value) division %%%%%%
                        if not(isempty(denominator(denominator<threshold)))
                            Pt = sqrt(M*N)*0.9472; % for 16QAM
                            Pr = mean(abs(Ydd_dft(:)));
                            Ph = (Pr - noiseVar*sqrt(M*N));
                            Ph_est = mean(abs(Hdd_dft(:)));
                            est_err = (Ph_est*Pt - Ph)/sqrt(M*N)/2;
                            denominator = denominator + est_err;
                        end
                        eq = conj(Hdd_dft) ./ denominator;

                        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        Xdd_tmp = ifft2(Ydd_dft .* eq);
                end
                Xdd(l+1,:) = Xdd_tmp(l+1,:);
            end
%             figure(1); clf; plot(Xdd, '.'); grid on; grid minor; xlim([-1.5 1.5]); ylim([-1.5 1.5]); % hold on; plot(rxIdealEqBlocks, '.'); hold off;
        end
        
        function Xdd = equalizeUsingPilotResponse(this, Ydd, Hdd)
            % ----- INPUT PARAMETERS -----
            %   - Ydd: received symbols in the delay Doppler domain
            %   - Hdd: channel impulse response that is captured by a pilot symbol in the delay Doppler domain

            M = this.NumSubcarriers;
            N = this.NumSymbols;
            Hdd = 1/sqrt(M*N)*Hdd;
 
            Xdd = zeros(M,N);
            for l = 0:M-1
                Xdd_tmp = ifft2(fft2(Ydd) ./ fft2(Hdd));
                Xdd(l+1,:) = Xdd_tmp(l+1,:);
            end
        end        
    end
end