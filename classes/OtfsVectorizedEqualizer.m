% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef OtfsVectorizedEqualizer < matlab.System
    properties
        NumSymbols;
        CyclicPrefixLength;
        NumSubcarriers;
        SamplingRate;
        EqualizationAlgorithm;
        OutputConditionNumber;
        
        spKronDftMat;
        kronDftMat;
        spKronIdftMat;
    end
    
    methods
        function this = OtfsVectorizedEqualizer(varargin) % Constructor
            setProperties(this, nargin, varargin{:}, ...
                'EqualizationAlgorithm', ...
                'NumSymbols', ...          % N
                'CyclicPrefixLength', ...  % Ncp
                'NumSubcarriers', ...      % M
                'SamplingRate', ...
                'OutputConditionNumber', ...
                'EqualizationAlgorithm' ...  % ZF or MMESE
                );

            % Fourier transform matrix
%             F_M = exp(-1i*2*pi*(0:this.NumSubcarriers-1).'*(0:this.NumSubcarriers-1)/(this.NumSubcarriers));
            F_N = exp(-1i*2*pi*(0:this.NumSymbols-1).'*(0:this.NumSymbols-1)/this.NumSymbols);

            switch this.EqualizationAlgorithm
                case {'ZF', 'MMSE'}
                    this.spKronDftMat = sparse(kron(F_N, eye(this.NumSubcarriers)));
                    this.spKronIdftMat = sparse(kron(F_N', eye(this.NumSubcarriers)));
                case {'blockedMMSE'}
                    this.kronDftMat = kron(F_N, eye(this.NumSubcarriers));
            end
            
            if isempty(this.OutputConditionNumber)
                this.OutputConditionNumber = false;
            end
        end
        
        function [x, conditionNumber] = equalize(this, Ydd, pathGains, pathDopplers, pathDelays, pathOffsets, noiseVar, varargin)
            % ----- DESCRIPTION -----
            % This function re-generates the channel matrix using path parameters.
            % Then, the received symbols are equalized by multiplying the inverse matrix of the generated channel matrix.
            % ----- INPUT PARAMETERS -----
            %   - Ydd: received symbols in the Delay doppler domain
            %   - pathGains, pathDopplers, pathDelays, pathOffsets: path parameters
            %   - noiseVar: noise variance for MMSE
            
            % Convert fractional delays into integer delays
            [fdfedPathGains, fdfedPathDopplers, fdfedPathDelays, fdfedPathOffsets] = fdfpathpars(pathGains, pathDopplers, pathDelays, pathOffsets, this.SamplingRate);
            
            % Initialize parameters to be short names
            N = this.NumSymbols;
            Ncp = this.CyclicPrefixLength;
            M = this.NumSubcarriers;
            P = length(fdfedPathGains);

            % Delay-Doppler domain grid representation
            fdfedPathDelaysInSample = fdfedPathDelays*this.SamplingRate;
            fdfedPathDopplersInSample = fdfedPathDopplers * N * ((M+Ncp)/this.SamplingRate);

            %% Re-generate OFDM-symbol-wise (i: symbol number) time domain channel matrix: Hi (size: M x M)
            % Initialize matrixes
            Hi = cell(N,1);
            for i = 1:N
                Hi{i} = zeros(M,M);
            end

            % Re-generate Hi
            for p = 1:P
                for i = 1:N
                    timeIndices = (((M+Ncp)*(i-1)+Ncp) : ((M+Ncp)*i-1)) - fdfedPathDelaysInSample(p);
                    phaseShift = exp(1i * (2*pi * fdfedPathDopplersInSample(p) / ((M+Ncp)*N) * timeIndices)); 
                    dopplerMat = diag(phaseShift);
                    Hi{i} = Hi{i} + fdfedPathGains(p) * exp(1i*fdfedPathOffsets(p)) * circshift(dopplerMat, -fdfedPathDelaysInSample(p), 2); 
                end
            end

            switch this.EqualizationAlgorithm
                case {'ZF', 'MMSE'}
                    % Combine channel matrixes into a block-wise channel matrix: H (size: MN x MN)
                    H = zeros(M*N, M*N);
                    for i = 1:N
                        H((i-1)*M+1:i*M,(i-1)*M+1:i*M) = Hi{i};
                    end
                    spH = sparse(H);  % convert into a sparse matrix to reduce the amount of computation and memory

                    % Convert H into the time domain channel matrix in the delay-Doppler domain: Hdd (a.k.a. phi in the document)
%                     Hdd = 1/N * kron(this.F_N, eye(M)) * H * kron(conj(this.F_N), eye(M));  % non-sparse matrix version
                    spHdd = 1/N * this.spKronDftMat * spH * this.spKronIdftMat;  % This is much faster than above.
            end
            
            %% Equalize
            switch this.EqualizationAlgorithm
                case 'ZF'
                    Hdd = full(spHdd);  % Sparse matrix -> non-sparse matrix
                    x = Hdd\Ydd(:);  % Zero forcing
                    if this.OutputConditionNumber
                        conditionNumber = cond(Hdd);
                    end
                case 'MMSE'  % Performance is worse than ZF. Why?
                    tmp1 = full(spHdd*spHdd'+noiseVar*speye(M*N));
                    tmp2 = full(spHdd');
                    eq = tmp1\tmp2;
                    x = eq*Ydd(:);
                    if this.OutputConditionNumber
                        conditionNumber = cond(tmp1);
                    end
                case 'blockedMMSE'
                    % This implementation is based on: S. S. Das, V. Rangamgari, S. Tiwari, and S. C. Mondal, ÅgTime Domain Channel Estimation and Equalization of CP-OTFS Under Multiple Fractional Dopplers and Residual Synchronization Errors,Åh IEEE Access, pp. 1?1, 2020, doi: 10.1109/ACCESS.2020.3046487.
                    % This can reduce the computational complexity of the matrix inversion.  It will also improve the condition number of the channel matrix.

                    r = ifft(permute(isfft2d(Ydd),[2,1,3])); % time-domain received signal
                    rce_i = cell(N,1);
                    for i = 1:N
                        Pt = 1/M;
                        Pr = diag(abs(r(:,i)));
                        Ph = (Pr - diag(noiseVar))/M;
                        Ph_est = diag(Hi{i}*Hi{i}');
                        cov = diag(Ph_est).*Pt - Ph;
                        rce_i{i} = Hi{i}' / (Hi{i}*Hi{i}'+cov+noiseVar*eye(M)) * r(:,i);
                    end
                    
                    rce = cell2mat(rce_i);
                    x = this.kronDftMat*rce;
                    x = x/rms(x);
                    if this.OutputConditionNumber
                        conditionNumber = cond(Hi{i}*Hi{i}'+noiseVar*(eye(M)));
                    end                    
                case 'Eigen'
                    % This implementation is based on :[1]G. D. Surabhi and A. Chockalingam, ìLow-Complexity Linear Equalization for OTFS Modulation,? IEEE Communications Letters, vol. 24, no. 2, pp. 330?334, Feb. 2020, doi: 10.1109/LCOMM.2019.2956709.
                    % However, currently this doesn't work when using the propose channel estimator.
                    Hdd = full(spHdd);
                    eigH = eig(Hdd);
                    lambda = diag(eigH);
                    psi = inv(conj(lambda)*lambda+noiseVar*eye(length(eigH)))*conj(lambda);
%                     psi = diag( conj(eigH(:))./(abs(eigH(:)).^2 + noiseVar));
                    F_N = exp(-1i*2*pi*(0:this.NumSymbols-1).'*(0:this.NumSymbols-1)/this.NumSymbols);
                    F_M = exp(-1i*2*pi*(0:this.NumSubcarriers-1).'*(0:this.NumSubcarriers-1)/this.NumSubcarriers);
                    eq = kron(F_M,F_N)' * psi * kron(F_M,F_N);
                    x = eq*Ydd(:);
            end
            
            if this.OutputConditionNumber == false
                conditionNumber = nan;
            end
            
        end
        
        function [x] = equalizeUsingPilotResponse(this, Ydd, Hdd)
            % ----- DESCRIPTION -----
            % This function uses values of the delay Doppler domain channel matrix to generate OFDM-symbol-wise channel matrix (Hi).
            % Then, the received symbols are equalized by multiplying the inverse matrix of the generated channel matrix.
            % Note: This basically doesn't work in fractional Doppler channels. Due to this, only ZF is implemented.
            % ----- INPUT PARAMETERS -----
            %   - Ydd: received symbols in the Delay doppler domain
            %   - pathGains, pathDopplers, pathDelays, pathOffsets: path parameters
            %   - noiseVar: noise variance for MMSE

            M = this.NumSubcarriers;
            N = this.NumSymbols;
            Ncp = this.CyclicPrefixLength;
            Hdd = 1/sqrt(M*N)*Hdd;
 
            %% Re-generate OFDM-symbol-wise (i: symbol number) time domain channel matrix: Hi (size: M x M)
            % Initialize matrixes
            Hi = cell(N,1);
            
            for i = 1:N
                Hi{i} = zeros(M,M);
                for l = 0:M-1
                    for k = 0:N-1
                        timeIndices = (((i-1)*(M+Ncp)) : (i*(M+Ncp)-Ncp-1)) - l;
                        t = timeIndices / this.SamplingRate;
                        phaseShift = exp(1i * (2*pi * k * this.SamplingRate / ((M+Ncp)*N) * t)); 
                        dopplerMat = diag(phaseShift);
                        Hi{i} = Hi{i} + Hdd(l+1,k+1) * circshift(dopplerMat, -l, 2); 
                    end
                end
            end
            
            % Combine channel matrixes into a block-wise channel matrix: H (size: MN x MN)
            H = zeros(M*N, M*N);
            for i = 1:N
                H((i-1)*M+1:i*M,(i-1)*M+1:i*M) = Hi{i};
            end
            spH = sparse(H);
            
            % Convert H into the time domain channel matrix in the delay-Doppler domain: Hdd (a.k.a. phi in the document)
            spHdd = 1/N * this.spKronDftMat * spH * this.spKronIdftMat;
            HddEst = full(spHdd);
            x = HddEst\Ydd(:);  % Zero forcing
        end        
    end
end