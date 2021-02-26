% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef FractionalDelayFilter < matlab.System
    properties
        FilterOrder; % With a high filter order (>=10 (as a rule of thumb)), the accuracy may degrade due to the quantization error.
        
        fdfCoeffs;
        nonZeroSampIndices;
    end
        
    methods
        function this = FractionalDelayFilter(varargin) % Constructor
            if nargin > 0
                setProperties(this, nargin, varargin{:}, 'FilterOrder');
            end
            
            if isempty(this.FilterOrder)  % The filter order should be less than 10 because the accuracy will be degraded due to the rounding error even in float64
                this.FilterOrder = 8;
            end
        end
        
        function init(this, delays)
            this.fdfCoeffs = cell(length(delays),1);
            this.nonZeroSampIndices = zeros(length(delays),1);
            
            for delayIndex = 1:length(delays)
                delay = delays(delayIndex);
                filterOrder = min(floor(2*delay+1), this.FilterOrder);  
                % Note: The filter must be causal. The filter order is selected based on the delay to be causal and the parameter.

                % Calc a modified coefficient matrix of modified Farrow structure
                % Note: This algorithm is based on a book below.
                % [1] Välimäki, V., and T. I. Laakso. “Fractional Delay Filters—Design and Applications.? In Nonuniform Sampling, edited by Farokh Marvasti, 835?95. Information Technology: Transmission, Processing, and Storage. Boston, MA: Springer US, 2001. https://doi.org/10.1007/978-1-4615-1229-5_20.
                N = filterOrder;
                U = fliplr(vander(0:N));
                Q = invcramer(U);  % solve inverse matrix using Cramer's rule
                T = zeros(N+1,N+1);
                for n = 0:N
                    for m = 0:N
                        if n >= m
                            T(m+1,n+1) = round(N/2)^(n-m) * nchoosek(n,m);
                        end
                    end
                end
                Qf = T*Q;
                
                % Calculate filter coefficients of the FIR filter
                fracDelay = mod(delay,1);
                delayVec = fracDelay.^(0:filterOrder);
                this.fdfCoeffs{delayIndex} = delayVec * Qf;

                % Calculate non zero sample index (=start index)
                if mod(filterOrder,2)==0  % even
                    this.nonZeroSampIndices(delayIndex) = round(delay+0.5) - filterOrder/2;  % non-zero sample index, written as equation (3.37) in [1]
                else  % odd
                    this.nonZeroSampIndices(delayIndex) = floor(delay) - (filterOrder-1)/2;
                end
            end
            
            % Normalize the filter coefficient when a high filter order (since the accuracy of the filter is degraded)
            if this.FilterOrder > 9
                sinewave = this.apply(sin(0:pi/40:4*pi));
                calib = rms(sinewave(this.nonZeroSampIndices(end)+maxFilterOrder+20:this.nonZeroSampIndices(end)+maxFilterOrder+100-1,:))' / sqrt(1/2);
                for delayIndex = 1:length(delays)
                    this.fdfCoeffs{delayIndex} = this.fdfCoeffs{delayIndex}/calib(delayIndex);
                end
            end
        end
        
        
        function out = apply(this, in)
            numDelays = length(this.nonZeroSampIndices);
            tailLength = this.nonZeroSampIndices(numDelays) + length(this.fdfCoeffs{numDelays}) - 1;
            out = zeros(length(in)+tailLength, numDelays);
            for delayIndex = 1:numDelays
                filterout = conv(this.fdfCoeffs{delayIndex},in);
                if this.nonZeroSampIndices(delayIndex) > 0
                    out(this.nonZeroSampIndices(delayIndex)+1:this.nonZeroSampIndices(delayIndex)+length(in)+length(this.fdfCoeffs{delayIndex})-2, delayIndex) = filterout(2:end);
                else
                    out(1:length(in)+length(this.fdfCoeffs{delayIndex})+this.nonZeroSampIndices(delayIndex)-2, delayIndex) = filterout(-this.nonZeroSampIndices(delayIndex)+2:end);
                end
            end
        end
    end
end