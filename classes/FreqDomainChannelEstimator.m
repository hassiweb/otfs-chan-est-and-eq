% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef FreqDomainChannelEstimator < matlab.System
    properties
        EstimationMethod;
        InputType;
        CPOFDMModulator;
        FFTSize;
        CyclicPrefix;
        NumGuardBandCarriers;
    end
    
    methods
        function this = FreqDomainChannelEstimator(varargin)
            setProperties(this, nargin, varargin{:}, ...
                'EstimationMethod', ...
                'InputType', ...
                'CPOFDMModulator', ...
                'FFTSize', ...
                'CyclicPrefix');
        end
        
        function H = estimate(this, varargin)
            switch this.EstimationMethod
                case 'Ideal'
                    switch this.InputType
                        case 'TimeDomainChannelImpulseResponse'
                            cir = varargin{1};
                            %Calculate frequency response
                            H = this.CPOFDMModulator.demodulate(cir); H = H(:);
                    end
            end
        end
    end
end