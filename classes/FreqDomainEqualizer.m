% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef FreqDomainEqualizer < matlab.System
    properties
        EqualizationAlgorithm;
    end
    
    methods
        function this = FreqDomainEqualizer(varargin)
            setProperties(this, nargin, varargin{:}, 'EqualizationAlgorithm');
        end
        
        function [X, noiseVarianceAfterEqualizer] = equalize(this, Y, H, noiseVariance)
            switch this.EqualizationAlgorithm
                case 'ZF'
                    eq = conj(H)./(conj(H).*H);
                case 'MMSE' 
                    eq = conj(H)./(noiseVariance + conj(H).*H);
            end

            % Apply equalizer and calculate enhanced noise variance
            X = Y.*eq;
            noiseVarianceAfterEqualizer  = (sqrt(noiseVariance) * abs(eq)).^2;
        end
    end
end