% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef OtfsPrecoder < matlab.System
    
    properties
        NumDelayBins;
        NumDopplerBins;
        NumTransmitAntennas;
        NumReceiveAntennas;
    end
    
    methods
        function this = OtfsPrecoder(varargin) % Constructor
            setProperties(this, nargin, varargin{:}, ...
                'NumDelayBins', ...
                'NumDopplerBins', ...
                'NumTransmitAntennas', ...
                'NumReceiveAntennas' ...
            );
        end
        
        function txFreqTime = encode(this, txDelayDoppler)
            % OTFS modulation
            txFreqTime = isfft2d(txDelayDoppler);
            txFreqTime = permute(txFreqTime, [2 1 3]);

            % Spatial precoder
            % -- No spatial precoder used for OTFS --
            Wn = 1/sqrt(this.NumTransmitAntennas);
            txFreqTime = Wn * txFreqTime;
        end
        
        function rxDelayDoppler = decode(this, rxFreqTime, varargin)
            % OTFS demodulation
            rxDelayDoppler = sfft2d(rxFreqTime);
            rxDelayDoppler = permute(rxDelayDoppler, [2 1 3]);
        end            
    end
end