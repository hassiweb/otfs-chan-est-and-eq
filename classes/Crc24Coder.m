% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef Crc24Coder
    
    properties
        Crc24Generator;
        Crc24Detector;
    end
    
    methods
        function this = Crc24Coder(polyType, vargin)
            switch polyType
                case 'A'
                    poly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
                case 'B'
                    poly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];
                otherwise
                    poly = vargin{1};
            end
            
            this.Crc24Generator = comm.CRCGenerator('Polynomial',poly);
            this.Crc24Detector = comm.CRCDetector('Polynomial',poly);
        end
        
        function out = generate(this, in)
            in(isnan(in)) = 0; % Dealt with NaN as 0 (zero)
            out = step(this.Crc24Generator, in);
        end
        
        function [out, error] = detect(this, in)
            [out, error] = step(this.Crc24Detector, in);
        end
    end
end