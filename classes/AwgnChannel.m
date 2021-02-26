% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef AwgnChannel
    properties
        inputType;
        AwgnObj;
    end
        
    methods
        function this = AwgnChannel(inputType) % Constructor
            this.inputType = inputType;
            switch inputType
                case 'variance'
                    this.AwgnObj = comm.AWGNChannel( ...
                        'NoiseMethod', 'Variance', ...
                        'VarianceSource', 'Input port', ...
                        'RandomStream', 'mt19937ar with seed', ...
                        'Seed', 1);
            end
        end
        
        function [out, noise] = add(this, in, varargin)
            switch this.inputType
                case 'variance'
                    noiseVar = varargin{1};
                    out = this.AwgnObj.step(in, noiseVar);
                    noise = out - in;
                case 'N0'
                    N0 = varargin{1};
                    noise = sqrt(N0/2) * complex(randn(size(in)), randn(size(in)));
                    out = in + noise;
            end
        end
    end
end
        