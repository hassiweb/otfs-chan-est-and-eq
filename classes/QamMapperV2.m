% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef QamMapperV2 < matlab.System
    properties
        ModOrder; % number of points in the signal constellation (e.g. 4, 16, 64, 256)
        InputType;  % 'integer' or 'bit': followed by MATLAB function 'qammod'
        OutputType;  % 'integer', 'bit', 'llr', 'approxllr': followed by MATLAB function 'qamdemod'
        UnitAveragePower;
        LLROverflowPrevention;
    end
        
    methods
        function this = QamMapperV2(varargin) % Constructor
            if nargin > 0
                setProperties(this, nargin, varargin{:}, ...
                    'ModOrder', ...
                    'InputType', ...
                    'OutputType', ...
                    'UnitAveragePower', ...
                    'LLROverflowPrevention');
            end
        end
        
        function out = map(this, in)
            out = qammod(double(in), this.ModOrder, getconstellation(this.ModOrder), ...
                'InputType', this.InputType, ...
                'UnitAveragePower', this.UnitAveragePower);
        end

        function out = demap(this, in, varargin)
            switch this.OutputType
                case {'bit', 'integer'}
                    out = qamdemod(in, this.ModOrder, getconstellation(this.ModOrder), ...
                        'OutputType', this.OutputType, ...
                        'UnitAveragePower', this.UnitAveragePower, ...
                        'PlotConstellation', false);  % This doesn't work well
                case {'llr', 'approxllr'}
                    noiseVariance = varargin{1};
                    out = qamdemod(in, this.ModOrder, getconstellation(this.ModOrder), ...
                        'OutputType', this.OutputType, ...
                        'UnitAveragePower', this.UnitAveragePower, ...
                        'NoiseVariance', noiseVariance, ...
                        'PlotConstellation', false);  % This doesn't work well

                    if this.LLROverflowPrevention
                        while 1
                            if sum(isnan(out)) > 0
                                noiseVariance = noiseVariance * 10;  % update noise variance for NaN values
                            else
                                break;
                            end

                            results = qamdemod(in, this.ModOrder, getconstellation(this.ModOrder), ...
                                'OutputType', 'llr', ...
                                'UnitAveragePower', true, ...
                                'NoiseVariance', noiseVariance, ...
                                'PlotConstellation', false);  % THis doesn't work well

                            out(isnan(out)) = results(isnan(out)); % Update values that are NaN
                        end

                    else  % Give warning when LLR is overflowed
                        if sum(isnan(out)) > 0
                            out(isnan(out)) = 0;
                            warning("[WARNING] The output of QAM demmaping comprises NaN value(s). NaN values is converted into '0' to prevent the error in the later processes,");
                        end
                    end
                    
            end
        end
    end
end
        