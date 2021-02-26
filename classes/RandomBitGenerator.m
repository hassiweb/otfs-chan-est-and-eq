% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef RandomBitGenerator < matlab.System
    
    properties 
        DataLength;
        DataBits;
    end
    
    methods
        function this = RandomBitGenerator(varargin)
            setProperties(this, nargin, varargin{:}, 'DataLength');
        end
        
        function dataBits = generate(this)
            dataBits = randi([0 1], this.DataLength, 1);
            this.DataBits = dataBits;
        end
        
        function [numErrors, ber] = detect(this, hardIn)
            numErrors = sum(abs(this.DataBits-hardIn));
            ber = numErrors / this.DataLength;
        end

        function [numErrors, ber] = detectSoftBits(this, softIn)
            numErrors = sum(abs(this.DataBits-this.soft2hard(softIn)));
            ber = numErrors / this.DataLength;
        end
        
        function hard = soft2hard(~, soft)
            hard = zeros(size(soft));
            hard(soft<0) = 1; 
            hard(soft>=0) = 0; 
        end

    end
end