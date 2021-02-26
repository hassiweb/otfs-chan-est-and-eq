% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

classdef Interleaver < matlab.System
    
    properties
        % Input properties
        InterleaverIndicesSource;
        InputLength;  % for random interleaver
        InterleaverRandomSeed;
        
        % Initialized properties
        InterleaverIndices;
    end
    
    methods
        function this = Interleaver(varargin)
            setProperties(this, nargin, varargin{:});
            
            switch this.InterleaverIndicesSource
                case "Random"
                    this.InterleaverIndices = randintrlv(1:this.InputLength, this.InterleaverRandomSeed);
            end
        end
        
        function out = interleave(this, in)
            out = intrlv(in, this.InterleaverIndices);
        end
        
        function out = deinterleave(this, in)
            out = deintrlv(in, this.InterleaverIndices);
        end
    end
end