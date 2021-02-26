% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

function [fdfedPathGains, fdfedPathDopplers, fdfedPathDelays, fdfedPathOffsets] = fdfpathpars(pathGains, pathDopplers, pathFractionalDelays, pathOffsets, samplingRate)
% ### DESCRIPTION ###
% This function converts path parameters with fractional delays to path parameters with integer 
% delays using the fractional delay filter. 
% 
% ### INPUT PARAMETERS ###
%   - pathGains: path gains
%   - pathDopplers: path Doppler shifts 
%   - pathFractionalDelays: path delays with a fractional delay
%   - pathOffsets: initial offset of each path
%   - varargin: if there is an another parameter

FDF = FractionalDelayFilter;
FDF.init(pathFractionalDelays*samplingRate);

totalNumPars = 0;
for p = 1:length(pathFractionalDelays)
    totalNumPars = totalNumPars + length(FDF.fdfCoeffs{p}) - 1;
end

fdfedPathGains = zeros(totalNumPars,1);
fdfedPathDopplers = zeros(totalNumPars,1);
fdfedPathDelays = zeros(totalNumPars,1);
fdfedPathOffsets = zeros(totalNumPars,1);

newPathIdx = 1;
for p = 1:length(pathFractionalDelays)
    for d = 1:length(FDF.fdfCoeffs{p})
        if d - 2 + FDF.nonZeroSampIndices(p) < 0 % ignore non causal delay index
            continue;
        end
        fdfedPathGains(newPathIdx) = pathGains(p) * FDF.fdfCoeffs{p}(d);
        fdfedPathDopplers(newPathIdx) = pathDopplers(p);
        fdfedPathDelays(newPathIdx) = (d - 2 + FDF.nonZeroSampIndices(p))/samplingRate;
        fdfedPathOffsets(newPathIdx) = angle( exp( 1i * ( 2*pi*pathDopplers(p)*fdfedPathDelays(newPathIdx) + pathOffsets(p)) ) );
        
        newPathIdx = newPathIdx + 1;
    end
end

% Filter very small path
thresh = 1e-5;
fdfedPathDopplers(abs(fdfedPathGains) < thresh) = [];
fdfedPathDelays(abs(fdfedPathGains) < thresh) = [];
fdfedPathOffsets(abs(fdfedPathGains) < thresh) = [];
fdfedPathGains(abs(fdfedPathGains) < thresh) = [];

end