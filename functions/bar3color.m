% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

function bar3color(gridSymbols, varargin)
    barplot = bar3((gridSymbols));

    % Colorize according to Z-axis
    for k = 1:length(barplot)
        zdata = barplot(k).ZData;
        barplot(k).CData = zdata;
        barplot(k).FaceColor = 'interp';
    end
    
    % Display labels if inputted
    if nargin >= 3
        ylabel(varargin{1})
        xlabel(varargin{2})

        % Display labels if inputted
        if nargin >= 5
            xticklabels(varargin{3})
            yticks(1:length(varargin{4}))
            yticklabels(varargin{4})
        end
    end
end

