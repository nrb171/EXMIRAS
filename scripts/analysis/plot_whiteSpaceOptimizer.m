function [figIn, axIn] = plot_WhiteSpaceOptimizer(figIn, varargin)
    % This function will take an existing axes and modify it to minimize the
    % white space around the plot.
    %
    % INPUTS:
    %   figIn: figure handle
    % Optional Inputs:
    %   FontSize: font size for all text in the figure
    %   FontName: font name for all text in the figure
    %
    % OUTPUTS:
    %   figIn: figure handle - modified
    %   axIn: axes handle - modified
    %
    % EXAMPLES:
    %   figIn = figure;
    %   axIn = axes;
    %   [figIn, axIn] = plot_WhiteSpaceOptimizer(figIn, 'FontSize', 12,
    %   'FontName', 'Arial')
    %
    % CHANGELOG:
    %   2023/08 - Initial function write - Nicholas Barron
    %


    % Parse the input
    p = inputParser;
    p.addRequired('figIn');

    %p.addOptional('axIn') ;
    p.addParameter('FontSize', 12);
    p.addParameter('FontName', 'Arial');

    p.parse(figIn, varargin{:});
    fontSize = p.Results.FontSize;
    fontName = p.Results.FontName;
    %figIn = p.Results.figIn;
    axIn = findall(figIn,'type','axes');




    %update font for figure
    set(findall(figIn,'type','text'),'FontSize',fontSize)
    set(findall(figIn,'type','text'),'FontName',fontName)

    %update font for axes
    for ii = 1:numel(axIn)
        set(axIn(ii), 'FontSize', fontSize)
        set(axIn(ii), 'FontName', fontName)
        set(axIn(ii), 'LooseInset', max(get(axIn(ii), 'TightInset'), 0.025))
    end
end