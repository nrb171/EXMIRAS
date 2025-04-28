function [] = print2(varargin)
    % I hate the default matlab print function, so I made my own.
    % This function takes a figure and a filename and saves the figure to the
    % specified filename. The file format is determined by the file extension.
    % Supported file formats are .pdf, .eps, .svg, and .png.
    %
    % INPUTS:
    % fig: figure handle
    % filename: string or char array containing the filename
    %
    % OPTIONAL INPUTS:
    % Verbose: boolean, if true, prints a message to the console where the file
    %       was saved. Default is true.
    % FontName: string or char array containing the font name. Default is
    %       'Times New Roman'.
    % FontSize: numeric value containing the font size. Default is 6.
    %
    % DEPENDENCIES:
    % plot_whiteSpaceOptimizer.m
    %
    % EXAMPLES:
    % print2(fig, 'myFigure.pdf')
    % print2(fig, 'myFigure.eps')
    % print2(fig, 'myFigure.png', 'Verbose', false, 'FontName', 'Arial',
    % 'FontSize', 8)


    %parse the input
    if numel(varargin) >0
        if isa(varargin{1}, 'string') || isa(varargin{1}, 'char')
            varargin{2} = varargin{1};
            varargin{1} = gcf;
        end
    end
    p = inputParser;
    p.addOptional('fig', gcf, @(x) ishandle(x) && strcmp(get(x, 'Type'), 'figure'));
    p.addOptional('filename', '.temp.png',@(x) ischar(x) | isstring(x));
    p.addParameter('Quality', '-r100', @(x) ischar(x) | isstring(x));
    p.addParameter('Verbose', true, @islogical);
    p.addParameter('FontName', 'Times New Roman', @(x) ischar(x) | isstring(x));
    p.addParameter('FontSize', 6, @isnumeric);

    p.parse(varargin{:});
    fig = p.Results.fig;
    filename = p.Results.filename;
    fontSize = p.Results.FontSize;
    fontName = p.Results.FontName;
    quality = p.Results.Quality;
    plot_whiteSpaceOptimizer(fig, 'FontSize', fontSize, 'FontName',fontName);

    %set axis layer to top
    axs = findall(fig, 'type', 'axes');
    for i = 1:length(axs)
        axs(i).Layer = 'top';
    end


    %calculate the format of the file from the filename
    [pathstr, name, ext] = fileparts(filename);
    if ext == ".pdf"
        exportgraphics(fig,filename,'ContentType','vector')
    elseif ext == ".eps"
        exportgraphics(fig,filename,'ContentType','vector')
    elseif ext == ".svg"
        saveas(fig,filename)
    elseif ext == ".png"
        print(fig,'-dpng',quality,filename)
        % f=getframe(fig);
        % f=frame2im(f);
        % imwrite(f, filename);
    else
        error("Unsupported file format")
    end

    if p.Results.Verbose
        fprintf("Saved figure to %s\n", filename)
    end