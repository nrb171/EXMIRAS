%% startup file, initializes path variable/variables.
%  This file is executed at the start of MATLAB session.

% Add the path of the current directory
%recursively add all subdirectories to path
addpath(genpath('.'))
addpath(genpath('/h/eol/nbarron/workshop/lrose-emerald/matlab'))
addpath(genpath('/h/eol/nbarron/workshop/apar-scripts'))
addpath(genpath('/h/eol/nbarron/workshop/EXMIRAS/scripts'))
% other startup functions
removeTemp()


%fprintf('startup loaded :) \n')

function removeTemp()
    % Remove temporary files
    delete('.temp*')
end