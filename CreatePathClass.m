%---------------------------------------------------------------------
% Tells MATLAB where the directories are.
% Run before any other script.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
%
% P. I. Morales Guzman
% Last update: 02/06/2020
%---------------------------------------------------------------------

% Create path where all studies and files are.

pathCell = regexp(path, pathsep, 'split');

paths = {'../AWAKE_Class',... all directories with analysis code
    'E:/AWAKE_Data/',... for the external hard drive (must be set to E:/)
    'D:/AWAKE_Data_HD',... for the hard drive in the work laptop
    'C:/AWAKE_Data_laptop',... for the hard drive in the personal laptop
    '../simulations'... for the server
    % add yours
    };

for pp = 1:length(paths)
    if ispc && ~any(strcmpi(paths{pp}, pathCell))  % Windows is not case-sensitive
        addpath(genpath(paths{pp}))
    elseif (~ispc) && ~any(strcmp(paths{pp}, pathCell))
        addpath(genpath(paths{pp}))
    end
end % paths

