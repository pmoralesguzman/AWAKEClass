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
% Last update: 07/01/2021
%---------------------------------------------------------------------

function CreatePathClass(varargin)
% Create path where all studies and files are.

if nargin == 0
    force = false;
else
    force = true;
end

pathCell = regexp(path, pathsep, 'split');

paths = {'../AWAKE_Class',... all directories with analysis code
    'E:/AWAKE_Data',... for the external hard drive (must be set to E:/)
    'D:/AWAKE_Data_HD',... for the hard drive in the work laptop
    'C:/AWAKE_Data_laptop',... for the hard drive in the personal laptop
    '../simulations'... for the server
    % add yours
    };


for pp = 1:length(paths)
    s = what(paths{pp});
    if ispc && ~isempty(s) && ~any(strcmpi(s.path, pathCell)) || force  % Windows is not case-sensitive
        addpath(genpath(paths{pp}))
    elseif (~ispc) && ~isempty(s) && ~any(strcmp(s.path, pathCell)) || force 
        addpath(genpath(paths{pp}))
    end
end % paths

end

