%__________________________________________________________________________
% Example for the classes in the AWAKE Osiris Analysis Matlab Package
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 04/06/2020
%__________________________________________________________________________


datadir     = 'g0';
plasmaden   = 1.81e14;
dump        = 97;
dataformat  = 'mat';
raw_dataset = 'tag';
direction   = 'z';


OD = OsirisDenormalizer('datadir',datadir,'plasmaden',plasmaden,...
    'dump',dump,'dataformat',dataformat,'property','raw',...
    'raw_dataset',raw_dataset,'direction',direction);
OD.getdata();
OD.assign_raw();
length(OD.tag)


