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
dump        = 50;
dataformat  = 'mat';
property    = 'density';
    field   = 'e';
    species = 'proton_beam';
raw_dataset = 'x';
direction   = 'z';


OD = OsirisDenormalizer('datadir',datadir,'plasmaden',plasmaden,...
    'dump',dump,'dataformat',dataformat,'property',property,...
    'field',field,'species',species,'direction',direction);
OD.getdata();
OD.assign_density();
OD.denorm_distance();

OD.propagation_distance