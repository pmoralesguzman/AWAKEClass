%________________________________________________________________________
% Plot lineout standard file
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 21/02/2021
%________________________________________________________________________


% data directory
datadir = ['g0']; % data directory
plots_dir = ['fft/rz/grads/eseed']; % save directory for the plot

% parameters
plasma_density = 1.81e14;

% properties
property_plot = 'wakefields';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump_list = 50:1:50;
dataformat = 'mat';
useAvg = false;

% plot parameters
lineout_point = 5;
xi_range = [15,0];


P = Plotty('datadir',datadir,'plots_dir',plots_dir,...
    'plasmaden',plasma_density,'property_plot',property_plot,...
    'species',species,'field',field,'direction',direction,...
    'lineout_point',lineout_point,'xi_range',xi_range,...
    'dump_list',dump_list,'dataformat',dataformat,'useAvg',useAvg);

P.plot_lineout();




