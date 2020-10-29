%________________________________________________________________________
% Plot the 2D wakefields together with the proton bunch, or each
% individually
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 05/08/2020
%________________________________________________________________________

close all;

% file location variables
datadir = 'gm20';
dataformat = 'h5';
useAvg = false;
dump_list = 102:1:110;

% saving data
save_flag = true;
% plots_dir = ['field_density/grads/',datadir];
save_format = {'png'};

% plasma properties
plasmaden = 7e14;

% choose fields to plot
wakefields_direction = 'long'; % trans, long

% choose species density to plot
species = 'electron_bunch'; 

% choose limits (in cm, must denormalize)
trans_range = [0 0.15];
% xi_range = [4.8 4.4];
xi_range = [6.5 0];

% choose property to plot
property_plot = 'wakefields'; % density, wakefields, both

% create movie or not
create_movie = true;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose scaling
plot_scale = 'linear';

% choose if make pause or not
make_pause = true;

% directory to save the plots
plots_dir = ['field_density/grads/',datadir,'/',...
    property_plot,'/',wakefields_direction,'/',...
    'xi',num2str(round(xi_range(1))),'xi',...
    num2str(round(xi_range(2)))];



P = Plotty('datadir',datadir,'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'fig_number',1);

P.field_density_plot();

% yline(0.868,'r','LineWidth',2)
% yline(-0.868,'r','LineWidth',2)
% ylim([-1.6,1.6])
P.save_plot();

