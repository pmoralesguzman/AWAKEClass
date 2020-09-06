%________________________________________________________________________
% Plot the 2D wakefields together with the proton bunch, or each
% individually. Optimized to plot JohnMix sims.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 31/08/2020
%________________________________________________________________________

close all;

% file location variables
datadir = 'johnmixthin';
dataformat = 'h5';
useAvg = false;
dump_list = 2:1:2;

% saving data
save_flag = false;
% plots_dir = ['field_density/grads/',datadir];
save_format = {'png','eps'};

% plasma properties
plasmaden = 7e14;

% choose fields to plot
wakefields_direction = 'long'; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 10];
xi_range = [100 0];

% choose property to plot
property_plot = 'wakefields'; % density, wakefields, both

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose if make pause or not
make_pause = false;

% directory to save the plots
plots_dir = ['lineouterase/',datadir,'/',...
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
    'make_pause',make_pause,'lineout_point',10);

P.lineout_plot();
P.save_flag = true;
P.plot_name = 'r0p1kp';
title(['propagation dist. = ' num2str(P.propagation_distance/100),' m',' (r = 0.1/k_p)']); 
P.save_plot();
writematrix(P.lineout','r0p1kp.dat');
writematrix(P.z','z.dat');

