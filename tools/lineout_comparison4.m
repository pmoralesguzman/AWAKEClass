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
datadirs = {'gm20','gm20d2','gm20d3'};
dataformat = 'h5';
useAvg = true;
dump_list = 1:1:100;

% saving data
save_flag = false;
% plots_dir = ['field_density/grads/',datadir];
save_format = {'png','eps','fog'};

% plasma properties
plasmaden = 1.81e14;

% choose fields to plot
wakefields_direction = 'long'; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.3];
xi_range = [21 0];
lineout_point = 5;

% choose property to plot
property_plot = 'wakefields'; % density, wakefields, both

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose if make pause or not
make_pause = false;



% directory to save the plots
plots_dir = ['lineouterase/',datadirs{1},'/',...
    property_plot,'/',wakefields_direction,'/',...
    'xi',num2str(round(xi_range(1))),'xi',...
    num2str(round(xi_range(2)))];

for n = 1:length(dump_list)

P = Plotty('datadir',datadirs{1},'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list(n),...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'lineout_point',lineout_point);

P.lineout_plot();

P = Plotty('datadir',datadirs{2},'dataformat',dataformat,...
    'useAvg',useAvg,'dump_list',dump_list(n),...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'lineout_point',lineout_point);


hold on
P.lineout_plot();
hold off
legend(datadirs{:})
P.save_flag = true;
P.plot_name = ['lineout_comp_gm20','n',num2str(n)];
P.fig_handle = gcf;
title(['propagation dist. = ' num2str(P.propagation_distance/100),' m',' (r = 0.05/k_p)']); 
P.save_plot();
end

