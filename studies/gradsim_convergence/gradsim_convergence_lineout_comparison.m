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
datadirs = {'g0z2','g0r2','g0'};
% datadirs = {'g0zh','g0rh','g0'};

dataformat = 'mat';
useAvg = false;
dump_list = 100:1:100;

% saving data
save_flag = false;
save_format = {'png','eps','fig'};

% plasma properties
plasmaden = 1.81e14;

% choose fields to plot
wakefields_direction = 'long'; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.0536];
xi_range = [21 0];
lineout_point = 5;

% choose property to plot
property_plot = 'density'; % density, wakefields, both

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose if make pause or not
make_pause = false;



% directory to save the plots
plots_dir = ['gradsim_convergence/','lineoutcompare/',datadirs{1},'/',...
    property_plot,'/',wakefields_direction,'/',...
    'xi',num2str(round(xi_range(1))),'xi',...
    num2str(round(xi_range(2)))];

for n = 1:length(dump_list)

P = Plotty('datadir',datadirs{1},'dataformat','h5',...
    'useAvg',useAvg,'dump_list',dump_list(n),...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'lineout_point',lineout_point);

P.lineout_plot();

P = Plotty('datadir',datadirs{2},'dataformat','h5',...
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
if length(datadirs) == 3
P = Plotty('datadir',datadirs{3},'dataformat',dataformat,...
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
end

legend(datadirs{:},'location','best')
P.save_flag = true;
P.plot_name = ['lineout_comp','n',num2str(n)];
P.fig_handle = gcf;
ylabel('sum of density on axis')
title(['propagation dist. = ' num2str(P.propagation_distance/100),' m',' (r = [0,0.536] mm)']); 
P.save_plot();
end

