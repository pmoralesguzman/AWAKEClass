%________________________________________________________________________
% Script to produce waterfall plots for different xi ranges of the proton bunch.
% Special version to produce the plot for the paper.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 18/06/2020
%________________________________________________________________________


datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'gm20'};
plasmaden = 1.81e14;
dump_list = 0:100;
useAvg = true;
dataformat = 'h5';
mean_weights_flag = true;

% load color order for 9 gradients
load('color_purple_to_green.mat');

OPA = OsirisPhaseAnalysis('datadir',datadirs{1},...
    'property','fields','wakefields_direction','long',...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,...
    'force_waterfall',false,...
    'mean_weights_flag',mean_weights_flag);

OPA.xi_list = OPA.plasma_wavelength:2*OPA.plasma_wavelength:20;


dephasing_all = cell(1,length(datadirs));
dephasing_std = cell(1,length(datadirs));
for d = 1:length(datadirs)
    
    datadir = datadirs{d};
    
    switch datadir
        case  'gm20'
            title_g = 'g = -20 %';
        case 'gm15'
            title_g = 'g = -15 %';
        case 'gm10'
            title_g = 'g = -10 %';
        case 'gm5'
            title_g = 'g = -5 %';
        case  'g0'
            title_g = 'g = 0 %';
        case 'gp5'
            title_g = 'g = +5 %';
        case 'gp10'
            title_g = 'g = +10 %';
        case 'gp15'
            title_g = 'g = +15 %';
        case 'gp20'
            title_g = 'g = +20 %';
    end
    
    OPA.datadir = datadir;
    
    OPA.dephasing_average();
    dephasing_all{d} = OPA.mean_dephasing;
    dephasing_std{d} = OPA.std_dephasing;
end
%% Plotting
dephasing_z = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),...
    length(OPA.mean_dephasing));

fig1 = figure(1);
% colororder(cc)
p1 = cell(1,length(datadirs));
for d = 1:length(datadirs)
    hold on
    fill([dephasing_z';flipud(dephasing_z')],...
        [dephasing_all{d}' - dephasing_std{d}'/2;flipud(dephasing_all{d}' + dephasing_std{d}'/2)],...
        cc(d,:),'linestyle','none','facealpha',0.3);
    p1{d} = plot(dephasing_z,...
        dephasing_all{d},'LineWidth',2,'color',cc(d,:));
    hold off
end

ylabel('dephasing (\lambda_p)')
xlabel('propagation distance (m)')
legend([p1{:}],'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20','Location','best')
xlim([0 10])
ylim([-5 1])
title('weighted average dephasing')
drawnow;

P = Plotty('fig_handle',fig1);
P.plots_dir = 'dephasing/gradsim';
P.plot_name = 'wmeanx';
P.save_plot();


