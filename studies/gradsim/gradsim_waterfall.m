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


datadirs = {'gp20'};
% datadirs = {'R2gap0_2e14'};
plasmaden = 2e14;
dump_list = 0:1:200;
useAvg = true;
dataformat = 'h5';
% xlims = {[18 20],[8 10],[0 2]};
title_position = {'back','midle','front'};


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
        otherwise
            title_g = '';
            
    end
    
    P = Plotty('datadir',datadir,...
        'wakefields_direction','long',...
        'plasmaden',plasmaden,...
        'dump_list',dump_list,'useAvg',useAvg,...
        'dataformat',dataformat,'do_plot',true);
    OPA = OsirisPhaseAnalysis('datadir',datadir,...
        'wakefields_direction','long',...
        'plasmaden',plasmaden,...
        'dump_list',dump_list,'useAvg',useAvg,...
        'dataformat',dataformat);
    
    OPA.build_waterfall();
    P.waterfall_z = OPA.waterfall_z;
    P.waterfall_plot();
    
    %% plot three xi ranges
    fig = figure(100);
    fig.Units = 'normalized';
    fig.Position = [0.0955, 0.1322, 0.8095, 0.45];
    
    subp = cell(1,1);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Normal');
    % t.InnerPosition = [0.0955, 0.1322, 0.8095, 0.7916];
    for s = 1:1
        nexttile
        imagesc(P.waterfall_xi,P.waterfall_z,rot90(P.waterfall_mat,2));
        ax = gca;
        ax.XDir = 'reverse';
        ax.YDir = 'normal';
%         xlim(xlims{s})
        ylim(P.waterfall_z)
        
        title(title_position{s})
    end
    
    title(t,title_g)
    xlabel(t,'\xi (cm)')
    ylabel(t,'propagation distance (m)')
    colormap(bluewhitered);
    cbar = colorbar;
    cbar.Label.String = 'E_z (MV/m)';
    P.plots_dir = 'waterfall/gradsim';
    P.plot_name = datadir;
%     P.save_plot(fig);
end
