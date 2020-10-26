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


datadirs = {'s0','s1','s2'};
plasmaden = 1.81e14;
dump_list = 0:100;
useAvg = false;
dataformat = 'h5';
xlims = {[14 21],[0 7]};
title_position = {'back','front'};


for d = 1:length(datadirs)
    
    datadir = datadirs{d};
    
    switch datadir
        case  's0'
            title_s = 'seed pos. = center';
        case 's1'
            title_s = 'seed pos. = center - 1 \sigma';
        case 's2'
            title_s = 'seed pos. = center - 2 \sigma';
    end
    
    P = Plotty('datadir',datadir,...
        'property','fields','wakefields_direction','long',...
        'plasmaden',plasmaden,...
        'dump_list',dump_list,'useAvg',useAvg,...
        'dataformat',dataformat,'do_plot',false);
    P.waterfall_plot();
    
    %% plot three xi ranges
    fig = figure(100);
    fig.Units = 'normalized';
    fig.Position = [0.0955, 0.1322, 0.8095, 0.45];
    
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Normal');
    for s = 1:2
        nexttile
        imagesc(P.waterfall_xi,P.waterfall_z,rot90(P.waterfall_build,2));
        ax = gca;
        ax.XDir = 'reverse';
        ax.YDir = 'normal';
        xlim(xlims{s})
        ylim(P.waterfall_z)
        
        title(title_position{s})
    end
    
    title(t,title_s)
    xlabel(t,'\xi (cm)')
    ylabel(t,'propagation distance (m)')
    colormap(bluewhitered);
    cbar = colorbar;
    cbar.Label.String = 'E_z (MV/m)';
    P.plots_dir = 'waterfall/seedpos';
    P.plot_name = datadir;
    P.save_plot(fig);
end
