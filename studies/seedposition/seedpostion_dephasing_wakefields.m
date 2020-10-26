%________________________________________________________________________
% Script to calculate, from the waterfall data, the dephasing of the fields
% close to some selected xi.
% Especially developed for the seeding position studies.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/07/2020
%________________________________________________________________________

clear;
close all;
datadirs = {'s0','s1','s2'};
% datadirs = {'s2'};

plasmaden = 1.81e14;
property = 'fields';
dump_list = 0:100;
useAvg = false;
dataformat = 'h5';

trans_range = [0 0.02]; % if density
dephasing_xis = [1,7,14,20]; % cm
dephasing_search = '0x'; % 0x, max
force_waterfall = false;

% initialize some plotty
P1 = Plotty('plots_dir','dephasing/seedposition/');
P2 = Plotty('plots_dir','dephasing/seedposition/');

for ph = 1:length(dephasing_xis)
    dephasing_xi = dephasing_xis(ph);
    close all;
    for d = 1:length(datadirs)
        
        datadir = datadirs{d};
        
        
        OPA = OsirisPhaseAnalysis('datadir',datadir,...
            'property',property,'species','proton_beam',...
            'wakefields_direction','long',...
            'trans_range',trans_range,...
            'plasmaden',plasmaden,...
            'dump_list',dump_list,'useAvg',useAvg,...
            'dataformat',dataformat,...
            'dephasing_xi',dephasing_xi,...
            'dephasing_search',dephasing_search,...
            'force_waterfall',force_waterfall,...
            'lineout_point',8);
        OPA.dephasing();
        
        % waterfall plot datas
        P1.waterfall_xi = OPA.waterfall_xi;
        P1.waterfall_z = OPA.waterfall_z;
        P1.waterfall_mat = OPA.waterfall_mat;
        P1.property = OPA.property;
        
        P1.waterfall_plot(10+d);
        x_range = dephasing_xi+[-1 1];
        xlim(x_range);
        ind_x = fliplr((OPA.waterfall_xi < x_range(2)) & (OPA.waterfall_xi > x_range(1)));
        colormap_low = min(OPA.waterfall_mat(:,ind_x),[],'all');
        colormap_upp = max(OPA.waterfall_mat(:,ind_x),[],'all');
        caxis([colormap_low colormap_upp]);
        ax = imgca(P1.fig_handle);
        
        % load line from the microbunches for double plotting
        phase_filename = ['save_files/phase_',datadir,'_',num2str(dephasing_xi),'.mat'];
        bunches = load(phase_filename);
        bunches_x = bunches.phase_x_plot;
        bunches_y = bunches.phase_y_plot;
        
        % add the dephasing line (zero-crossing or max) to the waterfall plot
        hold on
        plot((-OPA.dephasing_line)*OPA.plasma_wavelength+OPA.simulation_window-OPA.dephasing_first,...
            linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
            'LineWidth',2,'Parent',ax,'color','k')
        plot(bunches_x,bunches_y,'LineWidth',2,'Parent',ax,'color',[1, 0, 0])
        ax.XDir = 'reverse';
        ax.YDir = 'normal';
        title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos. (',datadir,')'])
        hold off
        drawnow;
        P1.plot_name = ['seedpos_field_',datadir,'_',num2str(dephasing_xi)]; % PLOT NAME
        P1.save_plot();
        
        % just the dephasing line, holds the figure for the next datadir
        fig1 = figure(1);
        hold on
        p1 = plot(linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
            OPA.dephasing_line,'LineWidth',2);
        hold off
        
    end
end

% labels for the dephasing lines
ylabel('phase (\lambda_p)')
xlabel('propagation distance (m)')
legend('seeding center','seeding 1 \sigma_r ahead of center','seeding 2 \sigma_r ahead of center','Location','best')
xlim([0 10])
ylim([-1.6 0.1]);
title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos.'])

P2.plot_name = ['seedpos_dephasing_fields_',num2str(dephasing_xi)]; % PLOT NAME
P2.fig_handle = fig1;
P2.save_plot();


