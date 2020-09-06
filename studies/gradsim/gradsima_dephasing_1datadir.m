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
datadir = 'gp20';

plasmaden = 1.81e14;
property = 'fields';
dump_list = 10:1:100;
useAvg = false;
dataformat = 'mat';

wakefields_direction = 'trans';
lineout_point = 50;

trans_range = [0 0.02]; % if density
dephasing_xis = [1,2,3,4,4.9]; % cm
label_leg = num2cell(dephasing_xis);
dephasing_search = '0x'; % 0x, max
force_waterfall = false;

% initialize some plotty
P = Plotty('plots_dir','dephasing/gradsim/wakefields','save_format',{'png','eps'},...
    'wakefields_direction',wakefields_direction);

cc = [237,248,177
    199,233,180
    127,205,187
    65,182,196
    29,145,192
    34,94,168
    37,52,148
    8,29,88
    6,24,77
    4,19,66
    2,10,44]/256;

for ph = 1:length(dephasing_xis)
    dephasing_xi = dephasing_xis(ph);
    
    OPA = OsirisPhaseAnalysis('datadir',datadir,...
        'property',property,'species','proton_beam',...
        'direction','r',...
        'wakefields_direction',wakefields_direction,...
        'trans_range',trans_range,...
        'plasmaden',plasmaden,...
        'dump_list',dump_list,'useAvg',useAvg,...
        'dataformat',dataformat,...
        'dephasing_xi',dephasing_xi,...
        'dephasing_search',dephasing_search,...
        'force_waterfall',force_waterfall,...
        'lineout_point',lineout_point);
    OPA.dephasing();
    
    % waterfall plot datas
    P.waterfall_xi = OPA.waterfall_xi;
    P.waterfall_z = OPA.waterfall_z;
    P.waterfall_mat = OPA.waterfall_mat;   
    P.property = OPA.property;
    
            P.waterfall_plot(10+ph);
            x_range = dephasing_xi+[-1 1];
            xlim(x_range);
            ind_x = fliplr((OPA.waterfall_xi < x_range(2)) & (OPA.waterfall_xi > x_range(1)));
            colormap_low = min(OPA.waterfall_mat(:,ind_x),[],'all');
            colormap_upp = max(OPA.waterfall_mat(:,ind_x),[],'all');
            caxis([colormap_low colormap_upp]);
            colormap(bluewhitered);
            ax = imgca(P.fig_handle);
    
    % load line from the microbunches for double plotting
    %         phase_filename = ['save_files/phases/phase_',datadir,'_',num2str(dephasing_xi),'.mat'];
    %         bunches = load(phase_filename);
    %         bunches_x = bunches.phase_x_plot;
    %         bunches_y = bunches.phase_y_plot;
    
    % add the dephasing line (zero-crossing or max) to the waterfall plot
            hold on
            plot((-OPA.dephasing_line)*OPA.plasma_wavelength+OPA.simulation_window-OPA.dephasing_first,...
                linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
                'LineWidth',2,'Parent',ax,'color','k')
%             plot(bunches_x,bunches_y,'LineWidth',2,'Parent',ax,'color',[1, 0, 0])
            ax.XDir = 'reverse';
            ax.YDir = 'normal';
            title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos. (',datadir,')'])
            hold off
            drawnow;
            P.plot_name = ['gradsim_field_',datadir,'_',num2str(dephasing_xi)]; % PLOT NAME
            P.save_plot();
    
    % just the dephasing line, holds the figure for the next datadir
    fig1 = figure(1);
    colororder(cc);
    hold on
    p1 = plot(linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
        OPA.dephasing_line,'LineWidth',2);
    hold off
    
    % labels for the dephasing lines
    ylabel('phase (\lambda_p)')
    xlabel('propagation distance (m)')
    legend('1','2','3','4','5','Location','best')
    xlim([0 10])
    % ylim([-6 1]);
    title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos.'])
    
end

P.plot_name = ['gradsim_dephasing_fields_',datadir,num2str(dephasing_xi)]; % PLOT NAME
P.fig_handle = fig1;
P.save_plot();

