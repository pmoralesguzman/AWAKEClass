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
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'g0'};

plasmaden = 1.81e14;
property = 'fields';
dump_list = 0:100;
useAvg = true;
dataformat = 'h5';

trans_range = [0 0.02]; % if density
% dephasing_xis = [1,3,6]; % cm
dephasing_search = '0x'; % 0x, max
force_waterfall = false;

% initialize some plotty
P1 = Plotty('plots_dir','dephasing/gradsim/wakefields');
P2 = Plotty('plots_dir','dephasing/gradsim/wakefields');

cc = [118,42,131
153,112,171
194,165,207
231,212,232
10,10,10
217,240,211
166,219,160
90,174,97
27,120,55]/256;

OPA = OsirisPhaseAnalysis('datadir','g0',...
    'property',property,'species','proton_beam',...
    'direction','z',...
    'wakefields_direction','long',...
    'trans_range',trans_range,...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,...
    'lineout_point',8);

dephasing_xis = OPA.plasma_wavelength:2*OPA.plasma_wavelength:10;

for ph = 1:length(dephasing_xis)
    dephasing_xi = dephasing_xis(ph);
    %     close all;
    for d = 1:length(datadirs)
        
        datadir = datadirs{d};
        
        OPA = OsirisPhaseAnalysis('datadir',datadir,...
            'property',property,'species','proton_beam',...
            'direction','z',...
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
        if ph == 1
            % waterfall plot datas
            P1.waterfall_xi = OPA.waterfall_xi;
            P1.waterfall_z = OPA.waterfall_z;
            P1.waterfall_mat = OPA.waterfall_mat;
            P1.property = OPA.property;
            
            P1.waterfall_plot(10+d);
            x_range = [0 10];
            xlim(x_range);
            ind_x = fliplr((OPA.waterfall_xi < x_range(2)) & (OPA.waterfall_xi > x_range(1)));
            colormap_low = min(OPA.waterfall_mat(:,ind_x),[],'all');
            colormap_upp = max(OPA.waterfall_mat(:,ind_x),[],'all');
            caxis([colormap_low colormap_upp]);
            colormap(bluewhitered);
            ax = imgca(P1.fig_handle);
        end
        
        % add the dephasing line (zero-crossing or max) to the waterfall plot
        hold on
        plot((-OPA.dephasing_line)*OPA.plasma_wavelength+OPA.simulation_window-OPA.dephasing_first,...
            linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
            'LineWidth',1,'Parent',ax,'color','k')
        ax.XDir = 'reverse';
        ax.YDir = 'normal';
        title(['dephasing at different xi for mean'])
        hold off
        drawnow;
        P1.plot_name = ['gradsim_meandephasing_',datadir]; % PLOT NAME
        
        
    end
    P1.save_plot();
end



