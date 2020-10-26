%________________________________________________________________________
% FFT of the proton or plasma electrons density distribution, or of a
% lineout of the wakefields.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 08/06/2020
%________________________________________________________________________
clear;

% data directory
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'gm20d2','gm15d2','gm10d2','gm5d2','g0d2','gp5d2','gp10d2','gp15d2','gp20d2'};
leg         = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};

% datadirs = {'g0','gm20'};
% simulation parameters
dump = 133;
dataformat = 'mat';
useAvg = false;
use_raw = false;

% properties
plasma_density = 1.81e14;
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';

% limits
xi_range = [18 0.0]; % cm
trans_lims = [0 0.0536]; % cm
% analysis parameters
scan_type = 'cumulative'; % slice, cumulative
on_axis = 'sum'; % int, sum, intw, lineout

% switches
saveplots = true;

AFFT = AwakeFFT(...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'field',field,'direction',direction,...
    'dump',dump,'dataformat',dataformat,'useAvg',useAvg,...
    'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis);

fig_long = figure;
fig_long.Units = 'normalized';
fig_long.OuterPosition = [0.1 0.1 0.5 0.8]; %[0.1 0.3 0.8 0.5]
tt = tiledlayout(length(datadirs),1,'TileSpacing','none','Padding','none');

for d = 1:length(datadirs)
    %     close;
    datadir = datadirs{d};
    
    AFFT.datadir = datadir;
    AFFT.scan_type = scan_type;
    AFFT.trans_lims = trans_lims;
    
    if use_raw
        AFFT.fft_rawdataload();
    else
        AFFT.fft_dataload();
    end
    prop_distance_m = AFFT.propagation_distance/100; % propagation distance in m
    z_plot = AFFT.dtime + AFFT.simulation_window - AFFT.z;
    downsampleaverage = 1;
    
    %     longprofile = mean(reshape(AFFT.fft_densitymatrix(end,:),downsampleaverage,[]));
    longprofile = AFFT.fft_densitymatrix(end,:);
    %% plotazo
    
    ax(d) = nexttile;
    % plot limits
    %     hold on
    p = plot(z_plot,longprofile,'Linewidth',1);
    
    %         hold off
    
    if d < length(datadirs)
        ax(d).XTick = [];
        ax(d).XTickLabel = [];
    end
    ax(d).XDir = 'reverse';
    ax(d).YTick = [];
    ax(d).YTickLabel = [];
    legend(leg{d},'Location','Northwest','FontSize',12);
    legend('boxoff')
    xlim([min(z_plot),max(z_plot)])
    if d == 1
        max_ylim = 1.05*max(longprofile);
    end
    ylim([0,max_ylim])
    
    
    
end % for datadirs

%     grid on
xlabel(tt,'\xi (cm)','FontSize',15)
ylabel(tt,'charge (a.u.)','FontSize',15)
%     legend('simulation','experiment','location','best');
drawnow;

if saveplots
    plots_dir = ['gradsim_paper/long_profile/all/',''];
    plot_name = ['longprofile','all','n',num2str(dump),'charge'];
    P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'plot_name',plot_name,...
        'fig_handle',fig_long);
    P.save_plot();
end




