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
% datadirs = {'gm20','gm20d2'};
% datadirs = {'gm20','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'g0','gm20'};
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
plasma_radius = 0.25; % cm
xi_range = [21 0.0]; % cm
trans_lims = [0 0.066]; % cm
% analysis parameters
scan_type = 'cumulative'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout

% switches
saveplots = true;

AFFT = AwakeFFT(...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'field',field,'direction',direction,...
    'dump',dump,'dataformat',dataformat,'useAvg',useAvg,...
    'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis);


for d = 1:length(datadirs)
%     close;
    datadir = datadirs{d};
    nslices = length(trans_lims);
%     if d == 2; AFFT.dataformat = 'h5'; end
        
    AFFT.datadir = datadir;
    AFFT.scan_type = scan_type;
    AFFT.trans_lims = trans_lims;
    
    if use_raw
        AFFT.fft_rawdataload();
    else
        AFFT.fft_dataload();
    end
    prop_distance_m = AFFT.propagation_distance/100; % propagation distance in m
    if d == 1
        z_plot = AFFT.z; 
    else
        delta_z = AFFT.z(1)-z_plot(1);
        z_plot = AFFT.z - delta_z;
    end
    charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
    
    %%
    fig_long = figure(1);
    % plot limits
    hold on
    p = plot(z_plot,AFFT.fft_densitymatrix(end,:));
    fig_long.Units = 'normalized';
    fig_long.OuterPosition = [0.1 0.3 0.8 0.3]; %[0.1 0.3 0.8 0.5]
    hold off

%     grid on
    xlim([min(z_plot),max(AFFT.z)])
    xlabel('propagation distance (cm)')
%     legend('simulation','experiment','location','best');
    drawnow;
   
    if saveplots
        plots_dir = ['long_profile/',datadir,'/',''];
        plot_name = ['long_profile_',datadir,'n',num2str(dump)];
        P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'plot_name',plot_name,...
            'fig_handle',fig_long);
        P.save_plot();
    end
    
end % for datadirs







