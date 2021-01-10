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
datadirs = {'g0zh','g0rh','g0','g0z2','g0r2','g0dt92'};
% datadirs = {'gm20','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
markers = ['+*oxsp'];

% simulation parameters
dump = 100;
dataformat = 'h5';
useAvg = false;
use_raw = false;

% properties
plasma_density = 1.81e14;
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';


% limits
xi_range = [21 0.75]; % cm

trans_lims = [0.018]*(1:8);

% analysis parameters
scan_type = 'slice'; % slice, cumulative
on_axis = 'sum'; % int, sum, intw, lineout
max_dotsize = 2*360;
showChargeinDotSize = false; % show dot size to reflect charge

% switches
plot_powerspectra = false;
saveplots = true;
normalized_frequency_switch = true;

AFFT = AwakeFFT(...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'field',field,'direction',direction,...
    'dump',dump,'dataformat',dataformat,'useAvg',useAvg,...
    'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis);


for d = 1:length(datadirs)
    %     close;
    
    switch datadirs{d}
        case 'g0'
            dataformat = 'mat';
        otherwise
            dataformat = 'h5';
    end
    AFFT.dataformat = dataformat;
    datadir = datadirs{d};
    nslices = length(trans_lims);
    AFFT.datadir = datadir;
    AFFT.scan_type = scan_type;
    AFFT.trans_lims = trans_lims;
    
    if use_raw
        AFFT.fft_rawdataload();
    else
        AFFT.fft_dataload(true);
    end
    prop_distance_m = AFFT.propagation_distance/100; % propagation distance in m
    
    % calculates fft and gives AFFT.fft_frequencies
    % and AFFT.fft_powerspectrum (_den or _fld)
    AFFT.get_fft();
    
    switch AFFT.property
        case 'density'
            charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
            data_in_slice = AFFT.fft_powerspectrum_den;
        case 'fields'
            data_in_slice = AFFT.fft_powerspectrum_fld;
    end % switch property
    
    % initialize variables
    peak_freqs = zeros(nslices,1);
    dotsize = max_dotsize/10*ones(nslices,1);
    
    for r = 1:nslices
        
        AFFT.fft_peaks(data_in_slice(r,:));
        % maxloc is the location in freq space of the ampltiude peak in
        % the power spectrum
        if isempty(AFFT.maxloc)
            peak_freqs(r) = nan;
        else
            peak_freqs(r) = AFFT.maxloc;
        end
        
    end % for nslices
    
%     % Calculate DFT for whole range
%     AFFT.scan_type = 'cumulative';
%     AFFT.on_axis = 'int';
%     AFFT.trans_lims = plasma_radius; %three experimental sigmas at the end of the plasma 3*0.455
%     if use_raw
%         AFFT.fft_rawdataload();
%     else
%         AFFT.fft_dataload(true);
%     end
%     AFFT.get_fft();
%     AFFT.fft_peaks(AFFT.fft_powerspectrum_den);
%     freq_widerange = AFFT.maxloc;
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize==0) = nan;
    end
    
    %%
    fig_fvsr = figure(1);
    % plot limits
    trans_lims_mm = (10*trans_lims)'; % trans lims in mm
    xlimits = [0 1.01*trans_lims_mm(end)];
    freq_lims = [0.98*min([peak_freqs]),1.02*max([peak_freqs])];

    hold on
    
    if normalized_frequency_switch
%         scatter(trans_lims_mm,...
%             peak_freqs/AFFT.plasmafreq_GHz,dotsize,'filled');
         plot(trans_lims_mm,...
            peak_freqs/AFFT.plasmafreq_GHz,['-.',markers(d)],'MarkerSize',8,'LineWidth',2);
%         plot(xlimits,freq_widerange*ones(1,2)/AFFT.plasmafreq_GHz,'--b','LineWidth',2);
        freq_lims = [0.98*min([peak_freqs]/AFFT.plasmafreq_GHz),...
            1.02*max([peak_freqs]/AFFT.plasmafreq_GHz)];
        ylim(freq_lims)
        ylabel({'transverse slice frequency /','plasma frequency at z = 0 m'})
    else
        scatter(trans_lims_mm,peak_freqs,dotsize,'filled'); % simulation
        plot(xlimits,freq_widerange*ones(1,2),'--b','LineWidth',2);
        freq_lims = [0.98*min([exp_freqs;peak_freqs]),1.02*max([exp_freqs;peak_freqs])];
        ylim(freq_lims);
        ylabel('freq. (GHz)')
    end % if normalized frequency switch
    
    hold off
    grid on
    xlim(xlimits)
    xlabel('r (mm)')
    drawnow;
   

    
end % for datadirs
legend(datadirs,'location','best')
if saveplots
    plots_dir = ['convergence/fft/r/grads/',num2str(dump),''];
    plot_name = ['fvsr','g0'];
    P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'plot_name',plot_name,...
        'fig_handle',fig_fvsr);
    P.save_plot();
end








