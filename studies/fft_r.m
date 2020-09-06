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


% data directory
datadirs = {'gm20'};

% parameters
plasma_density = 1.81e14;

% properties
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump = 100;
dataformat = 'h5';
useAvg = false;

% limits
plasma_radius = 0.15; % cm
nslices = 15;
xi_range = [20 0];



% analysis parameters
scan_type = 'slice'; % slice, cumulative
on_axis = 'sum'; % int, sum, intw, lineout
max_dotsize = 2*360;
showChargeinDotSize = false; % show dot size to reflect charge

% switches
plot_powerspectra = false;
saveplots = true;
normalized_frequency_switch = false;

% calculated parameters
trans_lims = (1:nslices)/nslices*plasma_radius;


% initialize variables
peak_freqs = zeros(nslices,1);
dotsize = max_dotsize/10*ones(nslices,1);
charge_in_slice = ones(nslices,1);

for d = 1:length(datadirs)
    datadir = datadirs{d};
    
    switch datadir
        case 'gm20'
            freq_lims = [108 124]; % GHz gm20
        case 'gp20'
            freq_lims = [120 128]; % GHz gp20
    end
    % directory to save plots
    plots_dir = ['plots/fft/r/grads/',datadir,''];
    
    AFFT = AwakeFFT('datadir',datadir,...
        'plasmaden',plasma_density,'property',property,'species',species,...
        'direction',direction,...
        'dump',dump,'dataformat',dataformat,'useAvg',useAvg,...
        'trans_lims',trans_lims,'xi_range',xi_range,...
        'scan_type',scan_type,'on_axis',on_axis);
    
    AFFT.fft_dataload();
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
    
    for r = 1:nslices
        
        AFFT.fft_peaks(AFFT.fft_frequencies,data_in_slice(r,:));
        % maxloc is the location in freq space of the ampltiude peak in
        % the power spectrum
        if isempty(AFFT.maxloc)
            peak_freqs(r) = 0;
        else
            peak_freqs(r) = AFFT.maxloc;
        end
        
    end % for nslices
    
    % Calculate DFT for whole range
    AFFT.scan_type = 'cumulative';
    AFFT.trans_lims = 0.13636; %three experimental sigmas at the end of the plasma 3*0.455
    AFFT.fft_dataload(); AFFT.get_fft();
    AFFT.fft_peaks(AFFT.fft_frequencies,AFFT.fft_powerspectrum_den);
    freq_widerange = AFFT.maxloc;
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize==0) = nan;
    end
    peak_freqs_plot = peak_freqs;
    
    %%
    figure(666)
    trans_lims_m = (trans_lims*10)'; % trans lims in mm
    hold on
    if normalized_frequency_switch
        scatter(trans_lims_m,...
            peak_freqs_plot/AFFT.plasmafreq_GHz,dotsize,'filled');
        ylim([min(peak_freqs_plot),max(peak_freqs_plot)*1.01]/AFFT.plasmafreq_GHz)
        ylabel('freq. / plasma freq.')
    else
        scatter(trans_lims_m,...
            peak_freqs_plot,dotsize,'filled');
        %         ylim([min(peak_freqs_plot)*0.99,max(peak_freqs_plot)*1.01])
        ylim(freq_lims);
        ylabel('freq. (GHz)')
    end % if normalized frequency switch
    ax_handle = gca;
    plot(3*0.2*sqrt(1+100/4.9^2)*ones(1,2),ax_handle.YLim,'--r','LineWidth',2) % dashed line at 3 sigma_r (bunch travelling in vacuum)
    plot(ax_handle.XLim,freq_widerange*ones(1,2),'--b','LineWidth',2);
    hold off
    grid on
    xlim([trans_lims_m(1) trans_lims_m(end)])
    xlabel('r (mm)')
    title('g = +20 %')
    
    
end % for datadirs







