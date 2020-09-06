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
dump_list = 0:1:100;
dataformat = 'mat';
useAvg = false;

% limits
plasma_radius = 0.15; % cm
nslices = 15;
xi_range = [20 0];

% analysis parameters
scan_type = 'slice'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout
low_freqrange = 1; % min frequency to plot (norm. to plasma freq)
upp_freqrange = 1.05; % max frequency to plot (norm. to plasma freq)
max_dotsize = 360;
showChargeinDotSize = false; % show dot size to reflect charge

% switches
plot_powerspectra = false;
saveplots = true;

% calculated parameters
trans_lims = (1:nslices)/nslices*plasma_radius;

% initialize variables
amp_peak_freqs = zeros(nslices,length(dump_list));
dotsize = max_dotsize/10*ones(nslices,length(dump_list));
charge_in_slice = ones(nslices,length(dump_list));
prop_distance_m = zeros(1,length(dump_list));

for d = 1:length(datadirs)
    datadir = datadirs{d};
    
    % directory to save plots
    plots_dir = ['plots/fft/rz/grads/',datadir,''];
    if ~isfolder(plots_dir); mkdir(plots_dir); end
    
    for n = 1:length(dump_list)
        
        AFFT = AwakeFFT('datadir',datadir,...
            'plasmaden',plasma_density,'property',property,'species',species,...
            'direction',direction,...
            'dump',dump_list(n),'dataformat',dataformat,'useAvg',useAvg,...
            'trans_lims',trans_lims,'xi_range',xi_range,...
            'scan_type',scan_type,'on_axis',on_axis,...
            'fft_low_lim',low_freqrange,'fft_upp_lim',upp_freqrange);
        
        AFFT.fft_dataload();
        prop_distance_m(n) = AFFT.propagation_distance/100; % propagation distance in m
        
        % calculates fft and gives AFFT.fft_frequencies
        % and AFFT.fft_powerspectrum (_den or _fld)
        AFFT.get_fft();
        
        switch AFFT.property
            case 'density'
                charge_in_slice(:,n) = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
                data_in_slice = AFFT.fft_powerspectrum_den;
            case 'fields'
                data_in_slice = AFFT.fft_powerspectrum_fld;
        end % switch property
        
        for r = 1:nslices
            
            AFFT.fft_peaks(AFFT.fft_frequencies,data_in_slice(r,:));
            % maxloc is the location in freq space of the peak in
            % the power spectrum
            if isempty(AFFT.maxpeak) 
                amp_peak_freqs(r,n) = nan;
            else
                amp_peak_freqs(r,n) = AFFT.maxpeak./trapz(AFFT.fft_frequencies,data_in_slice(r,:));
            end
            
        end % for nslices
        disp(n)
    end %  for dump list
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize==0) = nan;
    end
    
    amp_peak_plot = amp_peak_freqs;
    
    amp_peak_plot(amp_peak_freqs < 0.01*max(max(amp_peak_freqs))) = nan;
%     amp_peak_plot(amp_peak_freqs/(AFFT.plasmafreq_GHz) > upp_freqrange) = nan;
    
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
    %%
    figure(667)
    trans_lims_m = (trans_lims*10)'; % trans lims in mm
    prop_distance_m_mat = repmat(prop_distance_m,length(trans_lims),1);
    hold on
    for n = 1:length(dump_list)
        scatter(prop_distance_m_mat(:,n),trans_lims_m,...
            dotsize(:,n),log(amp_peak_plot(:,n)),'filled')
    end
    hold off
    colormap(jet)
    xlim([prop_distance_m(1),prop_distance_m(end)])
    ylim([trans_lims_m(1),trans_lims_m(end)])
    
    xlabel('z (m)')
    ylabel('r (mm)')
    title(title_g)
    cbar = colorbar;
    cbar.Label.String = 'highest amplitude of the fft in the freq. range';
    
    
end % for datadirs







