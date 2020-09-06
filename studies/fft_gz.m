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
property = 'fields';
species = 'proton_beam';
field = 'e';
direction = 'r';

% simulation parameters
dump_list = 0:1:100;
dataformat = 'h5';
useAvg = true;

% limits
plasma_radius = 0.15; % cm
nslices = 1;
xi_range = [20 0];

% analysis parameters
scan_type = 'slice'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout
low_freqrange = 0.8; % min frequency to plot (norm. to plasma freq)
upp_freqrange = 1.015; % max frequency to plot (norm. to plasma freq)
max_dotsize = 360;
showChargeinDotSize = false; % show dot size to reflect charge

% switches
plot_powerspectra = false;
saveplots = true;

% calculated parameters
trans_lims = 0.02;

% initialize variables
peak_freqs = zeros(nslices,length(dump_list));
dotsize = max_dotsize/10*ones(nslices,length(dump_list));
charge_in_slice = ones(nslices,length(dump_list));
prop_distance_m = zeros(1,length(dump_list));

for d = 1:length(datadirs)
    datadir = datadirs{d};
    
    % directory to save plots
    plots_dir = ['plots/fft/gz/grads/',datadir,''];
    if ~isfolder(plots_dir); mkdir(plots_dir); end
    
    for n = 1:length(dump_list)
        
        AFFT = AwakeFFT('datadir',datadir,...
            'plasmaden',plasma_density,'property',property,'species',species,...
            'direction',direction,'wakefields_direction','trans',...
            'dump',dump_list(n),'dataformat',dataformat,'useAvg',useAvg,...
            'trans_lims',trans_lims,'xi_range',xi_range,...
            'scan_type',scan_type,'on_axis',on_axis);
        
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
            % maxloc is the location in freq space of the ampltiude peak in
            % the power spectrum
            if isempty(AFFT.maxloc)
                peak_freqs(r,n) = 0;
            else
                peak_freqs(r,n) = AFFT.maxloc;
            end
            
        end % for nslices
        disp(n)
    end %  for dump list
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize==0) = nan;
    end
    
    peak_freqs_plot = peak_freqs;
    peak_freqs_plot(peak_freqs/(AFFT.plasmafreq_GHz) < low_freqrange) = nan;
    peak_freqs_plot(peak_freqs/(AFFT.plasmafreq_GHz) > upp_freqrange) = nan;
    
    
    %%
    figure(666)
    trans_lims_m = (trans_lims*10)'; % trans lims in mm
    prop_distance_m_mat = repmat(prop_distance_m,length(trans_lims),1);
    hold on
    scatter(prop_distance_m_mat,peak_freqs_plot,...
        'filled')
    hold off
    colormap(jet)
    xlim([prop_distance_m(1) prop_distance_m(end)])
    ylim([trans_lims_m(1),trans_lims_m(end)])
    
    xlabel('z (m)')
    ylabel('r (mm)')
    title('g = -20 %')
    cbar = colorbar;
    cbar.Label.String = 'frequency / plasma frequency';
    
    P = Plotty('plots_dir','fft/rz/grads/gm20','plot_name','fvsrz');
    P.save_plot(gcf);
    
end % for datadirs







