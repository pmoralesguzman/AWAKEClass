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
grad_sim = -0.02;
% properties
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump_list = 0:1:100;
dataformat = 'h5';
useAvg = true;

% limits
plasma_radius = 0.02; % cm
nslices = 1;
xi_range = [20 0];

% analysis parameters
scan_type = 'cumulative'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout
low_freqrange = 0.88; % min frequency to plot (norm. to plasma freq)
upp_freqrange = 1.01; % max frequency to plot (norm. to plasma freq)
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
                peak_freqs_den(r,n) = 0;
            else
                peak_freqs_den(r,n) = AFFT.maxloc;
            end
        end
        AFFT.property = 'fields';
        AFFT.trans_lims = 0.01;
        
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
                peak_freqs_fld(r,n) = 0;
            else
                peak_freqs_fld(r,n) = AFFT.maxloc;
            end
        end
        
        disp(n)
    end %  for dump list
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize==0) = nan;
    end
    
    peak_freqs_den_plot = peak_freqs_den;
    peak_freqs_den_plot(peak_freqs_den/(AFFT.plasmafreq_GHz) < low_freqrange) = nan;
    peak_freqs_den_plot(peak_freqs_den/(AFFT.plasmafreq_GHz) > upp_freqrange) = nan;
    
    peak_freqs_fld_plot = peak_freqs_fld;
    peak_freqs_fld_plot(peak_freqs_fld/(AFFT.plasmafreq_GHz) < low_freqrange) = nan;
    peak_freqs_fld_plot(peak_freqs_fld/(AFFT.plasmafreq_GHz) > upp_freqrange) = nan;
    
    
    %%
    figure(666)
    trans_lims_m = (trans_lims*10)'; % trans lims in mm
    prop_distance_m_mat = repmat(prop_distance_m,length(trans_lims),1);
    hold on
    scatter(prop_distance_m_mat,peak_freqs_den_plot./AFFT.plasmafreq_GHz,...
        'filled')
    scatter(prop_distance_m_mat,peak_freqs_fld_plot./AFFT.plasmafreq_GHz,...
        'filled')
    plasmafreq_ini = AFFT.plasmafreq_GHz;
    AFFT.plasmaden = 1.81e14*(1+grad_sim*prop_distance_m_mat);
    AFFT.den2freq();
    fpe_plot = AFFT.plasmafreq_GHz/plasmafreq_ini;
    plot(prop_distance_m_mat,fpe_plot,'LineWidth',2)
    AFFT.plasmaden = 1.81e14;
    AFFT.den2freq();
    hold off
    %     colormap(jet)
    xlim([0 prop_distance_m(end)])
    ylim([low_freqrange,upp_freqrange])
    
    xlabel('z (m)')
    ylabel('freq. / plasma freq.')
    title('g = +20 %')
    grid on
    
    legend('microbunch freq.','wakefields freq.','plasma freq.')
    
    P = Plotty('plots_dir',['fft/fz/grads/',datadirs{1}],'plot_name','fvsz');
    P.fig_handle = gcf;
    P.save_plot();
    
end % for datadirs







