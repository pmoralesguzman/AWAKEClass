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

% load color order for 9 gradients
load('color_purple_to_green.mat');
load('color_red_to_blue.mat');

markers = {'o','+','x','s','d','p','h','.','*'};
% data directory
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'gp20'};
% parameters
plasma_density = 1.81e14;

% properties
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump_list = 59:2:100;
dataformat = 'h5';
useAvg = true;

% limits
nslices = 1;
xi_range = [14 0.384];

% analysis parameters
scan_type = 'cumulative'; % slice, cumulative
on_axis = 'sum'; % int, sum, intw, lineout
low_freqrange = 0.85; %0.88; % min frequency to plot (norm. to plasma freq)
upp_freqrange = 1.15; %1.01 % max frequency to plot (norm. to plasma freq)
max_dotsize = 360;
showChargeinDotSize = false; % show dot size to reflect charge

% switches
plot_powerspectra = false;
saveplots = false;

% calculated parameters
trans_lims = 0.02;

% initialize variables
peak_freqs = zeros(nslices,length(dump_list));
dotsize = max_dotsize/10*ones(nslices,length(dump_list));
charge_in_slice = ones(nslices,length(dump_list));
prop_distance_m = zeros(1,length(dump_list));
peak_freqs_den = zeros(1,length(dump_list));
peak_freqs_fld = zeros(1,length(dump_list));


for d = 1:9 %length(datadirs)
    datadir = datadirs{d};
    
    switch datadir
        case {'gp20','gp20d2'}
            title_g = 'g = +2 %/m';
            grad_sim = 0.02;
        case 'gp15'
            title_g = 'g = +1.5 %/m';
            grad_sim = 0.015;
        case 'gp10'
            title_g = 'g = +1 %/m';
            grad_sim = 0.01;
        case 'gp5'
            title_g = 'g = +0.5 %/m';
            grad_sim = 0.005;
        case {'g0','g0d2'}
            title_g = 'g = 0 %/m';
            grad_sim = 0.0;
        case 'gm5'
            title_g = 'g = -0.5 %/m';
            grad_sim = -0.005;
        case 'gm10'
            title_g = 'g = -1 %/m';
            grad_sim = -0.01;
        case 'gm15'
            title_g = 'g = -1.5 %/m';
            grad_sim = -0.015;
        case 'gm20'
            title_g = 'g = -2 %/m';
            grad_sim = -0.02;
        case 'R2gap0_2e14'
            title_g = 'density step';
    end % switch datadir
    
    AFFT = AwakeFFT('datadir',datadir,...
        'plasmaden',plasma_density,'property','density','species',species,...
        'direction',direction,'wakefields_direction','long',...
        'dump',dump_list(1),'dataformat',dataformat,'useAvg',useAvg,...
        'trans_lims',trans_lims,'xi_range',xi_range,...
        'scan_type',scan_type,'on_axis',on_axis);
    
    AFFT.fft_dataload();
    AFFT.get_fft();
    
    freq_ind = ((AFFT.fft_frequencies*1e-9)/AFFT.plasmafreq_GHz > 0.89) & (AFFT.fft_frequencies*1e-9/AFFT.plasmafreq_GHz < 1.1);
    
    powerspectra_den = zeros(sum(freq_ind),length(dump_list));
    powerspectra_fld = powerspectra_den;
    
    for n = 1:length(dump_list)
        
        AFFT.dump = dump_list(n);
        AFFT.property = 'density';
        AFFT.trans_lims = trans_lims;
        
        AFFT.fft_dataload();
        prop_distance_m(n) = AFFT.propagation_distance/100; % propagation distance in m
        
        % calculates fft and gives AFFT.fft_frequencies
        % and AFFT.fft_powerspectrum (_den or _fld)
        AFFT.get_fft();
        
        charge_in_slice(:,n) = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
        powerspectra_den(:,n) = AFFT.fft_powerspectrum_den(freq_ind);
        
        for r = 1:nslices
            
            AFFT.fft_peaks(AFFT.fft_powerspectrum_den(r,:));
            % max_loc is the location in freq space of the ampltiude peak in
            % the power spectrum
            if isempty(AFFT.max_loc)
                peak_freqs_den(r,n) = 0;
            else
                peak_freqs_den(r,n) = AFFT.max_loc;
            end
        end
        
        
        
                AFFT.property = 'fields';
                AFFT.trans_lims = 0.005;
                if dump_list(n) <= 100 % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    AFFT.fft_dataload();
                    prop_distance_m(n) = AFFT.propagation_distance/100; % propagation distance in m
        
                    % calculates fft and gives AFFT.fft_frequencies
                    % and AFFT.fft_powerspectrum (_den or _fld)
                    AFFT.get_fft();
        
                    powerspectra_fld(:,n) = AFFT.fft_powerspectrum_fld(freq_ind);
        
                    for r = 1:nslices
        
                        AFFT.fft_peaks(AFFT.fft_powerspectrum_fld(r,:));
                        % max_loc is the location in freq space of the ampltiude peak in
                        % the power spectrum
                        if isempty(AFFT.max_loc)
                            peak_freqs_fld(r,n) = 0;
                        else
                            peak_freqs_fld(r,n) = AFFT.max_loc;
                        end
                    end
        
                end
        AFFT.progress_dump('frequency',n,length(dump_list))
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
    figure(1)
    trans_lims_m = (trans_lims*10)'; % trans lims in mm
    prop_distance_m_mat = repmat(prop_distance_m,length(trans_lims),1);
    hold on
    AFFT.plasmaden = plasma_density;
    AFFT.den2freq();
    plot(prop_distance_m_mat,peak_freqs_den_plot./AFFT.plasmafreq_GHz,...
        'Color',cc(d,:),'MarkerEdgeColor',cc(d,:),'MarkerFaceColor',cc(d,:),...
        'Marker',markers{d});
        plot(prop_distance_m_mat,peak_freqs_fld_plot./AFFT.plasmafreq_GHz,...
            'Color',ccrb(d,:),'MarkerEdgeColor',ccrb(d,:),'MarkerFaceColor',ccrb(d,:),...
            'Marker',markers{10-d})
    plasmafreq_ini = AFFT.plasmafreq_GHz;
    AFFT.plasmaden = plasma_density*(1+grad_sim*prop_distance_m_mat(prop_distance_m_mat < 10.2)); % !!!!!!!!!!!!!!!!!!!!
    
    AFFT.den2freq();
    fpe_plot = AFFT.plasmafreq_GHz/plasmafreq_ini;
    plot(prop_distance_m_mat(prop_distance_m_mat < 10.2),fpe_plot,'LineWidth',2,'Color',[0.2,0.2,0.2]) % !!!!!!!!!!!!!!!!!!!!
    AFFT.plasmaden = plasma_density;
    AFFT.den2freq();
    
    xlim([6 prop_distance_m(end)])
    ylim([0.89,1.1])
    
    xlabel('z (m)')
    ylabel('frequency / plasma frequency at z = 0')
    %     title([title_g,'  r = 0.1 mm'])
    grid on
    
        legend('microbunch freq.','wakefields freq.','plasma freq.','location','best')
    %     legend('microbunch freq.','plasma freq.')
    
    %     close
    
    %     figure % waterfall density
    %     imagesc(prop_distance_m,AFFT.fft_frequencies(freq_ind)*1e-9/AFFT.plasmafreq_GHz,powerspectra_den);
    %     title('microbunch density')
    %     ylabel('norm. frequencies')
    %     xlabel('z (m)')
    %     set(gca,'YDir','normal')
    %     colorbar;
    %
    %     figure % waterfall fields
    %     imagesc(prop_distance_m,AFFT.fft_frequencies(freq_ind)*1e-9/AFFT.plasmafreq_GHz,powerspectra_fld);
    %     title('longitudinal wakefields')
    %     ylabel('norm. frequencies')
    %     xlabel('z (m)')
    %     set(gca,'YDir','normal')
    %     colorbar;
    
    
end % for datadirs


drawnow;
P = Plotty('plots_dir',['gradsim_ipac/fft/fzipac/all',datadir],'plot_name','fvszboth100fld');
P.fig_handle = gcf;
P.save_plot();




