%________________________________________________________________________
% Frequency Map in r and z of 1 gradient
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
% Last update: 04/10/2020
%________________________________________________________________________


% data directory
datadirs    = {'R2gap0_2e14'};
dataformat  = 'h5';
useAvg      = false;
dump_list   = 0:1:200;

% save directory
plots_dir = ['test/fft/rz/',datadirs{1},''];
plot_name_suffix = [''];
save_format = {'png','eps','fig'};

% properties
plasma_density  = 2e14;
property        = 'density';
species         = 'proton_beam';
field           = 'e';
direction       = 'z';
wakefields_direction = 'trans';

% analysis
xi_range        = [6 0.0];
% trans_range is set by trans_lims in this analysis
plasma_radius   = 0.15; % cm
nslices         = 15;

scan_type       = 'slice'; % slice, cumulative
on_axis         = 'int'; % int, sum, intw, lineout
freq_range      = [0.88 1.1];
max_dotsize     = 250;
showChargeinDotSize = true; % show dot size to reflect charge


% switches
plot_powerspectra   = false;
save_all_plots      = false;
save_plot_flag      = true;

% calculated parameters
trans_lims = (1:nslices)/nslices*plasma_radius;

% initialize variables
peak_freqs = zeros(nslices,length(dump_list));
dotsize = max_dotsize/10*ones(nslices,length(dump_list));
charge_in_slice = ones(nslices,length(dump_list));
prop_distance_m = zeros(1,length(dump_list));

AFFT = AwakeFFT(...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'field',field,'direction',direction,'wakefields_direction',wakefields_direction,...
    'xi_range',xi_range,'trans_lims',trans_lims,...
    'scan_type',scan_type,'on_axis',on_axis); %...
%     'fft_low_lim',low_freqrange,'fft_upp_lim',upp_freqrange,...
%     'peak_distance',abs(diff([low_freqrange,upp_freqrange]))/10);

P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'save_format',save_format,'save_flag',save_plot_flag);

for d = 1:length(datadirs)
    AFFT.datadir = datadirs{d};
    if length(AFFT.datadir) > 4
       gradsim_datadir_length = 4;
    else
        gradsim_datadir_length =  length(AFFT.datadir);
    end % if length datadir
    
    if isempty(freq_range)
        switch AFFT.datadir(1:gradsim_datadir_length)
            case 'gp20'
                freq_range = [0.97 1.2];
            case 'gp15'
            case 'gp10'
            case 'gp5'
            case {'g0','g0d2'}
                freq_range = [0.98 1.05];
            case 'gm5'
            case 'gm10'
            case 'gm15'
            case 'gm20'
                freq_range = [0.89 1.02];
        end % switch datadir
        
    end % if isempty freq_range
    low_freqrange = freq_range(1);
    upp_freqrange = freq_range(2);
    
    for n = 1:length(dump_list)
        AFFT.dump = dump_list(n);
        
        AFFT.fft_dataload();
        prop_distance_m(n) = AFFT.propagation_distance/100; % propagation distance in m

        % calculates fft and gives AFFT.fft_frequencies
        % and AFFT.fft_powerspectrum (_den or _fld)
        AFFT.get_fft();
        
        switch AFFT.property
            case 'density'
                charge_in_slice(:,n) = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
                powerspectrum_in_slice = AFFT.fft_powerspectrum_den;
                fft_phase = AFFT.fft_phase_den;
                profile_for_phase = AFFT.fft_densitymatrix;
            case 'fields'
                powerspectrum_in_slice = AFFT.fft_powerspectrum_fld;
                fft_phase = AFFT.fft_phase_fld;
                profile_for_phase = AFFT.fft_fieldmatrix;

        end % switch property
        
        for r = 1:nslices
            
            AFFT.fft_peaks(powerspectrum_in_slice(r,:));
            
            %%% Plots of the power spectrum.
            if plot_powerspectra
                P2 = P;
                trans_lims_plot_title = [0,trans_lims];
                fft_freqs_plot = AFFT.fft_frequencies(AFFT.ind_pksearch)/1e9;
                amplitude_plot = powerspectrum_in_slice(r,AFFT.ind_pksearch);
                P2.plots_dir = [plots_dir,'/power_spectra'];
                P2.plot_name = ['powerspectrum_dump',num2str(AFFT.dump),...
                    '_slice',num2str(r)];
                
                fig_powerspectra = figure('visible','off');
                plot(fft_freqs_plot,amplitude_plot,'LineWidth',2);
                title(['power spec. (prop. dist. = ',...
                    num2str(AFFT.propagation_distance/100,2),'m; r = [',...
                    num2str(trans_lims_plot_title(r)*10,2),', ',num2str(trans_lims_plot_title(r+1)*10,2),']mm)']);
                ylabel('amplitude');
                xlabel('frequency (GHz)');
                xlim([fft_freqs_plot(1),fft_freqs_plot(end)]);
                drawnow;
                P2.fig_handle = fig_powerspectra;
                P2.save_plot();
                close all;
                AFFT.progress_dump('save slice fft',r,nslices);
                
                indz = AFFT.z - AFFT.dtime > 7;
                z_plot = AFFT.z(indz);
                
                density_plot = AFFT.fft_densitymatrix(r,indz);
                P2.plots_dir = [plots_dir,'/density_distributions'];
                P2.plot_name = ['density_long_dump',num2str(AFFT.dump),...
                    '_slice',num2str(r)];
                
                figy = figure('visible','off');
                plot(z_plot,density_plot,'LineWidth',2);
                title(['long. den. distributon (prop. dist. = ',...
                    num2str(AFFT.propagation_distance/100,2),'m; r = [',...
                    num2str(trans_lims_plot_title(r)*10,2),', ',num2str(trans_lims_plot_title(r+1)*10,2),']mm)']);
                ylabel('density (1/cm)');
                xlabel('z-ct (m)');
                xlim([z_plot(1),z_plot(end)]);
                figy.Units = 'normalized';
                figy.Position = [0.05 0.2 0.9 0.5];
                drawnow;
                P2.fig_handle = figy;
                P2.save_plot();
                close all;
                AFFT.progress_dump('save slice density',r,nslices); 
            end
            
            
            % maxloc is the location in freq space of the ampltiude peak in
            % the power spectrum
            if isempty(AFFT.maxloc)
                peak_freqs(r,n) = 0;
            else
                peak_freqs(r,n) = AFFT.maxloc;
            end
            
        end % for nslices
        
        AFFT.progress_dump('dump',n,length(dump_list));
    end %  for dump list
    
    if showChargeinDotSize
        
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        if strcmp(AFFT.on_axis,'sum')
            size_limit = median(charge_in_slice,'all') + 1*std(charge_in_slice,0,'all');
            charge_for_size = charge_in_slice;
            charge_for_size(charge_in_slice > size_limit) = size_limit;
            norm_size = charge_for_size/max(charge_for_size,[],'all');
            dotsize = max_dotsize*norm_size;
        end
        dotsize(dotsize == 0) = nan;
    end
    
    %% plot section
    peak_freqs_plot = peak_freqs;
    peak_freqs_plot(peak_freqs/(AFFT.plasmafreq_GHz) < low_freqrange) = nan;
    peak_freqs_plot(peak_freqs/(AFFT.plasmafreq_GHz) > upp_freqrange) = nan;
    
    fvsrz = figure(667);
    trans_lims_m = (trans_lims*10)'; % trans lims in mm
    prop_distance_m_mat = repmat(prop_distance_m,length(trans_lims),1);
    hold on
    for n = 1:length(dump_list)
        sc_handle = scatter(prop_distance_m_mat(:,n),trans_lims_m,...
            dotsize(:,n),peak_freqs_plot(:,n)/AFFT.plasmafreq_GHz,'filled');
    end
    hold off
    colormap(jet)
    xlim([prop_distance_m(1),prop_distance_m(end)])
    ylim([trans_lims_m(1),trans_lims_m(end)])
    fvsrz.CurrentAxes.FontSize = 15; 
    xlabel('z (m)')
    ylabel('r (mm)')
%     title('g = -20 %')
    cbar = colorbar;
    cbar.Label.String = {'transverse slice frequency /',' plasma frequency at z = 0 m'};
    fvsrz.Units = 'normalized';
    fvsrz.Position = [0.1 0.2 0.7 0.4];
    drawnow;
    
    P.plot_name = ['fvsrz',property,AFFT.datadir,'n',num2str(dump_list(end))];
    P.fig_handle = fvsrz;
    
    P.save_plot();
    
end % for datadirs


