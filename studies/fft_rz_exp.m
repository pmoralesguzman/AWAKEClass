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
datadirs = {'gp20'};
plots_dir = ['fft/rz/grads/',datadirs{1},'_2'];
P = Plotty();


% parameters
plasma_density = 1.81e14;

% properties
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump_list = 60:1:60;
dataformat = 'h5';
useAvg = false;

% limits
plasma_radius = 0.16; % cm
nslices = 8;
xi_range = [14 0];

% analysis parameters
scan_type = 'slice'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout
low_freqrange = 0.98; % min frequency to plot (norm. to plasma freq)
upp_freqrange = 1.05; % max frequency to plot (norm. to plasma freq)
max_dotsize = 360;
showChargeinDotSize = true; % show dot size to reflect charge

% switches
plot_powerspectra = false;
save_all_plots = true;

% calculated parameters
trans_lims = (1:nslices)/nslices*plasma_radius;

% initialize variables
peak_freqs = zeros(nslices,length(dump_list));
dotsize = max_dotsize/10*ones(nslices,length(dump_list));
charge_in_slice = ones(nslices,length(dump_list));
prop_distance_m = zeros(1,length(dump_list));

AFFT = AwakeFFT(...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'direction',direction,'wakefields_direction','trans',...
    'dataformat',dataformat,'useAvg',useAvg,...
    'trans_lims',trans_lims,'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis,...
    'fft_low_lim',0.5,'fft_upp_lim',1.5);

for d = 1:length(datadirs)
    AFFT.datadir = datadirs{d};
    
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
            if save_all_plots
                trans_lims_plot_title = [0,trans_lims];
                fft_freqs_plot = AFFT.fft_frequencies(AFFT.ind_pksearch)/1e9;
                amplitude_plot = powerspectrum_in_slice(r,AFFT.ind_pksearch);
                P.plots_dir = [plots_dir,'/power_spectra'];
                P.plot_name = ['powerspectrum_dump',num2str(AFFT.dump),...
                    '_slice',num2str(r)];
                
                figx = figure('visible','off');
                plot(fft_freqs_plot,amplitude_plot,'LineWidth',2);
                title(['power spec. (prop. dist. = ',...
                    num2str(AFFT.propagation_distance/100,2),'m; r = [',...
                    num2str(trans_lims_plot_title(r)*10,2),', ',num2str(trans_lims_plot_title(r+1)*10,2),']mm)']);
                ylabel('amplitude');
                xlabel('frequency (GHz)');
                xlim([fft_freqs_plot(1),fft_freqs_plot(end)]);
                drawnow;
                P.fig_handle = figx;
                P.save_plot();
                close all;
                AFFT.progress_dump('save slice fft',r,nslices);
                
                indz = AFFT.z - AFFT.dtime > 7;
                z_plot = AFFT.z(indz);
                
                density_plot = AFFT.fft_densitymatrix(r,indz);
                P.plots_dir = [plots_dir,'/density_distributions'];
                P.plot_name = ['density_long_dump',num2str(AFFT.dump),...
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
                P.fig_handle = figy;
                P.save_plot();
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
        dotsize(dotsize==0) = nan;
    end
    
    peak_freqs_plot = peak_freqs;
    peak_freqs_plot(peak_freqs/(AFFT.plasmafreq_GHz) < low_freqrange) = nan;
    peak_freqs_plot(peak_freqs/(AFFT.plasmafreq_GHz) > upp_freqrange) = nan;
    
    
    %%
    fvsrz = figure(666);
    trans_lims_m = (trans_lims*10)'; % trans lims in mm
    prop_distance_m_mat = repmat(prop_distance_m,length(trans_lims),1);
    hold on
    for n = 1:length(dump_list)
        scatter(prop_distance_m_mat(:,n),trans_lims_m,...
            dotsize(:,n),peak_freqs_plot(:,n)/AFFT.plasmafreq_GHz,'filled')
    end
    hold off
    colormap(jet)
    xlim([prop_distance_m(1) prop_distance_m(end)])
    ylim([trans_lims_m(1),trans_lims_m(end)])
    
    xlabel('z (m)')
    ylabel('r (mm)')
    title('g = +20 %')
    cbar = colorbar;
    cbar.Label.String = 'frequency / plasma frequency';
    
    P.plot_name = 'fvsrz';
    P.fig_handle = fvsrz;
    
    P.save_plot();
    
end % for datadirs







