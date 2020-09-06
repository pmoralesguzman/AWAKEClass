%__________________________________________________________________________
% Script that takes tracks according to some conditions and plots them.
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 08/07/2020
%__________________________________________________________________________
% comment: particle tracks are saved each 5 dumps, they start at dump 1,
% the second save is dump 4, third save dump 9 and so on. So the last point
% for the gradients which has 200 points corresponds to dump 99.5. To have
% corresponding values, one should use dump 99 with point 199. 

% clear;
close all;
% data location
datadirs = {'gp20'};
dataformat = 'mat';
useAvg = false;

dump_fft = 99; % must be at 99

% simulation parameters
plasmaden = 1.81e14; %cm^-3
species = 'proton_beam';

trans_ranges = (0:0.01:0.14)' + [0.0 0.01]; % cm

% FFT parameters
scan_type = 'slice';
on_axis = 'int';
xi_range = [21 0.74];
nslices = 1;

% plotting
P = Plotty();
% color_meanline = [
% 107,174,214
% 66,146,198
% 251,106,74
% 245,90,60
% 239,59,44
% 252,187,161
% 252,146,114
% 252,187,161
% 252,146,114
% 251,106,74
% 245,90,60
% 239,59,44
% 203,24,29
% 165,15,21
% 103,0,13]/256;
% color_meanline = ['kkkkkkkkkkkkkkk'];
color_meanline = [
229,245,224
199,233,192
161,217,155
116,196,118
65,171,93
35,139,69
0,109,44
0,68,27
253,141,60
241,105,19
217,72,1
252,146,114
251,106,74
239,59,44
203,24,29]/256;

% OPT = OsirisParticleTracking('datadir',datadir,'plasmaden',plasmaden,...
%         'dataformat','mat',...
%         'property','tracks');
% OPT.getdata();


for d = 1:length(datadirs)
    datadir = datadirs{d};
    P.plots_dir = ['tracking/gradsim/',datadir];
    for rr = 1:length(trans_ranges)
        trans_range = trans_ranges(rr,:);
        trans_lim = trans_range;
        OPT.trans_range = trans_range;

    AFFT = AwakeFFT('datadir',datadir,...
        'plasmaden',plasmaden,'property','density','species',species,...
        'dump',dump_fft,'dataformat',dataformat,'useAvg',useAvg,...
        'trans_lims',trans_lim,'xi_range',xi_range,...
        'scan_type',scan_type,'on_axis',on_axis);
    
    
    AFFT.fft_dataload();
    prop_distance_m = AFFT.propagation_distance/100; % propagation distance in m
    
    % calculates fft and gives AFFT.fft_frequencies
    % and AFFT.fft_powerspectrum (_den or _fld)
    AFFT.get_fft();
    
    charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
    AFFT.fft_densitymatrix(1,:) = [];
    AFFT.fft_powerspectrum_den(1,:) = [];
    AFFT.fft_phase_den(1,:) = [];
    data_in_slice = AFFT.fft_powerspectrum_den;
    
    peak_freqs = zeros(nslices,length(datadirs));
    peak_amplitude = zeros(nslices,length(datadirs));
    peak_phase = zeros(nslices,length(datadirs));
    
    for r = 1:nslices
        
        AFFT.fft_peaks(data_in_slice(r,:),AFFT.fft_phase_den(r,:),AFFT.fft_densitymatrix(r,:));
        % maxloc is the location in freq space of the ampltiude peak in
        % the power spectrum
        if isempty(AFFT.maxloc)
            peak_freqs(r,d) = 0;
            peak_amplitude(r,d) = 0;
            peak_phase(r,d) = 0;
        else
            peak_freqs(r,d) = AFFT.maxloc;
            peak_amplitude(r,d) = AFFT.maxpeak;
            peak_phase(r,d) = AFFT.maxphase;
        end
        
        wavenumber(r,d) = peak_freqs(r,d)*1e9/AFFT.c_cm;
        
        cos_signal{r,d} = cos(2*pi*AFFT.z*wavenumber(r,d) + peak_phase(r,d));
        z_cos_positive{r,d} = AFFT.z(cos_signal{r,d} > 0);
        
        ind_r = find((OPT.denorm_distance(OPT.tracks_r(:,end-1)) < trans_range(2)) ...
            & (OPT.denorm_distance(OPT.tracks_r(:,end-1)) > trans_range(1)));
        
        
        par_last_z = OPT.denorm_distance(OPT.tracks_z(ind_r,end-1));
        par_last_r = OPT.denorm_distance(OPT.tracks_r(ind_r,end-1));        
        diff_matrix = par_last_z - z_cos_positive{r,d};
        
        ind_diff = logical(sum(abs(diff_matrix) <= AFFT.dz/2,2)); %indices for particles which are on the positive zones of cosine
        ind_final = ind_r(ind_diff);
        
        
        par_z = OPT.denorm_distance(OPT.tracks_z(ind_final,1:199)); % z position for particles in the positive region of the cosine from the fft
        par_r = OPT.denorm_distance(OPT.tracks_r(ind_final,1:199)); % r position for particles in the positive region of the cosine from the fft
        par_q = OPT.tracks_q(ind_final,1);
        
        fig1 = figure(1);
        
%         ind_plot = randi([1 size(par_z,1)],[100,1]);
%         
%         p1 = plot(par_z(ind_plot,:)'/100,10*par_r(ind_plot,:)');
%         p1colors = get(p1,'color');
%         alpha_q = par_q(ind_plot)/max(par_q(ind_plot));
%         p1colorsalpha = mat2cell([cell2mat(p1colors),alpha_q],ones(1,length(alpha_q)));
%         set(p1,{'color'},p1colorsalpha)
        
        rave = 10*sum((par_q.*par_r))/sum(par_q);
        zave = sum((par_q.*par_z))/sum(par_q)/100;

        hold on
        pave = plot(zave,rave,'color',color_meanline(rr,:),'LineWidth',2);
        hold off

        
        xlabel('propagation direction (m)');
        ylabel('transverse direction (mm)');
        
        xlim([0 10.28])
        fig1.Units = 'normalized';
        fig1.Position = [0.2 0.2 0.7 0.5];
        
        
        P.plot_name = 'microbunches_mean_tracks';
        P.fig_handle = fig1;
        P.save_plot();
        
    end % for nslices
    
    
    
    low_lim = AFFT.fft_low_lim*AFFT.plasmafreq_GHz;
    upp_lim = AFFT.fft_upp_lim*AFFT.plasmafreq_GHz;
    
%     % transform freqs from Hz to GHz
    freqs = AFFT.fft_frequencies/1e9;
%     
    %% Find peaks
%     % look for peaks within the set limits
%     ind_pksearch = ((freqs > low_lim) & (freqs < upp_lim));
% 
%     norm_amplitude = max(AFFT.fft_densitymatrix(r,:))/peak_amplitude/2;
%     fig22 = figure(22);
%     hold on
%     plot(AFFT.z,AFFT.fft_densitymatrix(r,:),'Linewidth',2)
%     plot(AFFT.z,norm_amplitude*peak_amplitude*cos_signal{r,d},'Linewidth',2)
%     hold off
%     
%     ylabel('density profile on axis (1/cm)');
%     xlabel('prop. distance (cm)')
%     xlim([1025.5 1030.35]);
%     ylim([0 5e9]);
%     
%     fig22.Units = 'normalized';
%     fig22.Position = [0.1 0.3 0.8 0.3];
%     
%     P.plot_name = 'microbunches_cosaxis';
%     P.fig_handle = fig22;
%     P.save_plot();
    %
    
    %
    %     figure
    %     plot(AFFT.fft_frequencies(ind_pksearch),AFFT.fft_phase_den(ind_pksearch))

    end
end
% 
% fig888 = figure(888);
% plot(AFFT.fft_frequencies(ind_pksearch)/1e9,AFFT.fft_powerspectrum_den(ind_pksearch),'Linewidth',2)
% ylabel('Amplitude')
% xlabel('frequency (GHz)')
% xlim([60 180])
% P.plot_name = 'fftgm20';
% P.fig_handle = fig888;
% P.save_plot();
% 
% fig999 = figure(999);
% scatter(AFFT.fft_frequencies(ind_pksearch)/1e9,AFFT.fft_phase_den(ind_pksearch),'filled')
% ylabel('Phase')
% xlabel('frequency (GHz)')
% xlim([115 125])
% P.plot_name = 'phasegm20';
% P.fig_handle = fig999;
% P.save_plot();


