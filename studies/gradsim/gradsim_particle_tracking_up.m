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
% COMMENT ABOVE NOT VALID
% Apparently tracks for the up thing are saved in 1000 points

% clear;
close all;
% data location
datadirs = {'gm20'};
dataformat = 'mat';
dataformat_tracks = 'mat';
useAvg = false;

dump_fft = 60; % must be at 99

% simulation parameters
plasmaden = 1.81e14; %cm^-3
species = 'proton_beam';

% trans_range = 0.14 + [0.0 0.01]; % cm
trans_range = 0.08 + [0.0 0.01]; % cm
trans_lim = trans_range;
% FFT parameters
scan_type = 'slice';
on_axis = 'int';
xi_range = [21 0.74];
nslices = 1;

% plotting
P = Plotty();


for d = 1:length(datadirs)
    datadir = datadirs{d};
    P.plots_dir = ['gradsim_paper/tracking_up/',datadir];
%     
%     OPT = OsirisParticleTracking('datadir',datadir,'plasmaden',plasmaden,...
%         'dataformat',dataformat_tracks,...
%         'trans_range',trans_range,...
%         'property','tracks','trackfile_suffix','_up');
%     OPT.getdata();

    
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
    
    AFFT.fft_densitymatrix(1,:) = [];
    AFFT.fft_powerspectrum_den(1,:) = [];
    AFFT.fft_phase_den(1,:) = [];
    data_in_slice = AFFT.fft_powerspectrum_den;
    
    peak_freqs = zeros(nslices,length(datadirs));
    peak_amplitude = zeros(nslices,length(datadirs));
    peak_phase = zeros(nslices,length(datadirs));
    
    for r = 1:nslices
%         AFFT.maxloc = 115;
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
        
%         ind_r = find((OPT.denorm_distance(OPT.tracks_r(:,dump_fft*10+1)) < trans_range(2)) ...
%             & (OPT.denorm_distance(OPT.tracks_r(:,dump_fft*10+1)) > trans_range(1)));
%         
%         
%         par_dump_z = OPT.denorm_distance(OPT.tracks_z(ind_r,dump_fft*10+1));
%         par_dump_r = OPT.denorm_distance(OPT.tracks_r(ind_r,dump_fft*10+1));        
%         diff_matrix = par_dump_z - z_cos_positive{r,d};
%         
%         ind_diff = logical(sum(abs(diff_matrix) <= AFFT.dz/2,2)); %indices for particles which are on the positive zones of cosine
%         ind_final = ind_r(ind_diff);
%         
%         
%         par_z = OPT.denorm_distance(OPT.tracks_z(ind_final,1:1000)); % z position for particles in the positive region of the cosine from the fft
%         par_r = OPT.denorm_distance(OPT.tracks_r(ind_final,1:1000)); % r position for particles in the positive region of the cosine from the fft
%         par_q = OPT.tracks_q(ind_final,1);
%         par_q = repmat(par_q,1,size(par_r,2));
%         % nanino
%         par_q(par_r==0) = nan;
%         par_z(par_z==0) = nan;
%         par_r(par_r==0) = nan;
%         
%         
%         fig1 = figure(1);
%         
%         ind_plot = randi([1 size(par_z,1)],[150,1]);
%         
%         p1 = plot(par_z(ind_plot,:)'/100,10*par_r(ind_plot,:)');
%         p1colors = get(p1,'color');
%         alpha_q = par_q(ind_plot)/max(par_q(ind_plot));
%         p1colorsalpha = mat2cell([cell2mat(p1colors),alpha_q],ones(1,length(alpha_q)));
%         set(p1,{'color'},p1colorsalpha)
%         
%         rave = 10*sum((par_q.*par_r),'omitnan')./sum(par_q,'omitnan');
%         zave = sum((par_q.*par_z),'omitnan')./sum(par_q,'omitnan')/100;
% 
%         hold on
%         plot(zave,rave,'b','LineWidth',3);
%         hold off
% 
%         
%         xlabel('propagation direction (m)');
%         ylabel('transverse direction (mm)');
%         
%         xlim([0 10.28])
%         fig1.Units = 'normalized';
%         fig1.Position = [0.2 0.2 0.7 0.5];
%         
%         
%         P.plot_name = 'microbunches_tracks';
%         P.fig_handle = fig1;
%         P.save_plot();
%         pause
    end % for nslices
    
    
    
    low_lim = AFFT.fft_low_lim*AFFT.plasmafreq_GHz;
    upp_lim = AFFT.fft_upp_lim*AFFT.plasmafreq_GHz;
    
%     % transform freqs from Hz to GHz
    freqs = AFFT.fft_frequencies/1e9;
%     
    %% Find peaks
    % look for peaks within the set limits
    ind_pksearch = ((freqs > low_lim) & (freqs < upp_lim));

    norm_amplitude = max(AFFT.fft_densitymatrix(r,:))/peak_amplitude/2;
    fig22 = figure(22);
    hold on
    plot(AFFT.z,AFFT.fft_densitymatrix(r,:),'Linewidth',1)
    plot(AFFT.z,norm_amplitude*peak_amplitude*cos_signal{r,d},'Linewidth',2)
    hold off
    
    ylabel('density profile on axis (1/cm)');
    xlabel('prop. distance (cm)')
    xlim([629 633]);
    ylim([0 5e9]);
    
    fig22.Units = 'normalized';
    fig22.Position = [0.1 0.3 0.8 0.3];
    
    P.plot_name = 'microbunches_cosaxis';
    P.fig_handle = fig22;
    P.save_plot();
    %
    
    %
    %     figure
    %     plot(AFFT.fft_frequencies(ind_pksearch),AFFT.fft_phase_den(ind_pksearch))

    
end

fig888 = figure(888);
plot(AFFT.fft_frequencies(ind_pksearch)/1e9,AFFT.fft_powerspectrum_den(ind_pksearch),'Linewidth',2)
ylabel('Amplitude')
xlabel('frequency (GHz)')
xlim([60 180])
P.plot_name = 'fftgm20';
P.fig_handle = fig888;
P.save_plot();

fig999 = figure(999);
scatter(AFFT.fft_frequencies(ind_pksearch)/1e9,AFFT.fft_phase_den(ind_pksearch),'filled')
ylabel('Phase')
xlabel('frequency (GHz)')
xlim([115 125])
P.plot_name = 'phasegm20';
P.fig_handle = fig999;
P.save_plot();


