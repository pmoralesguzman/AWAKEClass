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
datadirs = {'gm20d2'};
% datadirs = {'gm20','gm10','gm5','g0','gp5','gp10','gp15','gp20'};

% simulation parameters
dump = 100;
dataformat = 'h5';
useAvg = true;

% properties
plasma_density = 1.81e14;
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';

% choose positive or negative side for the experimental data
posinega = 'p';

% limits
plasma_radius = 0.3; % cm
xi_range = [21 0.7414]; % cm

% analysis parameters
scan_type = 'slice'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout
max_dotsize = 2*360;
showChargeinDotSize = false; % show dot size to reflect charge

% switches
plot_powerspectra = false;
saveplots = false;
normalized_frequency_switch = false;

AFFT = AwakeFFT(...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'field',field,'direction',direction,...
    'dump',dump,'dataformat',dataformat,'useAvg',useAvg,...
    'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis);


for d = 1:length(datadirs)
    close;
    datadir = datadirs{d};
    
    % switch for datadir to choose freqs limits and
    switch datadir
        case {'gm20','gm20d2'}
            %         freq_lims = [107 124]; % GHz gm20
            exp_g = 'gm19';
        case 'gp20'
            %         freq_lims = [120 128]; % GHz gp20
            exp_g = 'gp20';
        case 'g0'
            %         freq_lims = [120 121.5]; % GHz gp20
            exp_g = 'g0';
        case 'gm10'
            %         freq_lims = [108 124]; % GHz gp20
            exp_g = 'gm9';
        case 'gm5'
            %         freq_lims = [108 124]; % GHz gp20
            exp_g = 'gm5';
        case 'gp5'
            %         freq_lims = [120 124]; % GHz gp20
            exp_g = 'gp4';
        case 'gp10'
            %         freq_lims = [120 128]; % GHz gp20
            exp_g = 'gp9';
        case 'gp15'
            %         freq_lims = [120 128]; % GHz gp20
            exp_g = 'gp13';
    end
    
    % load experimental data
    fexp = load('fvsgr_exp.mat');
    exp_freqsall = fexp.(exp_g);
    exp_radii = fexp.radius_mm;
    
    % choose transverse limits exactly as experiment
    switch posinega
        case 'p'
            ind_exp_radii = (abs(exp_radii) < plasma_radius*10) & (exp_radii > 0);
            exp_translims = abs(exp_radii(ind_exp_radii));
            trans_lims = exp_translims'/10;
        case 'n'
            ind_exp_radii = (abs(exp_radii) < plasma_radius*10) & (exp_radii < 0);
            exp_translims = abs(exp_radii(ind_exp_radii));
            trans_lims = fliplr(exp_translims'/10);
    end
    
    nslices = length(trans_lims);
    exp_freqs = exp_freqsall(ind_exp_radii);
    
    AFFT.datadir = datadir;
    AFFT.scan_type = scan_type;
    AFFT.trans_lims = trans_lims;

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
    
    % Calculate DFT for whole range
    AFFT.scan_type = 'cumulative';
    AFFT.trans_lims = plasma_radius; %three experimental sigmas at the end of the plasma 3*0.455
    AFFT.fft_dataload(); AFFT.get_fft();
    AFFT.fft_peaks(AFFT.fft_powerspectrum_den);
    freq_widerange = AFFT.maxloc;
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize==0) = nan;
    end
    
    %%
    fig_fvsr = figure(666);
    % plot limits
    trans_lims_mm = (10*trans_lims)'; % trans lims in mm
    xlimits = [0 1.01*trans_lims_mm(end)];
    freq_lims = [0.98*min([exp_freqs;peak_freqs]),1.02*max([exp_freqs;peak_freqs])];

    hold on
    
    if normalized_frequency_switch
        scatter(trans_lims_mm,...
            peak_freqs/AFFT.plasmafreq_GHz,dotsize,'filled');
        ylim([min(peak_freqs),max(peak_freqs)*1.01]/AFFT.plasmafreq_GHz)
        ylabel('freq. / plasma freq.')
    else
        scatter(trans_lims_mm,peak_freqs,dotsize,'filled'); % simulation
        scatter(exp_translims,exp_freqs,dotsize,'x','LineWidth',1); % experiment
        ylim(freq_lims);
        ylabel('freq. (GHz)')
    end % if normalized frequency switch
    
    plot(xlimits,freq_widerange*ones(1,2),'--b','LineWidth',2);
    hold off
    grid on
    xlim(xlimits)
    xlabel('r (mm)')
    legend('simulation','experiment','location','best');
    drawnow;
   
    if saveplots
        plots_dir = ['fft/r/grads/',num2str(dump),''];
        plot_name = ['fvsr',datadir,posinega];
        P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'plot_name',plot_name,...
            'fig_handle',fig_fvsr);
        P.save_plot();
    end
    
end % for datadirs







