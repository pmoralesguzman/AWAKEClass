%________________________________________________________________________
% Freq vs r for 1 gradient. 
% Comparison of simulation and experiment. Transverse limits are set by
% experimental data. 
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
% Last update: 05/10/2020
%________________________________________________________________________
clear;


% data directory
datadirs = {'gm20'};
% datadirs = {'gm20','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
dataformat  = 'mat';
useAvg      = false;
dump        = 133;

% save directory
plots_dir = ['gradsim_paper/fft/r/grads/',num2str(dump),''];

save_format = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';
field           = 'e';
direction       = 'z';
wakefields_direction = 'trans';

% analysis
xi_range        = [21 0.74]; % cm
% trans_range   is set by trans_lims in this analysis (calculated from exp
% data in this case)
trans_upperlimit   = 0.16; % cm
posinega        = 'n'; % choose positive or negative side for the experimental data)


scan_type       = 'slice'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
max_dotsize     = 150;

% switches
use_raw             = false;
plot_powerspectra   = false;
save_all_plots      = false;
save_plot_flag      = false;
showChargeinDotSize = false; % show dot size to reflect charge
normalized_frequency_switch = true;

% classes 
AFFT = AwakeFFT(...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'field',field,'direction',direction,'wakefields_direction',wakefields_direction,...
    'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis);

P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'save_format',save_format);


for d = 1:length(datadirs)
    
    %     close;
    datadir = datadirs{d};
    plot_name = ['fvsr',datadir,'r',num2str(trans_upperlimit*100),posinega];
    
    
    switch datadir
        case {'gm20'}
            %         freq_lims = [107 124]; % GHz gm20
            exp_g = 'gm19';
        case 'gp20'
            %         freq_lims = [120 128]; % GHz gp20
            exp_g = 'gp20';
        case {'g0','g0d2'}
            %         freq_lims = [120 121.5]; % GHz gp20
            exp_g = 'g0';
        case 'gm10'
            %         freq_lims = [108 124]; % GHz gp20
            exp_g = 'gm9';
        case {'gm5','gm5d2'}
            %         freq_lims = [108 124]; % GHz gp20
            exp_g = 'gm5';
        case {'gp5','gp5d2'}
            %         freq_lims = [120 124]; % GHz gp20
            exp_g = 'gp4';
        case 'gp10'
            %         freq_lims = [120 128]; % GHz gp20
            exp_g = 'gp9';
        case 'gp15'
            %         freq_lims = [120 128]; % GHz gp20
            exp_g = 'gp13';
        otherwise
            warning('gradient not found');
            exp_g = 'g0';
    end
    
    % load experimental data
    fexp = load('fvsgr_exp.mat');
    exp_freqsall = fexp.(exp_g);
    exp_radii = fexp.radius_mm;
    
    % choose transverse limits exactly as experiment
    switch posinega
        case 'p'
            ind_exp_radii = (abs(exp_radii) < trans_upperlimit*10) & (exp_radii > 0);
            exp_translims = abs(exp_radii(ind_exp_radii));
            exp_translims = exp_translims(1:2:end);
            trans_lims = exp_translims'/10;
        case 'n'
            ind_exp_radii = (abs(exp_radii) < trans_upperlimit*10) & (exp_radii < 0);
            exp_translims = abs(exp_radii(ind_exp_radii));
            exp_translims = exp_translims(1:2:end);
            trans_lims = fliplr(exp_translims'/10);
    end
    
    nslices = length(trans_lims);
    exp_freqs = exp_freqsall(ind_exp_radii);
    exp_freqs = exp_freqs(1:2:end);
    
    AFFT.datadir = datadir;
    AFFT.scan_type = scan_type;
    AFFT.trans_lims = trans_lims;
    
    if use_raw
        AFFT.fft_rawdataload();
    else
        AFFT.fft_dataload();
    end
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
    AFFT.on_axis = 'int';
    AFFT.trans_lims = trans_upperlimit; %three experimental sigmas at the end of the plasma 3*0.455
    if use_raw
        AFFT.fft_rawdataload();
    else
        AFFT.fft_dataload();
    end
    AFFT.get_fft();
    AFFT.fft_peaks(AFFT.fft_powerspectrum_den);
    freq_widerange = AFFT.maxloc;
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize == 0) = nan;
    else
        dotsize = max_dotsize*ones(length(charge_in_slice),1);
        dotsize(dotsize == 0) = nan;
    end % end if charge dot size
    
    %%
    fig_fvsr = figure;
    
    % plot limits
    trans_lims_mm = (10*trans_lims)'; % trans lims in mm
    xlimits = [0 1.01*trans_lims_mm(end)];
    
    hold on
    
    if normalized_frequency_switch
        scatter(trans_lims_mm,...
            peak_freqs/AFFT.plasmafreq_GHz,dotsize,'filled');
        scatter(exp_translims,exp_freqs/AFFT.plasmafreq_GHz,dotsize,'x','LineWidth',1); % experiment
        plot(xlimits,freq_widerange*ones(1,2)/AFFT.plasmafreq_GHz,'--b','LineWidth',2);
        freq_lims = [0.98*min([exp_freqs;peak_freqs]/AFFT.plasmafreq_GHz),...
            1.02*max([exp_freqs;peak_freqs]/AFFT.plasmafreq_GHz)];
        ylim(freq_lims)
        ylabel({'transverse slice frequency /','plasma frequency at z = 0 m'})
    else
        scatter(trans_lims_mm,peak_freqs,dotsize,'filled'); % simulation
        scatter(exp_translims,exp_freqs,dotsize,'x','LineWidth',1); % experiment
        plot(xlimits,freq_widerange*ones(1,2),'--b','LineWidth',2);
        freq_lims = [0.98*min([exp_freqs;peak_freqs]),1.02*max([exp_freqs;peak_freqs])];
        ylim(freq_lims);
        ylabel('freq. (GHz)')
    end % if normalized frequency switch
    
    hold off
    grid on
    xlim(xlimits)
    xlabel('r (mm)')
    legend('simulation','experiment','location','best');
    drawnow;
    
    P.plot_name = plot_name;
    P.fig_handle = fig_fvsr;
    P.save_plot();
    
end % for datadirs







