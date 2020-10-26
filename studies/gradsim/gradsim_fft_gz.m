%________________________________________________________________________
% freq in one slice in r, for all z, for all gradients, in one plot
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
% TODO (05/10/2020): add different options when plotting one or several gradients
% one: title with gradient, line with plasma frequency
% all: use tags jsut for all 
% Take code from fft_fz
%________________________________________________________________________

% load files
load('color_purple_to_green.mat');

% data directory
datadirs    = {'gm20','gm10'};
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
dataformat  = 'mat';
useAvg      = false;
dump_list   = 0:1:134;

% save directory
if length(datadirs) == 1
    plots_dir = ['fft/gz/',datadirs{1},''];
else
    plots_dir = ['fft/gz/','gradsim'];
end % if length datadirs
plot_name = ['fvsgz'];

plot_name_suffix = [''];
save_format = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';
field           = 'e';
direction       = 'z';
wakefields_direction = 'trans';

% analysis
xi_range        = [21 0.0];
% trans_range is set by trans_lims in this analysis
plasma_radius   = 0.15; % cm
trans_lims       = [0.066]; % cm

scan_type       = 'slice'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
freq_range      = [0.5 1.5]; % maybe norm. to plasma freq
max_dotsize     = 250;
showChargeinDotSize = true; % show dot size to reflect charge


% switches
plot_powerspectra   = false;
save_plot_flag      = true;

% calculated parameters
low_freqrange   = freq_range(1); % min frequency to plot (norm. to plasma freq)
upp_freqrange   = freq_range(2); % max frequency to plot (norm. to plasma freq)

% initialize variables
peak_freqs = zeros(length(datadirs),length(dump_list));
charge_in_slice = ones(length(datadirs),length(dump_list));
prop_distance_m = zeros(length(datadirs),length(dump_list));

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
    
    for n = 1:length(dump_list)
        
        AFFT.dump = dump_list(n);
        
        AFFT.fft_dataload();
        prop_distance_m(d,n) = AFFT.propagation_distance/100; % propagation distance in m
        
        % calculates fft and gives AFFT.fft_frequencies
        % and AFFT.fft_powerspectrum (_den or _fld)
        AFFT.get_fft();
        
        switch AFFT.property
            case 'density'
                charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
                powerspectrum_in_slice = AFFT.fft_powerspectrum_den;
            case 'fields'
                powerspectrum_in_slice = AFFT.fft_powerspectrum_fld;
        end % switch property
        
        AFFT.fft_peaks(powerspectrum_in_slice);
        % maxloc is the location in freq space of the ampltiude peak in
        % the power spectrum
        if isempty(AFFT.maxloc)
            peak_freqs(d,n) = 0;
        else
            peak_freqs(d,n) = AFFT.maxloc;
        end
        AFFT.progress_dump(['fft gz ',AFFT.datadir],n,length(dump_list));
    end %  for dump list
    
    
end % for datadirs

peak_freqs_plot = peak_freqs;
peak_freqs_plot(peak_freqs/(AFFT.plasmafreq_GHz) < low_freqrange) = nan;
peak_freqs_plot(peak_freqs/(AFFT.plasmafreq_GHz) > upp_freqrange) = nan;

%% Plotty
fig_gz = figure;
colororder(cc);
for d = 1:length(datadirs)
    hold on
    scatter(prop_distance_m(d,:),peak_freqs_plot(d,:)/AFFT.plasmafreq_GHz,'filled')
    hold off
end % for datadirs

grid on
xlabel('z (m)')
ylabel({'bunch core frequency /','plasma frequency at z = 0 m'});
title('trans. lim. = 0.66 mm')
legend(leg,'Location','southwest')

P.plot_name = plot_name;
P.save_plot(fig_gz);








