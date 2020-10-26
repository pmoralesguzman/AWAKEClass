%________________________________________________________________________
% Standard script for the AwakeFFT class.
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
datadirs    = {'gm20d2'};
dataformat  = 'h5';
useAvg      = true;
dump_list   = 0:1:134;

% save directory
plots_dir = ['gradsim_paper/fft/rz/',datadirs{1},''];
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
nslices         = 15;

scan_type       = 'slice'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
freq_range      = [0.86 1.025]; % maybe norm. to plasma freq 
max_dotsize     = 250;
showChargeinDotSize = true; % show dot size to reflect charge


% switches
plot_powerspectra   = false;
save_all_plots      = false;
save_plot_flag      = false;

% calculated parameters
trans_lims = (1:nslices)/nslices*plasma_radius;
low_freqrange   = freq_range(1); % min frequency to plot (norm. to plasma freq) 
upp_freqrange   = freq_range(2); % max frequency to plot (norm. to plasma freq) 

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
    
    for n = 1:length(dump_list)
        
    end % for dump list
end % for datadirs







