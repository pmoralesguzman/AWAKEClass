%________________________________________________________________________
% FFT of the proton or plasma electrons density distribution, or of a
% lineout of the wakefields. Special version to produce the plot for
% the paper, which shows the fft for a narrow and wide window.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 14/06/2020
%________________________________________________________________________

close all;
clear;

% load experimental data from Fabian
freqs_exp      = load('gradsim_freqs.txt');
grads_exp      = freqs_exp(:,1)/10;
freqs_exp(:,1) = [];
% 1: full green, 2: CTR black, 3: narrow red, 4: wide blue

% data directory
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'gm5'};
% parameters
plasma_density = 1.81e14;
grads_sim = [-20,-15,-10,-5,0,5,10,15,20]/10;

% properties
property = 'density';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump = 133;
dataformat = 'mat';
useAvg = false;

% limits
plasma_radius = 0.198; % cm (4*0.455, following exp, see below trans_lims)
nslices = 2;
xi_range = [21 0.744];
freq_lims = [107 128]; % GHz

% analysis parameters
scan_type = 'cumulative'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout
max_dotsize = 2*360;
showChargeinDotSize = true; % show dot size to reflect charge

% switches
plot_powerspectra = false;
saveplots = true;
normalized_frequency_switch = true;

% calculated parameters
trans_lims = [0.066, plasma_radius];
% measured radius = 0.66 mm, first exp limit = 0.868 mm
% last exp limit 2.604 mm
% theo radius at plasma exit = 0.0455, 0.0455*1.3152 = 0.0598 cm

% initialize variables
peak_freqs = zeros(nslices,1);
dotsize = max_dotsize/10*ones(nslices,1);
charge_in_slice = ones(nslices,1);

for d = 1:length(datadirs)
    datadir = datadirs{d};
    
    % directory to save plots
    plots_dir = ['fft/rg/',''];    
    
    AFFT = AwakeFFT('datadir',datadir,...
        'plasmaden',plasma_density,'property',property,'species',species,...
        'direction',direction,...
        'dump',dump,'dataformat',dataformat,'useAvg',useAvg,...
        'trans_lims',trans_lims,'xi_range',xi_range,...
        'scan_type',scan_type,'on_axis',on_axis);
    
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
    
    for r = 1:nslices
        
        AFFT.fft_peaks(data_in_slice(r,:));
        % maxloc is the location in freq space of the ampltiude peak in
        % the power spectrum
        if isempty(AFFT.maxloc)
            peak_freqs(r,d) = 0;
        else
            peak_freqs(r,d) = AFFT.maxloc;
        end
        
    end % for nslices
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize==0) = nan;
    end
    disp(d)
end % for datadirs

peak_freqs_plot = peak_freqs;

%% plotty
fig33 = figure(33);
hold on
if normalized_frequency_switch
    % simulation
    plot_sim(1) = plot(grads_sim,... % narrow trans. range
        peak_freqs_plot(1,:)/AFFT.plasmafreq_GHz,'o');
    set(plot_sim(1), 'markerfacecolor', get(plot_sim(1), 'color'));
    plot_sim(2) = plot(grads_sim,... % wide trans. range
        peak_freqs_plot(2,:)/AFFT.plasmafreq_GHz,'o');
    set(plot_sim(2), 'markerfacecolor', get(plot_sim(2), 'color'));
    
    % experiment
    plot_exp(1) = scatter(grads_exp,... % narrow trans. range
        freqs_exp(:,4)/AFFT.plasmafreq_GHz,'s','filled');
    plot_exp(2) = scatter(grads_exp,... % wide trans. range
        freqs_exp(:,3)/AFFT.plasmafreq_GHz,'d','filled');
    
    ylim(freq_lims/AFFT.plasmafreq_GHz)
    
    ylabel('freq. / plasma freq.')
    
    legend([plot_sim(1) plot_sim(2) plot_exp(1) plot_exp(2)],...
        'Sim. narrow','Sim. wide','Exp. narrow','Exp. wide',...
        'location','southeast','AutoUpdate','off');
    
    ax_handle = gca;
    plot([-2.1 2.1],[1,1],'b','LineWidth',2)
    plot(grads_sim,[0.894394264356599,0.921920505722216,0.948648374071357,0.974643553503204,0.999963186893636,1.02465735438711,1.04877023825847,1.07234105181787,1.09540478827280],'b','LineWidth',2)
    
else
    scatter(trans_lims_m,...
        peak_freqs_plot,dotsize,'filled');
    %         ylim([min(peak_freqs_plot)*0.99,max(peak_freqs_plot)*1.01])
    ylim(freq_lims);
    ylabel('freq. (GHz)')
end % if normalized frequency switch

%     ax_handle = gca;
%     plot(3*0.2*sqrt(1+100/4.9^2)*ones(1,2),ax_handle.YLim,'--r','LineWidth',2) % dashed line at 3 sigma_r (bunch travelling in vacuum)
%     plot(ax_handle.XLim,freq_widerange*ones(1,2),'--b','LineWidth',2);

xlim([-2.1 2.1])
hold off
grid on
xlabel('gradient (%/m)')

P = Plotty('fig_handle',fig33,'plots_dir',plots_dir,'plot_name',['fvsg','','']);

if saveplots
    P.save_plot();
end



