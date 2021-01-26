%________________________________________________________________________
% Freq vs r for 1 gradient.
% Comparison of simulation and experiment.
% FFT of the proton density distribution.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 11/01/2021
%________________________________________________________________________
clear;

% data directory
datadirs = {'gm20'};
% datadirs = {'gm20','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
dataformat  = 'mat';
useAvg      = false;
dump        = 133;

% save directory
plots_dir = ['gradsim/fft/r/',datadirs{1},'n',num2str(dump),''];

save_format = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';

% analysis
xi_range        = [14 0.38]; % cm % 2nd microbunch end = 0.656 cm, 1st 384
trans_lims_sim     = 0.02*(1:10); % cm 
trans_lims_exp     = 0.02*(-10:10); % cm
nslices_sim = length(trans_lims_sim); % does not include 0
nslices_exp = length(trans_lims_exp) - 1; % includes a 0

trans_upperlimit   = 0.16; % cm

scan_type       = 'slice'; % slice, cumulative
on_axis         = 'sum'; % int, sum, intw, lineout
max_dotsize     = 150;

% switches
save_plot_flag      = false;
showChargeinDotSize = false; % show dot size to reflect charge

% classes
% AwakeFFT for the simulation
AFFT = AwakeFFT('load_data_flag',0,'center_around_0_flag',0,...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims_sim,...
    'scan_type',scan_type,'on_axis',on_axis);

% AwakeFFT for the experimental data
AFFT2 = AwakeFFT('load_data_flag',0,'center_around_0_flag',0,...
    'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'xi_range',xi_range,'trans_lims',trans_lims_exp,...
    'scan_type',scan_type,'on_axis',on_axis,'center_around_0_points',151);

P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'save_format',save_format,...
    'datadir',datadirs{1});
plot_name = ['fvsr',datadirs{1},'r',num2str(trans_upperlimit*100)];

EDA = ExperimentalDataAnalyser('datadir',datadirs{1});


for d = 1:length(datadirs)
    
    datadir = datadirs{d};
    AFFT.datadir = datadir;
    
    AFFT.datadir = datadirs{d};
    
    timetemp = load([datadirs{d},'densitytimeprofile_slit.mat']);
    simtimeprofile = timetemp.densitymatrix;
    z_sim = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2));
    r_sim = linspace(0,timetemp.ylims(2),size(simtimeprofile,1));
    AFFT.r = r_sim;
    zz_sim = z_sim*1e-12*EDA.c_cm;
    
    z_ind = zz_sim > xi_range(2) & ... %large
        zz_sim <= xi_range(1); % small
    AFFT.z = zz_sim(z_ind);
    AFFT.proton_beam = simtimeprofile(:,z_ind);
    AFFT.dr = r_sim(2)-r_sim(1);
    AFFT.dz = zz_sim(2)-zz_sim(1);
    
    % load simulation data and build densitymatrix
    AFFT.fft_dataload();
    xx = AFFT.fft_densitymatrix;
    
    % load experimental data
    EDA.loadSCdata();
    % set dr in AFFT2 to build densitymatrix
    AFFT2.dr = EDA.SCI_r(2) - EDA.SCI_r(1);
    AFFT2.r = EDA.SCI_r;
    % trim
    z_exp = EDA.SCI_z*1e-12*P.c_cm;
    z_ind = z_exp > xi_range(2) & ... %large
        z_exp <= xi_range(1); % small
    AFFT2.z = z_exp(z_ind);
    
    % set experimental data to proton_beam for analysis
    AFFT2.proton_beam = EDA.SCI(:,z_ind);
    
    % build densitymatrix for the experimental data
    AFFT2.fft_dataload();
    
    % calculates fft and gives AFFT.fft_frequencies and AFFT.fft_powerspectrum
    AFFT.get_fft(); AFFT2.get_fft();
    
    charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
    exp_charge_in_slice = sum(AFFT2.fft_densitymatrix,2);
    sim_densitymatrix = AFFT.fft_densitymatrix; % for checking
    exp_densitymatrix = AFFT2.fft_densitymatrix; % for checking
    
    simdata_in_slice = AFFT.fft_powerspectrum_den;
    expdata_in_slice = AFFT2.fft_powerspectrum_den;
    
    % initialize variables
    simpeak_freqs = zeros(nslices_sim,1);
    exppeak_freqs = zeros(nslices_exp,1);
    for r = 1:nslices_sim
        AFFT.fft_peaks(simdata_in_slice(r,:));
        simpeak_freqs(r) = AFFT.max_loc;
    end % for nslices sim
    
    for r = 1:nslices_exp
        AFFT2.fft_peaks(expdata_in_slice(r,:));
        exppeak_freqs(r) = AFFT2.max_loc;
    end % for nslices exp
    
    % Calculate DFT for whole range
    AFFT.scan_type = 'cumulative';
    AFFT.on_axis = 'int';
    AFFT.trans_lims = trans_upperlimit;
    AFFT.fft_dataload();
    AFFT.get_fft();
    AFFT.fft_peaks(AFFT.fft_powerspectrum_den);
    freq_widerange = AFFT.max_loc;
    
    if showChargeinDotSize
        dotsize = max_dotsize*charge_in_slice/max(charge_in_slice,[],'all');
        dotsize(dotsize == 0) = nan;
    else
        dotsize = max_dotsize*ones(length(charge_in_slice),1);
        dotsize(dotsize == 0) = nan;
    end % end if charge dot size
    
end % for datadirs

fig_fvsr = figure();
fig_fvsr.OuterPosition = [336.1111  444.5556  862.2222  409.3333];
ax_fvsr = axes('Parent',fig_fvsr);
ax_fvsr.FontSize = 12;
% plot limits
trans_lims_mm = (10*trans_lims_sim)'; % trans lims in mm
xlimits = [-1.01*trans_lims_mm(end) 1.01*trans_lims_mm(end)];

hold on
scatter([flipud(-trans_lims_mm);trans_lims_mm],...
    [flipud(simpeak_freqs);simpeak_freqs]/AFFT.plasmafreq_GHz,[dotsize;dotsize],'filled');

scatter([flipud(-trans_lims_mm);trans_lims_mm],exppeak_freqs/AFFT.plasmafreq_GHz,...
    [dotsize;dotsize],'x','LineWidth',1); % experiment
plot(xlimits,freq_widerange*ones(1,2)/AFFT.plasmafreq_GHz,'--b','LineWidth',2);

hold off

freq_lims = [0.98*min(simpeak_freqs/AFFT.plasmafreq_GHz),...
    1.02*max(simpeak_freqs/AFFT.plasmafreq_GHz)];
ylim(freq_lims); xlim(xlimits);
% ylim([0.89 1.05]);
ylabel({'transverse slice','normalized frequency (a.u.)'})
xlabel('x (mm)')
legend('simulation','experiment','location','southeast');
P.plot_name = plot_name;
P.fig_handle = fig_fvsr;
P.save_plot();
