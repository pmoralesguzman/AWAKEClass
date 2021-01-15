%________________________________________________________________________
% Freq vs r for 1 gradient.
% Comparison of simulation space vs time profile
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
datadirs = {'gm10'};
% datadirs = {'gm20','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
dataformat  = 'mat';
useAvg      = false;
dump        = 133;

% save directory
plots_dir = ['AAC/fft/r/grads/',num2str(dump),''];

save_format = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;
property        = 'density';
species         = 'proton_beam';
field           = 'e';
direction       = 'z';
wakefields_direction = 'trans';

% analysis
xi_range        = [18 0.74]; % cm
trans_lims     = [0.018]*(1:9);
nslices = length(trans_lims);
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






ppnn = 'pn';

for d = 1:length(datadirs)
    for pn = 1:1
        
        % classes
        AFFT = AwakeFFT(...
            'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
            'plasmaden',plasma_density,'property',property,'species',species,...
            'field',field,'direction',direction,'wakefields_direction',wakefields_direction,...
            'xi_range',xi_range,'trans_lims',trans_lims,...
            'scan_type',scan_type,'on_axis',on_axis);
        
        AFFT2 = AwakeFFT(...
            'datadir',datadirs{1},'dataformat',dataformat,'useAvg',useAvg,'dump',dump,...
            'plasmaden',plasma_density,'property',property,'species',species,...
            'field',field,'direction',direction,'wakefields_direction',wakefields_direction,...
            'xi_range',xi_range,'trans_lims',trans_lims,...
            'scan_type',scan_type,'on_axis',on_axis);
        
        
        P = Plotty('plasmaden',plasma_density,'plots_dir',plots_dir,'save_format',save_format,...
            'datadir',datadirs{1});
        
        timetemp = load([datadirs{d},'densitytimeprofile.mat']);
        simtimeprofile = timetemp.densitymatrix;
        z_plot = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2));
        r_plot = linspace(0,timetemp.ylims(2),size(simtimeprofile,1));
        

        
        AFFT2.dr = r_plot(2) - r_plot(1);
        
        %     close;
        datadir = datadirs{d};
        plot_name = ['fvsrsimstudent'];
        
        AFFT.datadir = datadir;
        AFFT.scan_type = scan_type;
        
        if use_raw
            AFFT.fft_rawdataload();
        else
            AFFT.fft_dataload(true);
            AFFT2.r = r_plot;
            % trim
            z = z_plot*1e-12*P.c_cm;
            
            z_ind = z > xi_range(2) & ... %large
                z <= xi_range(1); % small
            AFFT2.z = z(z_ind);
            
            AFFT2.proton_beam = simtimeprofile(:,z_ind);
            AFFT2.fft_dataload(false);
        end
        prop_distance_m = AFFT.propagation_distance/100; % propagation distance in m
        
        % calculates fft and gives AFFT.fft_frequencies
        % and AFFT.fft_powerspectrum (_den or _fld)
        AFFT.get_fft(); AFFT2.get_fft();
        
        switch AFFT.property
            case 'density'
                charge_in_slice = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
                data_in_slice = AFFT.fft_powerspectrum_den;
                data_in_slice2 = AFFT2.fft_powerspectrum_den;
            case 'fields'
                data_in_slice = AFFT.fft_powerspectrum_fld;
        end % switch property
        
        % initialize variables
        peak_freqs = zeros(nslices,1);
        
        for r = 1:nslices
            
            AFFT.fft_peaks(data_in_slice(r,:));
            AFFT2.fft_peaks(data_in_slice2(r,:));
            % maxloc is the location in freq space of the ampltiude peak in
            % the power spectrum
            if isempty(AFFT.maxloc)
                peak_freqs(r) = nan;
            else
                peak_freqs(r) = AFFT.maxloc;
                peak_freqs2(r) = AFFT2.maxloc;
            end
            
        end % for nslices
        
        % Calculate DFT for whole range
        AFFT.scan_type = 'cumulative';
        AFFT.on_axis = 'int';
        AFFT.trans_lims = trans_upperlimit; %three experimental sigmas at the end of the plasma 3*0.455
        if use_raw
            AFFT.fft_rawdataload();
        else
            AFFT.fft_dataload(true);
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
        
        peak_plot{pn} = peak_freqs2;
        
        %     hold on
        
        %     if normalized_frequency_switch
        %         scatter(trans_lims_mm,...
        %             peak_freqs/AFFT.plasmafreq_GHz,dotsize,'filled');
        %         scatter(trans_lims_mm,peak_freqs2/AFFT.plasmafreq_GHz,dotsize,'x','LineWidth',1); % experiment
        %         plot(xlimits,freq_widerange*ones(1,2)/AFFT.plasmafreq_GHz,'--b','LineWidth',2);
        %         freq_lims = [0.98*min([peak_freqs2';peak_freqs]/AFFT.plasmafreq_GHz),...
        %             1.02*max([peak_freqs2';peak_freqs]/AFFT.plasmafreq_GHz)];
        %         ylim(freq_lims)
        %         ylabel({'transverse slice frequency /','plasma frequency at z = 0 m'})
        %         set(gca,'xdir','reverse')
        %     else
        %         scatter(trans_lims_mm,peak_freqs,dotsize,'filled'); % simulation
        %         scatter(exp_translims,exp_freqs,dotsize,'x','LineWidth',1); % experiment
        %         plot(xlimits,freq_widerange*ones(1,2),'--b','LineWidth',2);
        %         freq_lims = [0.98*min([peak_freqs2';peak_freqs]),1.02*max([peak_freqs2';peak_freqs])];
        %         ylim(peak_freqs);
        %         ylabel('freq. (GHz)')
        %     end % if normalized frequency switch
        %
        %     hold off
        %     grid on
        %     xlim(xlimits)
        %     xlabel('r (mm)')
        %     legend('simulation','experiment','location','best');
        %     drawnow;
        %
        %     P.plot_name = plot_name;
        %     P.fig_handle = fig_fvsr;
        %     P.save_plot();
        
    end %posinega
    
end % for datadirs

fig_fvsr = figure;

xlimits = [0,1.01*trans_lims_mm(end)];

hold on
scatter([trans_lims_mm],...
    [peak_freqs]/AFFT.plasmafreq_GHz,[dotsize],'filled'); %space-dependent
scatter([trans_lims_mm],[peak_plot{1}]/AFFT.plasmafreq_GHz,...
    [dotsize],'x','LineWidth',2); % time-dependent
plot(xlimits,freq_widerange*ones(1,2)/AFFT.plasmafreq_GHz,'--b','LineWidth',2);

hold off

freq_lims = [0.98*min([peak_freqs2';peak_freqs]/AFFT.plasmafreq_GHz),...
    1.02*max([peak_freqs2';peak_freqs]/AFFT.plasmafreq_GHz)];
ylim(freq_lims); xlim(xlimits);
ylabel({'normalized transverse slice frequency (a.u.)'})
xlabel('r (mm)')
legend('space-dependent / c','time-dependent','location','best');
P.plot_name = plot_name;
P.fig_handle = fig_fvsr;
P.save_plot();
