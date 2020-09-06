%________________________________________________________________________
%Plot lineout
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
datadirs = {'g0'};
plots_dir = 'fft/rz/grads/g0';
P = Plotty();


% parameters
plasma_density = 1.81e14;

% properties
property = 'fields';
species = 'proton_beam';
field = 'e';
direction = 'z';

% simulation parameters
dump_list = 1:1:1;
dataformat = 'h5';
useAvg = true;

% transverse limits
trans_lims = [0.02];

% analysis parameters
scan_type = 'slice'; % slice, cumulative
on_axis = 'int'; % int, sum, intw, lineout

% switches
plot_powerspectra = true;
save_all_plots = true;

% calculated parameters

% initialize variables
prop_distance_m = zeros(1,length(dump_list));

AFFT = AwakeFFT(...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'direction',direction,'wakefields_direction','long',...
    'dataformat',dataformat,'useAvg',useAvg,...
    'trans_lims',trans_lims,...
    'scan_type',scan_type,'on_axis',on_axis);

for d = 1:length(datadirs)
    AFFT.datadir = datadirs{d};
    
    for n = 1:length(dump_list)
        AFFT.dump = dump_list(n);
        
        
        AFFT.fft_dataload();
        prop_distance_m(n) = AFFT.propagation_distance/100; % propagation distance in m

        
        switch AFFT.property
            case 'density'
                charge_in_slice(:,n) = AFFT.dz*sum(AFFT.fft_densitymatrix,2);
            case 'fields'
                field_lineout(:,n) = AFFT.fft_fieldmatrix;
                data_in_slice = AFFT.fft_powerspectrum_fld;
        end % switch property
        
        for r = 1:1
                        
            %%% Plots of the power spectrum.
            if save_all_plots
                trans_lims_plot_title = [0,trans_lims];
                
                indz = AFFT.z - AFFT.dtime > 7;
                z_plot = AFFT.z(indz);



                field_plot = AFFT.fft_fieldmatrix(r,indz);
                P.plots_dir = [plots_dir,'/field_lineouts'];
                P.plot_name = ['field_long_dump',num2str(AFFT.dump),...
                    '_slice',num2str(r)];
                
                figy = figure('visible','on');
                plot(z_plot,field_plot,'LineWidth',2);
                title(['lineout of long. fields (prop. dist. = ',...
                    num2str(AFFT.propagation_distance/100,2),'m; r = ',...
                    num2str(trans_lims_plot_title(r+1)*10,2),' mm)']);
                ylabel('E_z (MV/m)');
                xlabel('z-ct (m)');
                xlim([z_plot(1),z_plot(end)]);
                figy.Units = 'normalized';
                figy.Position = [0.05 0.2 0.9 0.5];
                drawnow;
                P.fig_handle = figy;
                P.save_plot();
%                 close all;
                AFFT.progress_dump('save slice field',r,1);
                
%                 density_plot = AFFT.fft_densitymatrix(r,indz);
%                 P.plots_dir = [plots_dir,'/density_distributions'];
%                 P.plot_name = ['density_long_dump',num2str(AFFT.dump),...
%                     '_slice',num2str(r)];
%                 
%                 figy = figure('visible','on');
%                 plot(z_plot,density_plot,'LineWidth',2);
%                 title(['long. den. distributon (prop. dist. = ',...
%                     num2str(AFFT.propagation_distance/100,2),'m; r = [',...
%                     num2str(trans_lims_plot_title(r)*10,2),', ',num2str(trans_lims_plot_title(r+1)*10,2),']mm)']);
%                 ylabel('density (1/cm)');
%                 xlabel('z-ct (m)');
%                 xlim([z_plot(1),z_plot(end)]);
%                 figy.Units = 'normalized';
%                 figy.Position = [0.05 0.2 0.9 0.5];
%                 drawnow;
%                 P.fig_handle = figy;
%                 P.save_plot();
% %                 close all;
%                 AFFT.progress_dump('save slice density',r,1); 
            end

            
        end % for nslices
        
        AFFT.progress_dump('dump',n,length(dump_list));
    end %  for dump list

    
end % for datadirs







