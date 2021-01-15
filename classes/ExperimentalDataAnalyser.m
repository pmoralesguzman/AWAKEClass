%__________________________________________________________________________
% Subclass of OsirisDenormalizer
% Loads experimental data from the streak camera.
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 07/01/2021
%__________________________________________________________________________

% Input:
%
% Output:
%
% Methods
% - yyy: (function description)
%
%


classdef ExperimentalDataAnalyser < handle & ExperimentalDataLoader
    
    
    properties(Constant, Hidden)
        
    end % constant properties
    
    properties
        
    end % properties
    
    properties
        
        % output
        
        SCI;
        % output of gaussian fits
        SCI_IF_px; % ionization front position in pixels
        SCI_center_px; % bunch axis position in pixels
        SCI_center; % bunch axis position in cm
        SCI_width_px; % rms width in pixels
        SCI_width; % rms width in cm
        
        SCI_r; % transverse axis in cm
        SCI_z; % longitudinal axis in cm

        
    end % hidden properties
    
    methods
        
        % --CONSTRUCTOR--
        
        function obj = ExperimentalDataAnalyser(varargin)
            
            % obj.example_function_in_constructor()
            
            % parse input to load files
            p = inputParser;
            
            % comments
            p.addParameter('sigma_z', 7, @(x) isfloat(x) && x > 0);
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            if isempty(fieldnames(p.Unmatched))
                unmatched = {};
            else
                unmatched = reshape([fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            
            obj = obj@ExperimentalDataLoader(unmatched{:}); %Parse to superclass OsirisDenormalizer.m
            
            obj.sigma_z  = p.Results.sigma_z;
            
            
            
            
        end % constructor
        
        function obj = loadSCdata(obj)
            
            % load not aligned data
            obj.loadSCdata_nonaligned();
            %obj.SCI_nonaligned
            
            
            % calculate ionization front pixel
            laser_pulses = sum(obj.SCI_nonaligned(1:120,:));
            x_axis_for_peaks = 1:length(laser_pulses);
            max_peak_amplitude = max(laser_pulses);
            [peak_amplitudes,laser_peak_locs] = findpeaks(laser_pulses,x_axis_for_peaks,...
                'MinPeakDistance',350,'MinPeakHeight',0.1*max_peak_amplitude);
            
            % check if peaks are in the right position (plot laser pulse
            % and scatter peaks)
            plot_peaks_flag = 0;
            if plot_peaks_flag == 1
                figure;
                plot(x_axis_for_peaks,laser_pulses)
                hold on
                scatter(laser_peak_locs,peak_amplitudes)
                hold off
            end
            
            ionization_front_pulse_peak_location = laser_peak_locs(end);
            ionization_front_start = ionization_front_pulse_peak_location - 40;
            ionization_front_end = ionization_front_pulse_peak_location + 40;
            ionization_front_x = ionization_front_start:ionization_front_end;
            ionization_front_pulse = laser_pulses(ionization_front_x);
            
            gauss_eqn = 'a*exp(-((x-b)/c)^2)+d';
            start_points_laser = [peak_amplitudes(end),ionization_front_pulse_peak_location,...
                30,median(laser_pulses)];
            
            ionization_front_fit = fit(ionization_front_x',ionization_front_pulse',...
                gauss_eqn,'Start',start_points_laser);
            
            % check gaussian fit and pulse by plotting
            plot_fit_flag = 0;
            if plot_fit_flag == 1
                figure;
                plot(ionization_front_x,ionization_front_pulse)
                hold on
                plot(ionization_front_x,ionization_front_fit.a*...
                    exp(-((ionization_front_x-ionization_front_fit.b)./...
                    ionization_front_fit.c).^2) + ionization_front_fit.d);
                hold off
            end
             
            % position in px of ionization front
            obj.SCI_IF_px = round(ionization_front_fit.b);
            
            % place ionization front to match t = 0 from below
            obj.SCI = obj.SCI_nonaligned(:,obj.SCI_IF_px-1500+1:obj.SCI_IF_px+136);
            
            % calculate center
            obj.calculate_SCIcenter();
            
            % calculate width
            obj.calculate_SCIwidth();
            
            % calculate first microbunch position
            
            microbunch_lineout = sum(obj.SCI(obj.SCI_center_px+(-obj.SCI_width_px:obj.SCI_width_px),1:1500));
            [~,pos] = findpeaks(-microbunch_lineout(1:1500-5));
            obj.first_microbunch_position_px = 1500 - pos(end);
            obj.first_microbunch_position = obj.first_microbunch_position_px*obj.px2ps;

            % axes
            
            obj.SCI_r = ((1:1:672)-obj.SCI_center_px)*obj.px2cm; %cm
            obj.SCI_z = linspace(1500,-136,1636)*obj.px2ps; %ps
            
        end % loadSCdata
        
        function obj = calculate_SCIcenter(obj)

            % calculate center according to modulated part
            
            % get modulated part
            modulated_bunch = obj.SCI_nonaligned(:,1:obj.SCI_IF_px);
%             modulated_bunch = obj.SCI_nonaligned(:,obj.SCI_IF_px:end);
            modulated_bunch_sum = sum(modulated_bunch,2);
            bunch_center_y = 1:length(modulated_bunch_sum);
            
            start_points_bunch = [max(modulated_bunch_sum),length(modulated_bunch_sum)/2,...
                length(modulated_bunch_sum)/3,modulated_bunch_sum(end)];
            
            gauss_eqn = 'a*exp(-((x-b)/c)^2)+d';
            
            % gaussian fit for center
            bunch_center_fit = fit(bunch_center_y',modulated_bunch_sum,...
                gauss_eqn,'Start',start_points_bunch);
            
            % check gaussian fit and bunch sum by plotting
            plot_fit_flag = 0;
            if plot_fit_flag == 1
                figure;
                plot(bunch_center_y,modulated_bunch_sum)
                hold on
                plot(bunch_center_y,bunch_center_fit.a*...
                    exp(-((bunch_center_y-bunch_center_fit.b)./...
                    bunch_center_fit.c).^2) + bunch_center_fit.d);
                hold off
            end
            
            obj.SCI_center_px = round(bunch_center_fit.b);
%             obj.SCI_center_px = 338; % gm20 = 338;
            obj.SCI_center = obj.SCI_center_px*obj.px2cm;
                            
        end % calculate SCI center
        
        function obj = calculate_SCIwidth(obj)

            % calculate center according to modulated part
            
            % get modulated part
            unmodulated_bunch = obj.SCI_nonaligned(:,obj.SCI_IF_px:end);
            unmodulated_bunch_sum = sum(unmodulated_bunch,2);
            bunch_center_y = 1:length(unmodulated_bunch_sum);
            
            start_points_bunch = [max(unmodulated_bunch_sum),length(unmodulated_bunch_sum)/2,...
                length(unmodulated_bunch_sum)/3,unmodulated_bunch_sum(end)];
            
            gauss_eqn = 'a*exp(-((x-b)/c)^2)+d';
            
            % gaussian fit for center
            bunch_width_fit = fit(bunch_center_y',unmodulated_bunch_sum,...
                gauss_eqn,'Start',start_points_bunch);
            
            % check gaussian fit and bunch sum by plotting
            plot_fit_flag = 0;
            if plot_fit_flag == 1
                figure;
                plot(bunch_center_y,unmodulated_bunch_sum)
                hold on
                plot(bunch_center_y,bunch_width_fit.a*...
                    exp(-((bunch_center_y-bunch_width_fit.b)./...
                    bunch_width_fit.c).^2) + bunch_width_fit.d);
                hold off
            end
            
            obj.SCI_width_px = round(bunch_width_fit.c/sqrt(2));
            obj.SCI_width = obj.SCI_width_px*obj.px2cm;
            
        end % calculate SCI center

        
        
        
        
    end % ordinary methods
    
    
end % classdef
