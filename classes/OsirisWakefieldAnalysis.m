%__________________________________________________________________________
% Subclass of OsirisDenormalizer does various analyses on the amplitude of
% the wakefields.
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 14/07/2020
%__________________________________________________________________________

% Input:
%
% Output:
%
% Methods
% - yyy: (function description)
%
%


classdef OsirisWakefieldAnalysis < handle & OsirisDenormalizer
    
    properties(Constant, Hidden)
        
        
        
        
    end % constant properties
    
    properties
        
        sigma_z;
        bunch_center;
        search_type;
        dump_list;
        
    end % properties
    
    properties (Hidden)
        
        % output
        amplitude_z;
        pos_amplitude_z;
        propagation_z;
        
    end % hidden properties
    
    methods
        
        % --CONSTRUCTOR--
        
        function obj = OsirisWakefieldAnalysis(varargin)
            
            % obj.example_function_in_constructor()
            
            % parse input to load files
            p = inputParser;
            
            % comments
            p.addParameter('sigma_z', 7, @(x) isfloat(x) && x > 0);
            p.addParameter('bunch_center', 7, @(x) isfloat(x) && x > 0);
            p.addParameter('dump_list', 1:100, @(x) isfloat(x) && all(x >= 0));
            p.addParameter('search_type', 'max', @(x) ismember(x,{'max','mean'}));
            
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            if isempty(fieldnames(p.Unmatched))
                unmatched = {};
            else
                unmatched = reshape([fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            
            obj = obj@OsirisDenormalizer(unmatched{:}); %Parse to superclass OsirisDenormalizer.m
            
            obj.sigma_z                     = p.Results.sigma_z;
            obj.bunch_center                = p.Results.bunch_center;
            obj.dump_list                   = p.Results.dump_list;
            obj.search_type                 = p.Results.search_type;
            
            
            
        end % constructor
        
        
        function obj = amplitude_vs_z(obj)
            
            % allocate memory
            obj.propagation_z = zeros(1,length(obj.dump_list));
            obj.amplitude_z = zeros(1,length(obj.dump_list));
            obj.pos_amplitude_z = zeros(1,length(obj.dump_list));
            use_weights = false;
            if use_weights
                OPA = OsirisPhaseAnalysis('plasmaden',2,'useAvg',true,'dataformat','h5',...
                    'datadir',obj.datadir,'direction','z','dump_list',0:100);
                OPA.property = 'density';
                OPA.trans_range = [0 0.02];
                OPA.waterfall();
                
                weights_mat = [ones(3,length(OPA.waterfall_mat));flipud(OPA.waterfall_mat)];
            end
            for n = 1:length(obj.dump_list)
                
                obj.dump = obj.dump_list(n);
                obj.getdata();
                obj.assign_after_trim_flag = false; % gets ndataOut after trim without assigning and deleting
                obj.trim_data();
                obj.denorm_distance();
                obj.propagation_z(n) = obj.propagation_distance;
                
                use_envelope = false;
                search_xi = false;
                switch obj.search_type
                    case 'max'
                        % gives only upper end of envelope
                        %                         field_envelope = 0.0 + envelope(obj.ndataOut',round(length(obj.ndataOut)/50),'peak')';
                        if use_envelope
                            field_formax = envelope(obj.ndataOut',round(length(obj.ndataOut)/50),'rms')';
                        else
                            field_formax = obj.ndataOut;
                        end
                        max_temp = max(field_formax);
                        [obj.amplitude_z(n),ind_max] = max(max_temp);
                        obj.pos_amplitude_z(n) = (obj.z(ind_max) - obj.dtime - ...
                            (obj.simulation_window - obj.bunch_center))/obj.sigma_z;
                        if search_xi
                        obj.xi_range = (obj.simulation_window - (obj.z(ind_max) - obj.dtime)) +...
                            0.6*[obj.plasma_wavelength,-obj.plasma_wavelength];
                        end
                        obj.progress_dump('maximum field',n,length(obj.dump_list))
                    case 'mean'
                        wakefield_phase_for_mean = 'defoc'; %foc, defoc, accel, decel
                        
                        switch wakefield_phase_for_mean
                            case 'foc'
                                data_for_mean = -obj.ndataOut(obj.ndataOut<0);
                            case 'defoc' %%%---------------------------------------------------------
                                if use_weights
                                    obj.ndataOut = obj.ndataOut.*weights_mat(obj.dump,1:length(obj.ndataOut))/...
                                        (sum(weights_mat(obj.dump,:))*size(obj.ndataOut,1));
                                end
                                data_for_mean = obj.ndataOut(obj.ndataOut>0);
                            otherwise
                                data_for_mean = abs(obj.ndataOut);
                        end
                        
                        obj.amplitude_z(n) = mean(data_for_mean,'All');
                        
                        obj.progress_dump('mean field',n,length(obj.dump_list))
                end % switch search type
                
            end % for dump list

        end %amplitude_vs_z
        
    end % ordinary methods
    
    
    methods(Static)
        
        function xxx
            
        end % xxx
        
    end % static methods
    
end % classdef