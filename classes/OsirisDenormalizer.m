%__________________________________________________________________________
% Subclass of OsirisLoader that contains all the function related to the
% denormalization of variables.
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 04/06/2020
%__________________________________________________________________________

% Input: units, some methods have extra input
%
% Output:
%
% Methods
% - yyy: (function description)
%
%


classdef OsirisDenormalizer < handle & OsirisLoader
    
    properties(Constant,Hidden)
        
        % Physical constants
        
    end % constant properties
    
    properties(Hidden)
        
        % step size
        dr; dz;
        assign_after_trim_flag = true; % flag specifying if assign after trimming is needed (default = true)

    end
    
    properties
        
        % input, description in parser
        units;
        
        % output (all in physical units)
        dataOut; % matrix with the desired data
        time; % ntime of the dump in norm. units
        z; r; % long. and trans. axes
        simulation_window; % simulation window
        propagation_distance; % propagation distance (default in cm)
        propagation_distance_m; % propagation distance in m
        dtime; % time denormalized as distance
        
        % output
        plasmafreq; % plasma angular frequency in rad (omega_p)
        plasmafreq_GHz; %  plasma frequency in GHz
        plasma_wavelength; % plasma "wavelength" in cm
        
        % fields
        longfield; % longitudinal wakefields (Ez)
        transfield; % transverse wakefields (Er-Bth)
        
        % lineout
        lineout; 
        
        % species density
        proton_beam; % proton bunch density
        plasma_electrons; % plasma electrons density
        electron_bunch; % accelerated electron bunch density
        % more species to come ...
        
        % raw_data
        z_raw; % long. pos. of particle
        r_raw; % trans. pos. of particle
        pz_raw; % momentum in z of particle
        pr_raw; % momentum in r of particle
        pth_raw; % momentum in theta of particle
        E_raw; % energy of particle
        
        % data trimming (physical units)
        trans_range; % range for transverse trimming
        xi_range; % range for xi trimming
        
        
    end % properties
    
    methods
        
        % --CONSTRUCTOR--
        
        function obj = OsirisDenormalizer(varargin)
            
            % obj.example_function_in_constructor()
            
            % parse input to load files
            p = inputParser;
            
            % units for the physical units (each parameters has its default unit if this parameter is not specified)
            p.addParameter('units', 'cm', @(x) ismember(x,...
                {'cm','mm','um','m','ps','s',...
                'MV/m','MV/cm','GV/cm','V/cm','V/m'}));
            
            % trimming
            p.addParameter('trans_range',[0 1e10], @(x) isfloat(x));
            p.addParameter('xi_range',[1e10 -1], @(x) isfloat(x));
            
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            if isempty(fieldnames(p.Unmatched))
                unmatched = {};
            else
                unmatched = reshape([fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj = obj@OsirisLoader(unmatched{:}); %Parse to superclass OsirisLoader.m
            
            obj.units             = p.Results.units;
            obj.trans_range       = p.Results.trans_range;
            obj.xi_range          = p.Results.xi_range;
            
            obj.den2freq();
            
        end % constructor
        
        % --denormalizing functions--
        
        % for all denorm and norm functions, if they recieve no argument,
        % varargout becomes the obj itself (as usual). If they recieved an
        % extra argument (varargin) it considers it as the matrix to be
        % (de)normalized, and returns it as first argument. 
        
        function varargout = denorm_distance(obj,varargin)
            
            
            switch obj.units
                case 'um'
                    denorm_factor = obj.c_m*1e6/obj.plasmafreq;
                case 'mm'
                    denorm_factor = obj.c_m*1e3/obj.plasmafreq;
                case 'cm'
                    denorm_factor = obj.c_cm/obj.plasmafreq;
                case 'm'
                    denorm_factor = obj.c_m/obj.plasmafreq;
            end % switch units
            
            if isempty(varargin) % denormalize all distances (also time as distance)    
                obj.r = obj.nr*denorm_factor;
                obj.dr = obj.ndr*denorm_factor;
                obj.z = obj.nz*denorm_factor;
                obj.dz = obj.ndz*denorm_factor;
                obj.dtime = obj.ntime*denorm_factor;
                obj.simulation_window = obj.n_simulation_window*denorm_factor;
                obj.propagation_distance = obj.n_propagation_distance*denorm_factor;
                obj.propagation_distance_m = obj.n_propagation_distance*obj.c_m/obj.plasmafreq;
                varargout{1} = obj;
            else
                varargout{1} = varargin{1}*denorm_factor;
                varargout{2} = obj;
            end % if isempty
            
            
        end % denorm_distance
        
        function varargout = denorm_time(obj,varargin)
            
            switch obj.units
                case 'ps'
                    denorm_factor = 1e12/obj.plasmafreq;
                case 's'
                    denorm_factor = 1/obj.plasmafreq;
                otherwise % default is seconds (so if not set, it would be 'cm' and take this option)
                    denorm_factor = 1/obj.plasmafreq;
            end % switch units
            
            obj.time = obj.ntime*denorm_factor;
            
            if isempty(varargin)
                varargout{1} = obj;
            else
                varargout{1} = varargin{1}*denorm_factor;
                varargout{2} = obj;                
            end
            
        end % denorm_time
        
        function varargout = denorm_density(obj,varargin)
            
            if isempty(varargin)
                switch obj.species
                    case {'proton_beam','proton_bunch'}
                        obj.proton_beam = obj.plasmaden*obj.nproton_beam;
                    case 'plasma_electrons'
                        obj.plasma_electrons = obj.plasmaden*obj.nplasma_electrons;
                    otherwise
                        obj.electron_bunch = obj.plasmaden*obj.nelectron_bunch;
                end % switch species
                varargout{1} = obj;
            else
                varargout{1} = obj.plasmaden*varargin{1};
                varargout{2} = obj; 
            end % if isempty varargin
            
            % if there is some lineout
            if ~isempty(obj.nlineout)
                obj.lineout = obj.plasmaden*obj.nlineout;
            end % is empty lineout
            
        end % denorm_density
        
        function varargout = denorm_Efield(obj,varargin)
            % numerical values from Osiris User's Guide:
            % https://osirisdoc.wimpzilla.ist.utl.pt/dev/index.php/User_Guide
            
            switch obj.units
                
                case 'GV/cm'
                    denorm_factor = (1.704e-14)*obj.plasmafreq;
                case 'MV/cm'
                    denorm_factor = (1.704e-11)*obj.plasmafreq;
                case 'MV/m'
                    denorm_factor = (1.704e-9)*obj.plasmafreq;
                case 'V/cm'
                    denorm_factor = (1.704e-5)*obj.plasmafreq;
                case 'V/m'
                    denorm_factor = (1.704e-3)*obj.plasmafreq;
                otherwise % default is MV/m
                    denorm_factor = (1.704e-9)*obj.plasmafreq;
                    
            end % switch units
            
            obj.transfield = obj.ntransfield*denorm_factor;
            obj.longfield  = obj.nlongfield*denorm_factor;
            
            % if there is some lineout
            if ~isempty(obj.nlineout)
                obj.lineout = obj.nlineout*denorm_factor;
            end % is empty lineout
            
            if isempty(varargin)
                varargout{1} = obj;
            else
                varargout{1} = varargin{1}*denorm_factor;
                varargout{2} = obj;                
            end
            
        end % denorm_Efield
      
        function obj = denorm_Bfield(obj)
            error('NOT CODED (so far it had not been needed)');
        end % denorm_Bfield
         
        
        % --normalizing functions--
        
        function varargout = norm_distance(obj,varargin)
            
            switch obj.units
                case 'um'
                    norm_factor = obj.plasmafreq/(obj.c_m*1e6);
                case 'mm'
                    norm_factor = obj.plasmafreq/(obj.c_m*1e3);
                case 'cm'
                    norm_factor = obj.plasmafreq/obj.c_cm;
                case 'm'
                    norm_factor = obj.plasmafreq/obj.c_m;
            end
            
            if isempty(varargin)
                varargout{1} = obj;
            else
                varargout{1} = varargin{1}*norm_factor;
                varargout{2} = obj;
            end
                        
        end % norm_time
        
        function varargout = norm_time(obj,varargin)
            
            switch obj.units
                case 'ps'
                    denorm_factor = 1e12*obj.plasmafreq;
                case 's'
                    denorm_factor = obj.plasmafreq;
                otherwise % default is seconds (so if not set, it would be 'cm' and take this option)
                    denorm_factor = obj.plasmafreq;
            end % switch units
            
            
            if isempty(varargin)
                varargout{1} = obj;
            else
                varargout{1} = varargin{1}*denorm_factor;
                varargout{2} = obj;
            end
            
        end % norm_time
        
        % --other methods--
        
        function obj = den2freq(obj)
            % Small formula to calculate the plasma frequency from the electron
            % plasma density, for 'cold' electrons (ignored thermal motion).
            %
            % Protons considered static.
            %
            % Last update: 19/08/2019
            
            % give warning if plasmaden = 1 (so not set)
            if obj.plasmaden == 1
                warning('Plasma density set to 1/cm^3.')
            end
            
            % change density to m^-3
            plasmaden = obj.plasmaden*1e6; % m
            
            obj.plasmafreq = sqrt(plasmaden*(obj.e_charge_C^2)/(obj.e_mass_kg*obj.permittivity));
            obj.plasmafreq_GHz = obj.plasmafreq/(2*pi*1e9);
            
            obj.plasma_wavelength = 2*pi*obj.c_cm./obj.plasmafreq; % cm
            
        end %den2freq
        
        function obj = trim_data(obj)
            % trim the data in normalized units
            % don't forget the assign the data afterwards
            
            obj.denorm_distance(); % make sure to have normalized axis
            
            z_ind = obj.z > obj.dtime+obj.simulation_window - obj.xi_range(1) & ... %large
                obj.z <= obj.dtime+obj.simulation_window - obj.xi_range(2); % small
            obj.z = obj.z(z_ind);
            obj.nz = obj.nz(z_ind);
            
            r_ind = obj.r >= obj.trans_range(1) & ...
                obj.r < obj.trans_range(2);
            obj.r = obj.r(r_ind);
            obj.nr = obj.nr(r_ind);
                        
            % obj.dataOut(r_ind,z_ind);
            obj.ndataOut = obj.ndataOut(r_ind,z_ind);
            
            % load fields or density to corresponding variable
            if obj.assign_after_trim_flag
                switch obj.property
                    case 'fields'
                        obj.assign_fields();
                    case 'density'
                        obj.assign_density();
                end % switch property
            end % if assign after trim
        end % trim_data
        
        function obj = trim_rawdata(obj)
            % trim the data in normalized units
            % don't forget the assign the data afterwards
            
            obj.denorm_distance();
                        
            z_ind = obj.z_raw > max(obj.z) - obj.xi_range(1) & ... %large
                obj.z_raw <= max(obj.z) - obj.xi_range(2); % small
            
            r_ind = obj.r_raw >= obj.trans_range(1) & ...
                obj.r_raw < obj.trans_range(2);
            
            rz_ind = r_ind & z_ind;
            
            obj.r_raw = obj.r_raw(rz_ind);
            obj.z_raw = obj.z_raw(rz_ind);
            obj.q_raw = obj.q_raw(rz_ind);
            
        end % trim_rawdata
        
        function obj = getlineout(obj)
            
            obj.getdata();
            obj.assign_after_trim_flag = false;
            obj.trim_data();
            
            % AXIS
            switch obj.lineout_direction
                case 'long'
                    obj.nlineout = obj.ndataOut(obj.lineout_point,:);
                case 'trans'
                    obj.nlineout = obj.ndataOut(:,obj.lineout_point);
            end % switch lo dir
            
        end % getlineout
        
    end % ordinary methods
    
    
    methods(Static)
        
        
        % function out = name(in)
        % end
        
    end % static methods
    
end % classdef