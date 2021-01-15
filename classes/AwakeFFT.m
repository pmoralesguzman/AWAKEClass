%__________________________________________________________________________
% Subclass of OsirisLoader which performs FFT on data with different
% options.
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

% Input:
%
% Output:
%
% Methods:


classdef AwakeFFT < handle & OsirisDenormalizer
    
    properties
        
        % input
        on_axis;
        scan_type;
        trans_lims;
        xi_lims;
        n_zeropadding;
        % parameter to find the next power of 2, for n times the length of the time
        % domain data, for the zero padding. It makes finding the peak
        % of the power spectrum more accurate.
        
        % output
        % fft_dataload
        fft_densitymatrix;
        fft_fieldmatrix;
        % get_ftt
        fft_frequencies;
        fft_powerspectrum_den;
        fft_powerspectrum_fld;
        fft_phase_den;
        fft_phase_fld;
        % fft peaks
        fft_low_lim;
        fft_upp_lim;
        peak_distance;
        maxpeak;
        maxloc;
        maxphase;
        
        
    end % standard properties
    
    properties(Hidden)
        
        ind_pksearch;
        phase_starting_point;
        
        locs;
        pks;
        
        center_around_0 = 0;
        
    end % hidden properties
    
    
    methods
        
        function obj = AwakeFFT(varargin)
            
            p = inputParser;
            
            % selection options
            p.addParameter('on_axis','int', @(x) ismember(x,{'int','intw','sum','lineout'}));
            p.addParameter('scan_type','cumulative', @(x) ismember(x,{'cumulative','slice','slice_abs'}));
            p.addParameter('trans_lims',1e10, @(x) isfloat(x));
            p.addParameter('xi_lims',1e10, @(x) isfloat(x));
            p.addParameter('n_zeropadding', 5, @(x) isfloat(x) && (x > 1));
            
            % fft peaks (peak selection)
            p.addParameter('fft_low_lim', 0.5, @(x) isfloat(x) && (x > 0));
            p.addParameter('fft_upp_lim', 1.5, @(x) isfloat(x) && (x > 0));
            p.addParameter('peak_distance', 0.005, @(x) isfloat(x) && (x > 0));
            
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            if isempty(fieldnames(p.Unmatched))
                unmatched = {};
            else
                unmatched = reshape([fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj = obj@OsirisDenormalizer(unmatched{:}); %Parse to superclass OsirisDenormalizer.m
            
            obj.on_axis        = p.Results.on_axis;
            obj.scan_type      = p.Results.scan_type;
            obj.trans_lims     = p.Results.trans_lims;
            obj.xi_lims        = p.Results.xi_lims;
            obj.n_zeropadding  = p.Results.n_zeropadding;
            obj.fft_low_lim    = p.Results.fft_low_lim;
            obj.fft_upp_lim    = p.Results.fft_upp_lim;
            obj.peak_distance  = p.Results.peak_distance;
            
        end % constructor
        
        function obj = fft_dataload(obj,load_switch)
            
            if load_switch
                %load the data according to input
                obj.getdata();
                
                % trim the data
                obj.trim_data(); % zr axes already denorm; data is automatically assigned
                
                obj.denorm_Efield();
                obj.denorm_density();
            end
            
            switch obj.scan_type
                case {'slice','slice_abs'}
                    if obj.trans_lims(1) > 0
                        trans_lims_slice = [0,obj.trans_lims];
                    elseif obj.trans_lims(1) <= 0
                        trans_lims_slice = obj.trans_lims;
                    end
                    r_max = length(trans_lims_slice) - 1;
                case 'cumulative'
                    r_max = length(obj.trans_lims);
            end
            
            obj.fft_densitymatrix = zeros(r_max,length(obj.z));
            obj.fft_fieldmatrix = zeros(length(obj.trans_lims),length(obj.z));
            
            for r = 1:r_max
                
                switch obj.scan_type
                    case 'cumulative'
                        ir = abs(obj.r) < obj.trans_lims(r);
                    case 'slice'
                        ir = (obj.r >= trans_lims_slice(r)) & (obj.r < trans_lims_slice(r+1));
                    case 'slice_abs'
                        ir = (abs(obj.r) >= trans_lims_slice(r)) & (abs(obj.r) < trans_lims_slice(r+1));
                end
                                
                switch obj.property
                    case 'density'
                        switch obj.on_axis
                            case 'int'
                                obj.fft_densitymatrix(r,:) = 2*pi*obj.dr*sum((obj.r(ir)').*obj.(obj.species)(ir,:),1);
                                %                                 obj.fft_densitymatrix(r,:) = dr*sum((4*atan(0.004./obj.r(ir)').*obj.r(ir)').*obj.(obj.species)(ir,:),1);
                            case 'intw' % transform to sum ?
                                obj.fft_densitymatrix(r,:) = 2*pi*trapz(obj.r(ir),(obj.r(ir)').*obj.(obj.species)(ir,:))./obj.r(ir);
                            case 'sum'
                                obj.fft_densitymatrix(r,:) = obj.dr*sum(obj.(obj.species)(ir,:));
                            case 'lineout'
                                ir = find(ir);
                                if ir(end)<5
                                    ir(end) = 5;
                                    warning('Changed trans. limit to cell 5 to avoid noise near axis')
                                end
                                obj.fft_densitymatrix(r,:) = obj.(obj.species)(ir(end),:);
                        end % switch on axis
                        
                        
                    case 'fields'
                        if (~strcmp(obj.on_axis,'lineout') && (r == 1)); warning('Picking lineout'); end
                        ir = find(ir);
                        if ir(end) < 5
                            ir(end) = 5;
                            warning('Changed trans. limit to 5 to avoid noise near axis')
                        end
                        switch obj.wakefields_direction
                            case 'trans'
                                obj.fft_fieldmatrix(r,:) = obj.transfield(ir(end),:);
                            case 'long'
                                obj.fft_fieldmatrix(r,:) = obj.longfield(ir(end),:);
                        end % switch direction
                end % switch property
                
            end % for trans lims
            
            if obj.center_around_0 == 1
                fft_densitymatrix_smooth = smoothdata(obj.fft_densitymatrix,2,'gaussian',101);
                aa = obj.fft_densitymatrix;
                                obj.fft_densitymatrix = obj.fft_densitymatrix - fft_densitymatrix_smooth;
%                 bb = obj.fft_densitymatrix - fft_densitymatrix_smooth;
                a = 1;
            end % center_around_0
            
        end % fft_dataload
        
        function obj = fft_rawdataload(obj)
            
            
            obj.property = 'raw';
            obj.raw_dataset = 'q'; obj.getdata(); obj.assign_raw();
            obj.raw_dataset = 'ene'; obj.getdata(); obj.assign_raw();
            obj.direction = 'r';
            obj.raw_dataset = 'x'; obj.getdata(); obj.assign_raw();
            obj.raw_dataset = 'p'; obj.getdata(); obj.assign_raw();
            obj.direction = 'z';
            obj.raw_dataset = 'x'; obj.getdata(); obj.assign_raw();
            
            obj.property = 'density';
            obj.getdata();
            
%             % denormalize z
            obj.z_raw = obj.denorm_distance(obj.nz_raw);
%             
%             % trim in z 
%             iz = (obj.z_raw >= obj.xi_range(2)) & (obj.z_raw < obj.xi_range(1));
            

            
            
            % push charges and denormalize
            distance = 350;
            new_r = obj.charge_pusher(obj,distance);
            obj.r_raw = obj.denorm_distance(new_r);
           
            % trim the data
            obj.trim_rawdata(); % zr axes already denorm; data is automatically assigned
            obj.trim_data();
                        
            obj.fft_densitymatrix = zeros(length(obj.trans_lims),length(obj.z));
            obj.fft_fieldmatrix = zeros(length(obj.trans_lims),length(obj.z));
            
            for r = 1:length(obj.trans_lims)
                
                switch obj.scan_type
                    case 'cumulative'
                        ir = obj.r_raw < obj.trans_lims(r);
                    case 'slice'
                        if (r == 1) %&& (obj.trans_lims(1)~=0)
                            trans_lims_slice = [obj.trans_lims];
                            %                         elseif obj.trans_lims(1) == 0
                            %                         trans_lims_slice = obj.trans_lims;
                        end
                        ir = (obj.r_raw >= trans_lims_slice(r)) & (obj.r_raw < trans_lims_slice(r+1));
                end
                
                % select only those particles inside the transverse limits
                q_r = obj.q_raw(ir);
                
                % bin
                ind_bin = ceil((obj.z_raw(ir)-obj.dtime)/obj.dz);
                
                % sort the indeces of the binning and then sort the
                % elements in q
                [ind_sort,ind_order] = sort(ind_bin); 
                q_r_sort = q_r(ind_order);
                
                % say how many there is in each bin in a cumulative way,
                % to use that as indices
                % for q, to quickly build up the indexed density along z
                ind_sum = [0,0,histcounts(ind_sort,1:1:length(obj.z),'Normalization','cumcount')]; 

                % go through each element in bz and arrange that charge
                % according to the bins, using the indeces from histconts
                for bz = 1:length(obj.z)
                    obj.fft_densitymatrix(r,bz) = sum(q_r_sort(ind_sum(bz)+1:ind_sum(bz+1)));  
                end
             
            end % for trans lims
            
        end % fft_rawdataload

        
        function obj = get_fft(obj)
            %% Prepare variables for FFT in Matlab
            
            L = length(obj.z); % length of signal
            nzero = 2^nextpow2(obj.n_zeropadding*L); % nearest power of 2 for an efficient implementation of FFT
            
            switch obj.property
                case 'density'
                    fft_input = obj.fft_densitymatrix;
                case 'fields'
                    fft_input = obj.fft_fieldmatrix;
            end % switch property
            
            % Fast Fourier Transform (read Matlab documentation)
            LO = fft(fft_input,nzero,2);
            P2 = abs(LO/nzero);
            
            %  Matlab takes the positive and negative frequencies, so we take just
            %  the first half, which are the positive ones
            powerspectrum = P2(:,1:nzero/2+1);
            powerspectrum(:,2:end-1) = 2*powerspectrum(:,2:end-1);
            
            fft_phase = zeros(size(powerspectrum));
            for ph = 1:size(LO,1)
                LO_phase = LO(ph,1:nzero/2+1)/nzero;
                tol = max(abs(LO_phase))/100;
                ind_tol = abs(LO_phase) < tol;
                LO_phase(ind_tol) = 0;
                fft_phase(ph,:) = (angle(LO_phase));
            end
            
            
            switch obj.property
                case 'density'
                    obj.fft_powerspectrum_den = powerspectrum;
                    obj.fft_phase_den = fft_phase;
                case 'fields'
                    obj.fft_powerspectrum_fld = powerspectrum;
                    obj.fft_phase_fld = fft_phase;
            end % switch property
            
            % Transform density oscillations from space domain to time domain,
            % considering protons traveling at the speed of light.
            % Sampling frequency would be then how many points in the time ellapsed
            % in one window (already trimmed).
            % Only valid when proton bunch is not evolving rapidly.
            if ~strcmp(obj.units,'cm')
                warning('Units not cm. get_fft works correctly only for cm. Double check.')
            end % if units
            
            Fs = L/(max(obj.z)-min(obj.z))*obj.c_cm; %
            obj.fft_frequencies = Fs*(0:(nzero/2))/nzero;
            
            
        end % get_fft
        
        function obj = fft_peaks(obj,varargin)
            if length(varargin) == 1
                amplitude = varargin{1};
            else
                amplitude = varargin{1};
                fft_phase = varargin{2};
                long_profile = varargin{3};
            end
            % output (maxloc) is given already in GHz
            peak_dis = obj.peak_distance*obj.plasmafreq_GHz;
            low_lim = obj.fft_low_lim*obj.plasmafreq_GHz;
            upp_lim = obj.fft_upp_lim*obj.plasmafreq_GHz;
            
            % transform freqs from Hz to GHz
            freqs = obj.fft_frequencies/1e9;
            
            %% Find peaks
            % look for peaks within the set limits
            obj.ind_pksearch = ((freqs > low_lim) & (freqs < upp_lim));
            
            [obj.pks,obj.locs] = findpeaks(amplitude(obj.ind_pksearch),freqs(obj.ind_pksearch),...
                'MinPeakDistance',peak_dis);
            [obj.maxpeak,maxlocs] = max(obj.pks);
            if isempty(obj.maxloc) || ~(obj.maxloc == 115)  % some harcoded thing for tracking, no meaning
                obj.maxloc = obj.locs(maxlocs);
            end
            
            % calculate phase of dominant frequency
            
            if isempty(obj.maxloc)
                obj.maxphase = nan;
            elseif length(varargin) == 3
                %                 ind_maxloc = obj.maxloc == freqs;
                %                 switch obj.property
                %                     case 'density'
                %                         obj.maxphase = obj.fft_phase_den(ind_maxloc);
                %                     case 'fields'
                %                         obj.maxphase = obj.fft_phase_fld(ind_maxloc);
                %                 end % switch property
                
                ind_maxloc = obj.maxloc == freqs;
                
                obj.phase_starting_point = fft_phase(ind_maxloc);
                fo = fitoptions('Method','NonLinearLeastSquares',...
                    'Lower',[0],...
                    'Upper',[2*pi],...
                    'StartPoint',[obj.phase_starting_point]);
                wavenumber = obj.maxloc*1e9/obj.c_cm;
                norm_amplitude = max(amplitude)/2;
                ft = fittype('(cos(wavenumber*(2*pi*z) + phi) > 0)*cos(wavenumber*(2*pi*z) + phi)','problem','wavenumber',...
                    'independent','z','dependent','density','options',fo);
                fit_results = fit(obj.z',long_profile',ft,'problem',wavenumber);
                obj.maxphase = fit_results.phi;
            end % if isempty maxloc
            
            
        end % fft_peaks
        
    end % standard methods
    
    
end % AwakeFFT