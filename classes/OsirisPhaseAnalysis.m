%__________________________________________________________________________
% Subclass of OsirisDenormalizer that does the matrix for the waterfall
% plot from trans. or long. wakfields, as well as from proton bunch density
% data (microbunches). It also calculates the dephasing close to some xi by
% using the zero-crossing or maximum value of each row of the waterfall
% matrix. A 3° polynomial fit is done to the region of interest, to
% counteract the effec of noise.
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
% Methods
% - yyy: (function description)
%
%


classdef OsirisPhaseAnalysis < handle & OsirisDenormalizer
    
    properties(Constant, Hidden)
        
        
        
        
    end % constant properties
    
    properties
        
        % input, description in parser
        force_waterfall;
        dump_list;
        dephasing_xi;
        dephasing_search;
        xi_list;
        mean_weights_flag;
        
        % output
        waterfall_mat; % output example
        waterfall_xi;
        waterfall_z;
        dephasing_line;
        dephasing_first;
        weight_line;
        mean_dephasing;
        std_dephasing;
        
    end % properties
    
    properties (Hidden)
        
        waterfall_filename; % file name for the waterfall matrix
        density_weights_filename; % file name for the density weights for average of dephasing
        % --------- out:
        
        % count_microbunches
        microbunch_peaks; % peaks of the microbunches for a selected dump
        microbunch_locs; % locations of the peaks of the microbunches for a selected dump
        microbunch_number; % number of the microbunch close to the selected xi
        microbunch_loc; % position of a selected microbunch
        microbunch_peak; % peak value of the microbunch at position microbunch_loc
        
    end % hidden properties
    
    methods
        
        % --CONSTRUCTOR--
        
        function obj = OsirisPhaseAnalysis(varargin)
            
            % obj.example_function_in_constructor()
            
            % parse input to load files
            p = inputParser;
            
            % comments (describing variable)
            p.addParameter('any_text', 'baseline', @(x) ischar(x));
            % calculate the depahsing looking for the zero crossing or the
            % maximum/minimum value of the fields (minimum not added yet)
            p.addParameter('dephasing_search', '0x', @(x) ismember(x,{'0x','max'}));
            % comments
            p.addParameter('force_waterfall', false, @(x) islogical(x));
            % comments
            p.addParameter('dephasing_xi', 7, @(x) isfloat(x) && x > 0);
            p.addParameter('dump_list', 1:100, @(x) isfloat(x) && all(x >= 0));
            
            % flag if the weights are needed for the average phase
            % calculation (mean is the charge around the axis close to that
            % dephasing value)
            p.addParameter('mean_weights_flag', false, @(x) islogical(x));
            
            % comments
            p.addParameter('xi_list', [], @(x) isfloat(x) && all(x >=0));
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            if isempty(fieldnames(p.Unmatched))
                unmatched = {};
            else
                unmatched = reshape([fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            
            obj = obj@OsirisDenormalizer(unmatched{:}); %Parse to superclass OsirisDenormalizer.m
            
            obj.force_waterfall             = p.Results.force_waterfall;
            obj.dephasing_xi                = p.Results.dephasing_xi;
            obj.dephasing_search            = p.Results.dephasing_search;
            obj.dump_list                   = p.Results.dump_list;
            obj.xi_list                     = p.Results.xi_list;
            obj.mean_weights_flag           = p.Results.mean_weights_flag;
            
        end % constructor
        
        
        function obj = dephasing(obj)
            
            % out: dephasing_line
            
            obj.build_waterfall(); % set obj.waterfall_mat
            
            dephasing_point = obj.simulation_window - obj.dephasing_xi;
            waterfall_dephasing = flipud(obj.waterfall_mat);
            if obj.mean_weights_flag
                save_property = obj.property;
                save_trans_range = obj.trans_range;
                obj.property = 'density';
                obj.trans_range = [0 0.02];
                obj.build_waterfall();
                weights_mat = flipud(obj.waterfall_mat);
                obj.property = save_property;
                obj.trans_range = save_trans_range;
            end
            obj.weight_line = zeros(1,size(waterfall_dephasing,1));
            for ii = 1:size(waterfall_dephasing,1)
                
                dephasing_range_xi = dephasing_point + [-obj.plasma_wavelength/3,obj.plasma_wavelength/3];
                
                ind_xi = (obj.waterfall_xi > dephasing_range_xi(1)) & (obj.waterfall_xi < dephasing_range_xi(2));
                poly_xi = obj.waterfall_xi(ind_xi);
                field_line = waterfall_dephasing(ii,ind_xi);
                
                if obj.mean_weights_flag
                    dephasing_range_xi = dephasing_point + [-obj.plasma_wavelength/2,obj.plasma_wavelength/2];
                    ind_xi = (obj.waterfall_xi > dephasing_range_xi(1)) & (obj.waterfall_xi < dephasing_range_xi(2));
                    obj.weight_line(ii) = sum(weights_mat(ii,ind_xi));
                end
                
                [pfit,~,mu]  = polyfit(poly_xi,field_line,4);
                
                poly_xi10 = linspace(poly_xi(1),poly_xi(end),10*length(poly_xi));
                poly_data = polyval(pfit,poly_xi10,[],mu);
                
                switch obj.dephasing_search
                    case '0x'
                        up_zerox = @(v) find(v(1:end-1) <= 0 & v(2:end) > 0);	% Returns Up Zero-Crossing Indices
                        down_zerox = @(v) find(v(1:end-1) >= 0 & v(2:end) < 0);    % Returns Down Zero-Crossing Indices
                        fun_zerox = @(x0,y0,x1,y1) x0 - (y0.*(x0 - x1))./(y0 - y1); % Interpolated x value for Zero-Crossing
                        ind_zerox = sort([up_zerox(poly_data),down_zerox(poly_data)]);
                        zero_x = fun_zerox(poly_xi10(ind_zerox),poly_data(ind_zerox),poly_xi10(ind_zerox+1),poly_data(ind_zerox+1));
                        
                        % Checking for zero at the ignored value
                        if poly_data(end)==0
                            zero_x(end+1) = poly_xi10(end);
                        end
                        dephasing_line_temp = zero_x;
                        
                    case 'max'
                        [~,max_pos] = max(poly_data);
                        dephasing_line_temp = poly_xi10(max_pos);
                        
                end
%                 if ii > 89 && ii < 96
%                     figure(10)
%                 hold on
%                 plot(poly_xi,field_line)
%                 hold off
%                 figure(11)
%                 hold on
%                 plot(poly_xi10,poly_data)
%                 hold off
%                 end
                % Ideally there should be one zero crossing, specially
                % since an interpolation of the fields is used to find it,
                % but in the case there are many (range to large, or data too noisy),
                % I choose the one closer to xi.
                if length(dephasing_line_temp) == 1
                    obj.dephasing_line(ii) = dephasing_line_temp;
                elseif length(dephasing_line_temp) <= 5 && ~isempty(dephasing_line_temp)
                    [~,close_to_xi] = min(abs(dephasing_line_temp - dephasing_point));
                    obj.dephasing_line(ii) = dephasing_line_temp(close_to_xi);
                else
                    obj.dephasing_line(ii) = nan;
                end % if lenght zero x
                
                if ii >= 2
                    dephasing_line_nonnan = obj.dephasing_line(1:ii-1);
                    dephasing_line_nonnan(isnan(dephasing_line_nonnan)) = [];
                    if ~isempty(dephasing_line_nonnan) && (obj.dephasing_line(ii) - dephasing_line_nonnan(end) >= obj.plasma_wavelength/2)
                        obj.dephasing_line(ii) = obj.dephasing_line(ii) - obj.plasma_wavelength/2;
                    end
                end % if ii 2
                
                % change xi to current xi
                if ~isnan(obj.dephasing_line(ii))
                    dephasing_point = obj.dephasing_line(ii);
                end
            end
            
            % start dephasing at 0, and normalize to the plasma wavelength
            % at the beginning of the plasma.
            ind_nonnan = ~isnan(obj.dephasing_line);
            dephasing_line_nonnan = obj.dephasing_line(ind_nonnan);
            obj.dephasing_first = dephasing_line_nonnan(1);
            obj.dephasing_line = (obj.dephasing_line-obj.dephasing_first)/obj.plasma_wavelength;
            
            
        end % dephasing
        
        function obj = dephasing_average(obj)
            obj.dephasing_xi = obj.xi_list(1);
            obj.dephasing();
            dephasing_lines = zeros(length(obj.xi_list),size(obj.waterfall_mat,1));
            dephasing_lines(1,:) = obj.dephasing_line;
            
            if obj.mean_weights_flag
                density_weights = zeros(size(dephasing_lines));
                density_weights(1,:) = obj.weight_line;
            else
                density_weights = ones(size(dephasing_lines));
            end
            
            for nxi = 2:length(obj.xi_list)
                obj.dephasing_xi = obj.xi_list(nxi);
                obj.dephasing();
                dephasing_lines(nxi,:) = obj.dephasing_line;
                
                if obj.mean_weights_flag
                    density_weights(nxi,:) = obj.weight_line;
                end
                obj.progress_dump('dephasing',nxi,length(obj.xi_list));
            end % for nxi xi_list
            
            temp_mean_dephasing = zeros(1,size(obj.waterfall_mat,1));
            temp_std_dephasing = zeros(1,size(obj.waterfall_mat,1));
            
            temp_mean_dephasing(1) = sum(density_weights(:,1).*dephasing_lines(:,1),'omitnan')./sum(density_weights(:,1));
            temp_std_dephasing(1) = std(dephasing_lines(:,1),density_weights(:,1),'omitnan');
            diff_dephasing = diff(dephasing_lines,1,2);
            for w = 2:size(obj.waterfall_mat,1)
                temp_mean_dephasing(w) = temp_mean_dephasing(w-1) + sum(density_weights(:,w).*diff_dephasing(:,w-1),'omitnan')./sum(density_weights(:,w));
                temp_std_dephasing(w) = std(diff_dephasing(:,w-1),density_weights(:,w),'omitnan');
            end
            
            obj.mean_dephasing = temp_mean_dephasing;
            obj.std_dephasing = temp_std_dephasing;
            
            
%             obj.mean_dephasing = sum(density_weights.*dephasing_lines,'omitnan')./sum(density_weights);
%             obj.std_dephasing = zeros(1,size(obj.waterfall_mat,1));
%             for w = 1:size(obj.waterfall_mat,1)
% %                 obj.std_dephasing(w) = sqrt(sum((dephasing_lines(:,w)-obj.mean_dephasing(w)).^2.*density_weights(:,w),'omitnan')/...
% %                     sum(density_weights(:,w),'omitnan'));
%                 obj.std_dephasing(w) = std(dephasing_lines(:,w),density_weights(:,w),'omitnan');
%             end
%             pause
        end % phase_average
        
        function obj = build_waterfall(obj)
            if obj.useAvg
                avg = 'avg';
            else
                avg = '';
            end
            % check if there is already a waterfall plot to load
            obj.waterfall_filename = ['save_files/waterfall/waterfall_mat_',...
                obj.property,'_',obj.datadir,'_',obj.wakefields_direction,...
                num2str(obj.dump_list(1)),'to',num2str(obj.dump_list(end)),...
                avg,'.mat'];
            
            if isfile(obj.waterfall_filename) && ~obj.force_waterfall
                
                load_waterfall = load(obj.waterfall_filename);
                obj.waterfall_mat = load_waterfall.waterfallmat;
                obj.waterfall_z = load_waterfall.waterfallz;
                obj.waterfall_xi = load_waterfall.waterfallxi;
                obj.simulation_window = load_waterfall.simulationwindow;
                
            else
                
                switch obj.property
                    case 'fields'
                        % get the norm. lineout
                        obj.getlineout();
                    case 'density'
                        obj.getdata(); obj.denorm_distance(); obj.assign_density();
                        ir = (obj.r >= obj.trans_range(1)) & (obj.r < obj.trans_range(2));
                        obj.nlineout = 2*pi*obj.dr*sum((obj.r(ir)').*obj.(['n',obj.species])(ir,:),1);
                end % switch property
                
                % initialize variables
                skip = 0;
                propagation_distance_ini = [];
                obj.waterfall_mat = zeros(length(obj.dump_list),length(obj.nlineout));
                
                for n = 1:length(obj.dump_list)
                    
                    obj.dump = obj.dump_list(n);
                    obj.denorm_distance();
                    
                    switch obj.property
                        case 'fields'
                            obj.getlineout();
                        case 'density'
                            obj.getdata(); obj.denorm_distance(); obj.assign_density();
                            obj.nlineout = 2*pi*obj.dr*sum((obj.r(ir)').*obj.(['n',obj.species])(ir,:),1);
                            
                    end
                    
                    % check if first window (with no plasma) is not present
                    % anymore in the simulation window
                    if obj.ntime < obj.n_simulation_window
                        skip = skip + 1;
                        continue
                    elseif isempty(propagation_distance_ini)
                        propagation_distance_ini = obj.propagation_distance; % get initial value of the propagation distance for the plot
                    end
                    
                    obj.waterfall_mat(n - skip,:) = obj.nlineout;
                    
                    % if somehow the dumps have different lengths, activate the
                    % fix below and comment the line above.
                    
                    %     if length(waterfall_mat(im_ind,:)) == length(lineout)
                    %         waterfall_mat(im_ind,:) = lineout;
                    %     else
                    %         waterfall_mat(im_ind,:) =  waterfall_mat(im_ind+1,:);
                    %     end
                    
                    fprintf('waterfall: %.d / %.d \n',n,length(obj.dump_list));
                    
                    
                end % for dump list
                
                obj.waterfall_mat = fliplr(obj.waterfall_mat(1:n-skip,:));
                obj.waterfall_xi = linspace(obj.z(1)-obj.dtime,obj.z(end)-obj.dtime,size(obj.waterfall_mat,2));
                
                obj.waterfall_z = [propagation_distance_ini obj.propagation_distance]/100;
                
                switch obj.property
                    case 'fields'
                        obj.waterfall_mat = obj.denorm_Efield(obj.waterfall_mat);
                    case 'density'
                        obj.waterfall_mat = obj.denorm_density(obj.waterfall_mat);
                end
                
                obj.waterfall_mat = rot90(obj.waterfall_mat,2);
                waterfallmat = obj.waterfall_mat;
                waterfallxi = obj.waterfall_xi;
                waterfallz = obj.waterfall_z;
                simulationwindow = obj.simulation_window;
                save(obj.waterfall_filename,...
                    'waterfallmat','waterfallxi','waterfallz','simulationwindow');
                
                
                
            end
            
        end % waterfall
        
        function obj = count_microbunches(obj)
            obj.getdata(); obj.denorm_distance(); obj.assign_density();
            ir = (obj.r >= obj.trans_range(1)) & (obj.r < obj.trans_range(2));
            obj.nlineout = 2*pi*obj.dr*sum((obj.r(ir)').*obj.(['n',obj.species])(ir,:),1); % integrated lineout
            obj.denorm_density(); % to denormalize the lineout (proton density on axis)
            
            
            [obj.microbunch_peaks,obj.microbunch_locs] = findpeaks(obj.lineout,obj.waterfall_xi,...
                'MinPeakDistance',obj.plasma_wavelength*0.85);
            
            dephasing_xi_flipped = obj.simulation_window - obj.dephasing_xi; % flipped to match correctly
            
            [~,min_loc] = min(abs(obj.microbunch_locs - dephasing_xi_flipped));
            
            obj.microbunch_number = length(obj.microbunch_locs) - min_loc + 1;
            obj.microbunch_loc = obj.simulation_window - obj.microbunch_locs(min_loc);
            obj.microbunch_peak = obj.microbunch_peaks(min_loc);
            
            
        end % count_microbunches
        
    end % ordinary methods
    
    
    methods(Static)
        
        function xxx
            
        end % xxx
        
    end % static methods
    
end % classdef