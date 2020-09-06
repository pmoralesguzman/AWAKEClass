%__________________________________________________________________________
% Subclass of OsirisDenormalizer that selects and creates tags for particle
% tracking.
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


classdef OsirisParticleTracking < handle & OsirisDenormalizer
    
    properties(Constant, Hidden)
        
        
    end % constant properties
    
    properties
        
        % input, description in parser
        reference_dump;
        dump_list;
        ntag;
        in_taglist;
        
        % output
        selected_tags; % list of tags from particles that meet the conditions
        partrack; % tracks of tags that meet the conditions
        ind_tag;
        
    end % properties
    
    methods
        
        % --CONSTRUCTOR--
        
        function obj = OsirisParticleTracking(varargin)
            
            % obj.example_function_in_constructor()
            
            % parse input to load files
            p = inputParser;
            
            % comments (describing variable)
            p.addParameter('any_text', 'baseline', @(x) ischar(x));
            % comments
            p.addParameter('options_list', 'opt1', @(x) ismember(x,{'opt1','opt2'}));
            % comments
            p.addParameter('boolean', false, @(x) islogical(x));
            % range to look for particles in the reference dump in cm
            p.addParameter('trans_range', [0 1], @(x) isfloat(x));
            % number of particles to follow
            p.addParameter('ntag', 100, @(x) isfloat(x) && x > 0);
            % input list of tags
            p.addParameter('in_taglist', [], @(x) isfloat(x));
            
            
            
            % dump list for tracking
            p.addParameter('dump_list',0:100, @(x) isfloat(x));
            
            % dump in which tags will be searched (normally last dump)
            p.addParameter('reference_dump', 100, @(x) isfloat(x) && x >=0);
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            if isempty(fieldnames(p.Unmatched))
                unmatched = {};
            else
                unmatched = reshape([fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj = obj@OsirisDenormalizer(unmatched{:}); %Parse to superclass OsirisDenormalizer.m
            
            
            obj.trans_range             = p.Results.trans_range;
            obj.ntag                    = p.Results.ntag;
            obj.dump_list               = p.Results.dump_list;
            obj.in_taglist              = p.Results.in_taglist;
            obj.reference_dump          = p.Results.reference_dump;
            
        end % constructor
        
        function obj = select_tags(obj)
            
            obj.raw_dataset = 'tag';
            obj.getdata(); obj.assign_raw();
            obj.raw_dataset = 'x'; obj.direction = 'r';
            obj.getdata(); obj.assign_raw();
            obj.direction = 'z';
            obj.getdata(); obj.assign_raw();
            obj.r_raw = obj.denorm_distance(obj.nr_raw);
            obj.z_raw = obj.denorm_distance(obj.nz_raw);
            
            r_ind = (obj.r_raw > obj.trans_range(1)) & (obj.r_raw < obj.trans_range(2));
            
            % indices taken for xi = z - ct (usually it is the other way
            % around)
            if obj.xi_range(1) ~= 1e10
                obj.denorm_distance();
                z_ind = obj.z_raw > obj.dtime + obj.xi_range(1) & ... %small
                    obj.z_raw <= obj.dtime + obj.xi_range(2); % large
                
                z_selection = obj.z_raw(z_ind);
                r_selection = obj.r_raw(z_ind);
                unique_z = unique(z_selection);
                lunique_z = length(unique_z);
                
                rz_ind = r_ind & z_ind;
                
                if  lunique_z ~= 1
                    new_z_ind = obj.z_raw == max(unique_z);
                    % BEWARE: number below hardcoded for benchmark_etc
                    if sum(new_z_ind) >= 1494
                        rz_ind = r_ind & z_ind;
                    else
                        
                        z_in_unique = zeros(1,lunique_z);
                        
                        for zz = 1:lunique_z
                            z_in_unique(zz) = sum(obj.z_raw == unique_z(zz));
                        end
                        
                        % if all unique values in z have the same number of
                        % elements, we have an extra column
                        % (it can be that values in z change slightly)
                        if length(unique(z_in_unique)) == 1
                            % find the one farthest from front and eliminate its
                            % indices
                            for zz = 1:(lunique_z-1)
                                rz_ind = rz_ind & ~(obj.z_raw == min(unique_z));
                                unique_z(unique_z==min(unique_z)) = [];
                            end
                        else
                            pause;
                        end
                    end
                    
                    %                 if all(diff(unique_z) > abs(diff(obj.xi_range))/3)
                    %
                    %                 end
                end
                
                
                
                
                
            else
                rz_ind = r_ind;
            end
            tags_of_interest = obj.tag(1:2,rz_ind);

            if obj.ntag <= length(tags_of_interest)
                random_picker = randperm(length(tags_of_interest),obj.ntag);
                obj.selected_tags = tags_of_interest(1:2,random_picker); % pick 100 random particle from tag
            else
                obj.ntag = length(tags_of_interest);
                warning('More tags than particles, selecting all particles')
                obj.selected_tags = tags_of_interest;
            end
            
            
        end % select_tags
        
        function obj = find_tags(obj)
            
            obj.raw_dataset = 'tag';
            obj.getdata(); obj.assign_raw();
            obj.ind_tag = ismember(obj.tag',obj.in_taglist,'rows');
            %             ind_tag_temp = (obj.tag(1,:) == obj.in_taglist(ii,1)) & (obj.tag(2,:) == obj.in_taglist(ii,2)); % look for the selected tags in the tags of the current dump
            %             if ii == 1; obj.ind_tag = ind_tag_temp; end
            %             obj.ind_tag = or(obj.ind_tag,ind_tag_temp);
            %             fprintf('find %.d/%.d \n',ii,length(obj.in_taglist))
            
            %             a = 1;
        end % find_tags
        
        function obj = create_tracks(obj)
            
            %             r_pos = cell(length(obj.dump_list),1);
            obj.partrack = zeros(length(obj.dump_list),100);
            for n = 1:length(obj.dump_list)
                obj.dump = obj.dump_list(n);
                obj.raw_dataset = 'tag';
                obj.getdata(); obj.assign_raw();
                obj.raw_dataset = 'x'; obj.direction = 'r';
                obj.getdata(); obj.assign_raw();
                
                for rr = 1:length(rtags_oi)
                    ind_tag = (obj.tag(1,:) == rtags_oi(1,rr)) & (obj.tag(2,:) == rtags_oi(2,rr)); % look for the selected tags in the tags of the current dump
                    ind_check = sum(ind_tag);
                    obj.partrack(n,rr) = obj.nr_raw(ind_tag);
                    
                end % length tags
                
            end % for n dump_list
            
            obj.partrack = obj.denorm_distance(obj.partrack);
            
            
            
            
        end % create_tracks
        
        function obj = track_particles(obj)
            
        end
        
        function obj = plot_particles(obj)
            
            
        end
        
    end % ordinary methods
    
    
    methods(Static)
        
        function xxx
            
        end % xxx
        
    end % static methods
    
end % classdef