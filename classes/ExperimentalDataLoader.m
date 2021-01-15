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


classdef ExperimentalDataLoader < handle & OsirisDenormalizer
    
    
    properties(Constant, Hidden)
        px2cm = 0.0217/10;
        px2ps = 0.4121398108414873;
        
    end % constant properties
    
    properties
        
    end % properties
    
    properties (Hidden)
        
        % debugging
        sigma_z;
        
        % output
        SCI_nonaligned; % Streak Camera Image without information on the center or ionization front

        
    end % hidden properties
    
    methods
        
        % --CONSTRUCTOR--
        
        function obj = ExperimentalDataLoader(varargin)
            
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
            
            
            
            
        end % constructor
        
        
        
        function loadSCdata_nonaligned(obj)
            exp_data_format = 'mat'; % 'txt';
            switch exp_data_format
                case 'txt'
                    SCdata_filename = ['loading_files/',obj.datadir,'_stitch','.txt'];
                    fileID = fopen(SCdata_filename,'r');
                    formatSpec = '%f';
                    temp_SC = fscanf(fileID,formatSpec);
                    save(['loading_files/',obj.datadir,'_stitch.mat'],'temp_SC')
                case 'mat'
                    SCdata_filename = ['loading_files/',obj.datadir,'_stitch','.mat'];
                    load_temp = load(SCdata_filename);
                    temp_SC = load_temp.temp_SC;
            end
            verticalRes = 672;
            
            % laser positions above
            
            switch obj.datadir
                case 'gm20'
                    horizontalRes_gm20 = 1696;
                    obj.SCI_nonaligned = flipud(rot90(reshape(temp_SC,horizontalRes_gm20,verticalRes),3));
%                     obj.SCI = temp_SC2(:,44:1679);
%                     obj.SCcenter_px = 346;
                case 'gm10'
                    horizontalRes_gm10 = 1682;
                    obj.SCI_nonaligned = flipud(rot90(reshape(temp_SC,horizontalRes_gm10,verticalRes),3)); % 672,1682
%                     obj.SCI = temp_SC2(:,23:1658);
%                     obj.SCcenter_px = 342;
                case 'gm5'
                    horizontalRes_gm5 = 1657;
                    obj.SCI_nonaligned = flipud(rot90(reshape(temp_SC,horizontalRes_gm5,verticalRes),3));
%                     obj.SCI = temp_SC2(:,9:1644);
%                     obj.SCcenter_px = 336;
                case 'g0'
                    horizontalRes_gm0 = 1636;
                    obj.SCI_nonaligned = flipud(rot90(reshape(temp_SC,horizontalRes_gm0,verticalRes),3));
%                     obj.SCI = temp_SC2(:,1:1636);
%                     obj.SCcenter_px = 335;
                case 'gp5'
                    horizontalRes_gp5 = 1662;
                    obj.SCI_nonaligned = flipud(rot90(reshape(temp_SC,horizontalRes_gp5,verticalRes),3));
%                     obj.SCI = temp_SC2(:,12:1647);
%                     obj.SCcenter_px = 334;
                case 'gp10'
                    horizontalRes_gp10 = 1650;
                    obj.SCI_nonaligned = flipud(rot90(reshape(temp_SC,horizontalRes_gp10,verticalRes),3));
%                     obj.SCI = temp_SC2(:,13:1648);
%                     obj.SCcenter_px = 334;
                case 'gp15'
                    horizontalRes_gp15 = 1671;
                    obj.SCI_nonaligned = flipud(rot90(reshape(temp_SC,horizontalRes_gp15,verticalRes),3));
%                     obj.SCI = temp_SC2(:,21:1656);
%                     obj.SCcenter_px = 334;
                case 'gp20'
                    horizontalRes_gp20 = 1869;
                    obj.SCI_nonaligned = flipud(rot90(reshape(temp_SC,horizontalRes_gp20,verticalRes),3));
%                     obj.SCI = temp_SC2(:,200:1835);
%                     obj.SCcenter_px = 335;
            end % switch datadir
            


        end %loadSCdata_nonaligned
        
    end % ordinary methods
    
    
end % classdef