%__________________________________________________________________________
% Class with the main plotting function for the AWAKE Osiris Analysis
% Matlab Package
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/07/2020
%__________________________________________________________________________

% Input:
% - plots directory, plot name, saving instructions
%
% Output:
% - saved plot (no variable)
%
% Methods
% - save_plot: created saving directory and saves plot according to
% instrucions
%
%


classdef Plotty < handle & OsirisDenormalizer
    
    properties(Constant, Hidden)
        
        % Physical constants
        
        
        % values from: https://physics.nist.gov/cuu/Constants/index.html
        
    end % constant properties
    
    
    properties(Hidden)
        
        % waterfall input
        dump_list;
        do_plot;
        
        % waterfall output
        waterfall_xi;
        waterfall_z;
        
        % field density plot input
        save_movie_struct;
        
        % flags
        title_flag;
        make_pause;
        
        % plot scale
        plot_scale;
        
        % figure number
        fig_number;
        
        % field density plot
        trans_lines_position;
        mirror_flag;
        
        ax_field;
        ax_density;
        tile_handle;
        include_lineout;
        trans_lines_2D_flag
        lines_2D_flag;
        
        % movie
        frame_counter = 0;
        
    end % hidden properties
    
    properties
        
        % input, description in parser
        plots_dir;
        plot_name;
        save_flag;
        save_format;
        
        % plot field density input
        property_plot;
        denormalize_flag;
        create_movie;
        
        % output
        waterfall_mat;
        waterfall_handle;
        fig_handle;
        plot_handle;
        
        struct_movie;
        
    end % properties
    
    methods
        
        % --CONSTRUCTOR--
        
        function obj = Plotty(varargin)
            
            
            % parse input to load files
            p = inputParser;
            
            % directory where the plot should be saved
            p.addParameter('plots_dir', 'plots', @(x) ischar(x));
            p.addParameter('plot_name', 'plot', @(x) ischar(x));
            
            
            % for the case of plots that need a dump list (ex.: waterfall)
            p.addParameter('dump_list', 1:1:100, @(x) isfloat(x) && all(x >= 0));
            
            % Specify which file formats to save (setting to empty ('')
            % will not save anything)
            p.addParameter('save_format', {'png','fig','eps'}, @(x) any(ismember(x,{'png','fig','eps','pdf',''})));
            % Specify if save or not
            p.addParameter('save_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            % Specify if actually plot or not (the plot in question varies from function to function)
            p.addParameter('do_plot', true, @(x) islogical(x) || x == 0 || x == 1);
            % Specify if the units of the plot should be denormalized
            p.addParameter('denormalize_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            % (plot_field_density) Specifiy which property to plot, or if both
            p.addParameter('property_plot', 'both', @(x) any(ismember(x,{'wakefields','density','both'})));
            % (plot_field_density) Specify if movie should be created
            p.addParameter('create_movie', false, @(x) islogical(x) || x == 0 || x == 1);
            % (plot_field_density) Specify if struct with movie frames
            % should be saved
            p.addParameter('save_movie_struct', false, @(x) islogical(x) || x == 0 || x == 1);
            
            % handle of the figure to save
            p.addParameter('fig_handle',[], @(x) ishghandle(x,'figure'));
            
            % flag to put a title in the figures
            p.addParameter('title_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            
            % flag to make a pause in a series of plots
            p.addParameter('make_pause', true, @(x) islogical(x) || x == 0 || x == 1);
            
            % plot scale
            p.addParameter('plot_scale', 'linear', @(x) any(ismember(x,{'linear','log'})));
            
            % figure number
            p.addParameter('fig_number', 0, @(x) isfloat(x) && x > 0);
            
            % field density plot options
            p.addParameter('mirror_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            if isempty(fieldnames(p.Unmatched))
                unmatched = {};
            else
                unmatched = reshape([fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj = obj@OsirisDenormalizer(unmatched{:}); %Parse to superclass OsirisDenormalizer.m
            
            
            obj.plots_dir             = p.Results.plots_dir;
            obj.plot_name             = p.Results.plot_name;
            obj.dump_list             = p.Results.dump_list;
            obj.save_format           = p.Results.save_format;
            obj.save_flag             = p.Results.save_flag;
            obj.do_plot               = p.Results.do_plot;
            obj.fig_handle            = p.Results.fig_handle;
            obj.denormalize_flag      = p.Results.denormalize_flag;
            obj.property_plot         = p.Results.property_plot;
            obj.create_movie          = p.Results.create_movie;
            obj.save_movie_struct     = p.Results.save_movie_struct;
            obj.title_flag            = p.Results.title_flag;
            obj.make_pause            = p.Results.make_pause;
            obj.plot_scale            = p.Results.plot_scale;
            obj.fig_number            = p.Results.fig_number;
            obj.mirror_flag           = p.Results.mirror_flag;
        end % constructor
        
        function obj = save_plot(obj,varargin)
            
            if nargin > 1
                obj.fig_handle = varargin{1};
            end
            
            plotsdir = 'plots/';
            if ~isfolder([plotsdir,obj.plots_dir]); mkdir([plotsdir,obj.plots_dir]); end
            
            if obj.save_flag
                % save as png image
                if any(ismember(obj.save_format,{'png'}))
                    if ~isfolder([plotsdir,obj.plots_dir,'/png/']); mkdir([plotsdir,obj.plots_dir,'/png/']); end
                    exportgraphics(obj.fig_handle,[plotsdir,obj.plots_dir,'/png/',obj.plot_name,'.png'])
                end
                
                % save as matlab figure
                if any(ismember(obj.save_format,{'fig'}))
                    if ~isfolder([plotsdir,obj.plots_dir,'/fig/']); mkdir([plotsdir,obj.plots_dir,'/fig/']); end
                    savefig(obj.fig_handle,[plotsdir,obj.plots_dir,'/fig/',obj.plot_name,'.fig'])
                end
                
                % save as vector image
                if any(ismember(obj.save_format,{'eps'}))
                    if ~isfolder([plotsdir,obj.plots_dir,'/eps/']); mkdir([plotsdir,obj.plots_dir,'/eps/']); end
                    exportgraphics(obj.fig_handle,[plotsdir,obj.plots_dir,'/eps/',obj.plot_name,'.eps'],'ContentType','vector')
                end
                
                % save as vector image
                if any(ismember(obj.save_format,{'pdf'}))
                    if ~isfolder([plotsdir,obj.plots_dir,'/pdf/']); mkdir([plotsdir,obj.plots_dir,'/pdf/']); end
                    exportgraphics(obj.fig_handle,[plotsdir,obj.plots_dir,'/pdf/',obj.plot_name,'.pdf'],'ContentType','vector')
                end
                
            end % save flag
        end % save_plot
        
        function obj = waterfall_plot(obj)
            
            if obj.fig_number > 0
                fig_waterfall = figure(obj.fig_number);
            else
                fig_waterfall = figure;
            end
            
            obj.fig_handle = fig_waterfall;
            obj.waterfall_handle = imagesc(obj.waterfall_xi,...
                obj.waterfall_z,rot90(obj.waterfall_mat,2));
            switch obj.property
                case 'fields'
                    colormap(bluewhitered);
                case 'density'
                    c = gray;
                    c = flipud(c);
                    colormap(c);
            end
            drawnow;
            ax = gca;
            ax.XDir = 'reverse';
            ax.YDir = 'normal';
            cbar = colorbar;
            switch obj.property
                case 'fields'
                    switch obj.wakefields_direction
                        case 'long'
                            cbar.Label.String = 'E_z (MV/m)';
                        case 'trans'
                            cbar.Label.String = 'W_r (MV/m)';
                        otherwise
                            cbar.Label.String = 'E_z (MV/m)';
                    end
                    
                case 'density'
                    cbar.Label.String = 'charge (e)';
            end
            ylim([obj.waterfall_z(1) obj.waterfall_z(2)])
            ylabel('z (m)');
            xlabel('\xi (cm)');
            
            
        end % waterfall_plot
        
        function obj = plot_field(obj,varargin)
            
            % parse input to load files
            p = inputParser;
            
            % directory where the plot should be saved
            p.addParameter('field_plot', 0, @(x) isfloat(x));
            p.addParameter('r_plot', 0, @(x) isfloat(x));
            p.addParameter('z_plot', 0, @(x) isfloat(x));
            
            p.parse(varargin{:});
            
            field_plot         = p.Results.field_plot;
            r_plot             = p.Results.r_plot;
            z_plot             = p.Results.z_plot;
            
            if length(field_plot) <= 1
                [field_plot,~,r_plot,z_plot] = obj.load_data_field_density_plot();
            elseif all(r_plot == 0) || all(z_plot == 0)
                error('Please also give r_plot and z_plot.')
            end
            
            % mirroring for a better image, due to cylindrical symmetry
            
            field_plot_mirrored = [flip(field_plot);field_plot];
            
            % Establish the maximum field value
            % the standard deviation is used as a measure the avoid
            % noisy peaks that sets a wrong scale for the opaqueness
            % get 3 times the std deviation with no weights
            meanstd_field = 5*std(abs(field_plot),[],'all')+eps;
            
            
            % set the limits to the axes
            obj.ax_field = axes('parent',obj.tile_handle,'NextPlot','add');
            
            
            obj.ax_field.XLim = [min(z_plot),max(z_plot)];
            obj.ax_field.YLim = r_plot;
            obj.ax_field.XDir = 'reverse';
            imagesc(obj.ax_field,'XData',z_plot,'YData',r_plot,'CData',field_plot_mirrored,[-meanstd_field meanstd_field]);
            c_field = colorbar('location','eastoutside');
            colorbar_string = [obj.wakefields_direction,'. wakefields (MV/m)'];
            c_field.Label.String = colorbar_string;
            colormap(obj.ax_field,bluewhitered);
        end % plot field
        
        function obj = plot_density(obj,varargin)
            
            % parse input to load files
            p = inputParser;
            
            % directory where the plot should be saved
            p.addParameter('density_plot', 0, @(x) isfloat(x));
            p.addParameter('r_plot', 0, @(x) isfloat(x));
            p.addParameter('z_plot', 0, @(x) isfloat(x));
            
            p.parse(varargin{:});
            
            density_plot       = p.Results.density_plot;
            r_plot             = p.Results.r_plot;
            z_plot             = p.Results.z_plot;
            
            switch obj.species
                case {'plasma_positrons','plasma_electrons'}
                    density_plot = density_plot - median(density_plot,'all');
            end
            
            switch obj.plot_scale
                case 'log'
                    density_plot_l = log(density_plot+1);
                case 'linear'
                    density_plot_l = density_plot;
            end % switch plot scale
            %  clear density_plot
            
            if obj.mirror_flag
                density_plot_mirrored = [flip(density_plot_l);density_plot_l];
            else
                density_plot_mirrored = density_plot_l;
            end
            clear density_plot_l
            
            % stablish the opaque index
            % the standard deviation is used as a measure the avoid
            % noisy peaks that sets a wrong scale for the opaqueness
            % get 3 times the std deviation with no weights
            meanstd_density = 3*std(density_plot_mirrored,[],'all')+eps;
            max_opaqueness = 1;
            ind_opaque = max_opaqueness*density_plot_mirrored;
            ind_opaque(density_plot_mirrored > meanstd_density) = max_opaqueness*meanstd_density;
            ind_opaque = ind_opaque/max(ind_opaque,[],'all');
            if isempty(obj.tile_handle)
                obj.tile_handle = tiledlayout(1,1);
            end
            obj.ax_density = axes('parent',obj.tile_handle,'NextPlot','add','color','none');
            % set limits to the axis
            obj.ax_density.XLim = [min(z_plot),max(z_plot)];
            obj.ax_density.YLim = [min(r_plot),max(r_plot)];
            obj.ax_density.XDir = 'reverse';
            
            
            switch obj.plot_scale
                case 'log'
                    imagesc(obj.ax_density,'XData',z_plot,'YData',r_plot,'CData',density_plot_mirrored);
                case 'linear'
                    imagesc(obj.ax_density,'XData',z_plot,'YData',r_plot,'CData',double(density_plot_mirrored>0),'alphadata',ind_opaque);
%                     grad = colorGradient([1 1 1],[0 0 0],2);
                    grad = [1 1 1; 0 0 0];
                    colormap(obj.ax_density,grad);
            end % switch plot scale
        end % plot density
        
        function obj = plot_field_density(obj,varargin)
            
            % parse input to load files
            p = inputParser;
            
            % directory where the plot should be saved
            p.addParameter('field_plot', 0, @(x) isfloat(x));
            p.addParameter('density_plot', 0, @(x) isfloat(x));
            p.addParameter('r_plot', 0, @(x) isfloat(x));
            p.addParameter('z_plot', 0, @(x) isfloat(x));
            p.addParameter('mirror_flag', true, @(x) islogical(x) || x == 0 || x == 1);
            p.addParameter('trans_lines_2D_flag', false, @(x) islogical(x) || x == 0 || x == 1);
            p.addParameter('trans_lines_position', [0,0], @(x) isfloat(x));
            p.addParameter('ymax_density_profile', 8e14, @(x) isfloat(x));
            p.addParameter('include_lineout', 'no', @(x) any(ismember(x,{'no','both','field_lineout','density_profile'})));
                        
            p.parse(varargin{:});
            
            field_plot         = p.Results.field_plot;
            density_plot       = p.Results.density_plot;
            r_plot             = p.Results.r_plot;
            z_plot             = p.Results.z_plot;
            
            ymax_density_profile ...
                = p.Results.ymax_density_profile;
            obj.mirror_flag    = p.Results.mirror_flag;
            obj.trans_lines_2D_flag ...
                               = p.Results.trans_lines_2D_flag;
            obj.trans_lines_position ...
                = p.Results.trans_lines_position;
            obj.include_lineout  = p.Results.include_lineout;

            
            if obj.create_movie
                [~,video] = obj.setup_movie();
            end % create movie
            
            for n = 1:length(obj.dump_list)
                
                obj.dump = obj.dump_list(n);
                
                if ((numel(field_plot) == 1) && (numel(density_plot) == 1)) && (field_plot == 0 && density_plot == 0)
                    [field_plot,density_plot,r_plot,z_plot] = obj.load_data_field_density_plot();
                    obj.mirror_flag  = true;
                end
                
                % begin the plot
                
                if obj.fig_number > 0
                    fig_double = figure(obj.fig_number);
                else
                    fig_double = figure;
                end
                
                switch obj.include_lineout
                    case 'density_profile'
                        obj.tile_handle = tiledlayout(fig_double,2,1);
                        plot_position = [0.0471    0.2560    0.9081    0.5111];
                    case 'field_lineout'
                        obj.tile_handle = tiledlayout(fig_double,2,1);
                        plot_position = [0.0471    0.2560    0.9081    0.5111];
                    case 'both'
                        obj.tile_handle = tiledlayout(fig_double,3,1);
                        plot_position = [0.0346    0.1231    0.9008    0.7356];
                    otherwise
                        obj.tile_handle = tiledlayout(fig_double,1,1);
                        plot_position = [0.1 0.1 0.8 0.5];
                end
                obj.tile_handle.TileSpacing = 'compact';
                obj.tile_handle.Padding = 'compact';
                
                if n == 1
                    fig_double.Units = 'normalized';
                    fig_double.OuterPosition = plot_position;
                end
                
                if ismember(obj.property_plot,{'wakefields','both'})
                    obj.plot_field('field_plot',field_plot,'r_plot',r_plot,'z_plot',z_plot);
                    obj.ax_field.FontSize = 12;                
                end
                
                if ismember(obj.property_plot,{'density','both'})
                    obj.plot_density('density_plot',density_plot,'r_plot',r_plot,'z_plot',z_plot);
                    obj.ax_density.FontSize = 12;
                end
                
                if strcmp(obj.property_plot,'both')
                    obj.ax_field.Position = obj.ax_density.Position;
                    linkaxes([obj.ax_field,obj.ax_density],'xy')
                end
                
                if obj.title_flag
                    title(['propagation dist. = ',num2str(obj.propagation_distance/100,3),' m',''])
                end
                
                if obj.trans_lines_2D_flag
                    hold on
                    yline(obj.trans_lines_position(1)*10,'r','LineWidth',1)
                    yline(obj.trans_lines_position(2)*10,'r','LineWidth',1)
                    hold off
                end
                
                ylabel('r (mm)')
                
                if ismember(obj.include_lineout,{'density_profile','both'})
                    
                    obj.ax_density.XTickLabel = [];
                    obj.ax_field.XTickLabel = [];
                                        
                    r_lineplot = linspace(0,max(r_plot),size(density_plot,1));
                    ir = (r_lineplot < obj.trans_lines_position(2)*10) & (r_lineplot > obj.trans_lines_position(1)*10);
                    long_profile = obj.cylindrical_radial_integration(r_lineplot(ir),density_plot(ir,:),'just_sum'); %just sum
                    
                    ax_longprofile = nexttile;
                    z_lineplot = linspace(max(z_plot),min(z_plot),length(long_profile));
                    plongprofile = plot(z_lineplot,long_profile);
                    ax_longprofile.XDir = 'reverse';
                    
                    if strcmp(obj.include_lineout,'density_profile')
%                         ax_longprofile.XTick = [z_lineplot];
%                         ax_longprofile.XTickLabel = [z_lineplot];
                        xlabel('ct - z (cm)');
                    end
                    
                    xlim(obj.ax_density.XLim)
                    ylim([0 ymax_density_profile])
                    
                    ylabel({'charge','density (a. u.)'});
                    ax_longprofile.FontSize = 12;
                    
                end % if include lineout
                
                if ismember(obj.include_lineout,{'field_lineout','both'})
                    
                    ax_lineout = nexttile;
                    % HARDCODED 1111
                    lineout = field_plot(10,:);
                    plineout = plot(obj.dtime + obj.simulation_window - obj.z,lineout);
                    ax_lineout.XDir = 'reverse';
                    xlim(obj.ax_density.XLim);
                                        
                    if strcmp(obj.include_lineout,'both')
%                         ax_longprofile.XTickLabel = [];
%                         ax_longprofile.XTick = [z_lineplot];
                        
                    end
                    
                    ylabel({'E_z (MV/m)'});
                    xlabel('ct - z (cm)');
                    %                     ylim([-200 200])
                end % if include long profile
                
                drawnow;
                
                if obj.make_pause && (n < length(obj.dump_list))
                    pause;
                end
                
                obj.fig_handle = fig_double;
                
                obj.plot_name = [obj.datadir,obj.property_plot,'n',num2str(obj.dump),...
                    'xi',num2str(round(obj.xi_range(1))),'xi',...
                    num2str(round(obj.xi_range(2))),'t',num2str(round(obj.trans_range(1))),...
                    't',num2str(round(obj.trans_range(2)))];
                
                obj.save_plot();
                
                
                if obj.create_movie
                    obj.struct_movie(obj.frame_counter+1) = getframe(gcf);
                    obj.frame_counter = obj.frame_counter + 1;
                end % create movie
                
                if n < length(obj.dump_list)
                    clf
                    density_plot = 0;
                    field_plot = 0;
                    r_plot = 0;
                    z_plot = 0;
                end % clear figure
                obj.progress_dump('Plotting 2D',n,length(obj.dump_list))
            end % length dump list
            
            if obj.create_movie
                writeVideo(video,obj.struct_movie);
                struct_movie_save = obj.struct_movie;
                save(['save_files/field_density/struct',obj.datadir,...
                    '_',obj.wakefields_direction,'.mat'],'struct_movie_save')
                close(video);
            end % if create movie
            
        end % plot_field_density
        
        function [obj,video] = setup_movie(obj)

            movie_dir = ['movies/field_density',obj.include_lineout,'/'];
            struct_dir = ['save_files/field_density',obj.include_lineout,'/'];
            if ~isfolder(movie_dir)
                mkdir(movie_dir);
            end
            if ~isfolder(struct_dir)
                mkdir(struct_dir);
            end
            video = VideoWriter([movie_dir,...
                obj.wakefields_direction,obj.datadir,obj.property_plot,...
                'xi',num2str(round(obj.xi_range(1))),'xi',...
                num2str(round(obj.xi_range(2))),'t',num2str(round(obj.trans_range(1))),...
                't',num2str(round(obj.trans_range(2))),'.avi']);
            video.FrameRate = 4;
            open(video);
            struct_movie_temp(length(obj.dump_list)) = struct('cdata',[],'colormap',[]);
            obj.struct_movie = struct_movie_temp;
        end
        
        function [field_plot,density_plot,r_plot,z_plot] = load_data_field_density_plot(obj)
            %initialize variables to make life simpler afterwards (less
            %switches and ifs)
            field_plot = 0; 
            density_plot = 0;
            
            switch obj.property_plot
                case 'wakefields'
                    obj.property = 'fields';
                    obj.getdata();
                    if obj.denormalize_flag
                        obj.trim_data();
                        obj.denorm_Efield();
                        field_plot = obj.([obj.wakefields_direction,'field']);
                    else
                        field_plot = obj.ndataOut;
                    end
                    
                case 'density'
                    obj.property = 'density';
                    obj.getdata();
                    if obj.denormalize_flag
                        obj.trim_data();
                        obj.denorm_density();
                        density_plot = obj.(obj.species);
                    else
                        density_plot = obj.ndataOut;
                    end % denorm flag
                case 'both'
                    obj.property = 'fields';
                    obj.getdata();
                    if obj.denormalize_flag
                        obj.trim_data();
                    else
                        field_plot = obj.ndataOut;
                    end
                    
                    obj.property = 'density';
                    obj.getdata();
                    if obj.denormalize_flag
                        obj.trim_data();
                    else
                        density_plot = obj.ndataOut;
                    end
                    
                    if obj.denormalize_flag
                        obj.denorm_Efield();
                        obj.denorm_density();
                        field_plot = obj.([obj.wakefields_direction,'field']);
                        density_plot = obj.(obj.species);
                    end % denorm flag
                    
            end % switch property_plot
            
            %             intnear1500 = floor(length(density_plot)/1500);
            %             density_plot = density_plot(:,end-intnear1500*1500+1:end);
            %             density_plot = squeeze(sum(reshape(density_plot,size(density_plot,1),intnear1500,[]),2));
            
            r_plot = [-max(obj.r),max(obj.r)]*10; % in mm
            % z_plot = ([min(obj.z),max(obj.z)]-obj.simulation_window)/100; % in m
            z_plot_temp = obj.dtime+obj.simulation_window-obj.z;
            z_plot = [z_plot_temp(1),z_plot_temp(end)];
        end % load data field density plot
        % ------------------------------------------------------------------------------
        function obj = plot_lineout(obj) 
            
            if obj.create_movie
                [~,video] = obj.setup_movie();
            end % create movie
            
            for n = 1:length(obj.dump_list)
                %initialize variables to make life simpler afterwards (less
                %switches and ifs)
                lineout_plot = 0; ax_field = 0;
                density_plot = 0; ax_density = 0;
                obj.dump = obj.dump_list(n);
                
                switch obj.property_plot
                    case 'wakefields'
                        obj.property = 'fields';
                        obj.getlineout();
                        if obj.denormalize_flag
                            obj.denorm_Efield();
                            lineout_plot = obj.lineout;
                        else
                            lineout_plot = obj.nlineout;
                        end
                    case 'density'
                        obj.property = 'density';
                        obj.getdata();
                        if obj.denormalize_flag
                            obj.trim_data();
                            obj.denorm_density();
                            density_plot = obj.(obj.species);
                        else
                            density_plot = obj.ndataOut;
                        end % denorm flag
                        lineout_plot = obj.cylindrical_radial_integration(obj.r,density_plot,'just_sum');
                        
                end % switch property_plot
                
                % select the axis
                
                if obj.denormalize_flag
                    z_plot = (obj.dtime + obj.simulation_window - obj.z); % cm
                else
                    z_plot = (obj.ntime + obj.n_simulation_window - obj.nz)/(2*pi);
                end
                
                % begin the plot
                
                
                fig_lineout = figure(1);
                if n == 1
                    fig_lineout.Units = 'normalized';
                    fig_lineout.OuterPosition = [0.1 0.3 0.8 0.5]; %[0.1 0.3 0.8 0.5]
                end
                
                obj.plot_handle = plot(z_plot,lineout_plot,'k','LineWidth',1.2);
                xlim(([min(z_plot),max(z_plot)]));
                ax_lo = gca;
                ax_lo.XDir = 'reverse';
              
                if obj.denormalize_flag
                    if obj.title_flag
                        title(['propagation dist. = ',num2str(obj.propagation_distance/100,2),' m',''])
                    end
                    ylabel([obj.wakefields_direction,'. fields (MV/m)'])
                    xlabel('ct - z (cm)');
                else
                    if obj.title_flag
                        title(['propagation dist. = ',num2str(obj.n_propagation_distance,2),'']) 
                    end
                    ylabel([obj.wakefields_direction,'. fields'])
                    xlabel('ct - z (\lambda_p)');

                end
                
                drawnow;
                
                if obj.make_pause && (n < length(obj.dump_list))
                    pause;
                end
                
                obj.fig_handle = fig_lineout;
                
                obj.plot_name = [obj.datadir,obj.property_plot,'n',num2str(obj.dump),...
                    'xi',num2str(round(obj.xi_range(1))),'xi',...
                    num2str(round(obj.xi_range(2))),'t',num2str(round(obj.trans_range(1))),...
                    't',num2str(round(obj.trans_range(2)))];
                obj.save_plot();
                
               if obj.create_movie
                    obj.struct_movie(obj.frame_counter+1) = getframe(gcf);
                    obj.frame_counter = obj.frame_counter + 1;
                end % create movie
                
                if n < length(obj.dump_list)
                    clf
                end % clear figure
                obj.progress_dump('Plotting lineout',n,length(obj.dump_list))
            end % length dump list
            
        end % lineout plot
        
    end % ordinary methods
    
    
    methods(Static)
        
        
    end % static methods
    
end % classdef