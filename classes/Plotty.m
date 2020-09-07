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
        
    end % constant properties
    
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
            p.addParameter('save_format', {'png','fig','eps'}, @(x) any(ismember(x,{'png','fig','vector',''})));
            % Specify if save or not
            p.addParameter('save_flag', true, @(x) islogical(x));
            % Specify if actually plot or not (the plot in questions varies from fucntion to function)
            p.addParameter('do_plot', true, @(x) islogical(x));
            % Specify if the units of the plot should be denormalized
            p.addParameter('denormalize_flag', true, @(x) islogical(x));
            % (plot_field_density) Specifiy which property to plot, or if both
            p.addParameter('property_plot', 'both', @(x) any(ismember(x,{'wakefields','density','both'})));
            % (plot_field_density) Specify if movie should be created
            p.addParameter('create_movie', false, @(x) islogical(x));
            % (plot_field_density) Specify if struct with movie frames
            % should be saved
            p.addParameter('save_movie_struct', false, @(x) islogical(x));
            
            % handle of the figure to save
            p.addParameter('fig_handle',[], @(x) ishghandle(x,'figure'));
            
            % flag to put a title in the figures
            p.addParameter('title_flag', true, @(x) islogical(x));
            
            % flag to make a pause in a series of plots
            p.addParameter('make_pause', true, @(x) islogical(x));
            
            
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
            
            
        end % constructor
        
        function obj = save_plot(obj)
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
                
            end % save flag
        end % save_plot
        
        function obj = waterfall_plot(obj,varargin)
            
            if isfloat(varargin{1})
                fig_number = varargin{1};
            else
                fig_number = randi;
            end
            
            obj.fig_handle = figure(fig_number);
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
            ylabel('propagation distance (m)');
            xlabel('\xi (cm)');
            
            
        end % waterfall_plot
        
        function obj = field_density_plot(obj)
            
            if obj.create_movie
                movie_dir = ['movies/field_density/'];
                struct_dir = ['save_files/field_density/'];
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
                f = 0;
                struct_movie(length(obj.dump_list)) = struct('cdata',[],'colormap',[]);
            end % create movie
            
            for n = 1:length(obj.dump_list)
                %initialize variables to make life simpler afterwards (less
                %switches and ifs)
                field_plot = 0; ax_field = 0;
                density_plot = 0; ax_density = 0;
                obj.dump = obj.dump_list(n);
                
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
                
                
                % mirroring for a better image, due to cylindrical symmetry
                field_half = field_plot;
                field_plot = [flip(field_plot);field_plot];
                density_plot = [flip(density_plot);density_plot];
                r_plot = [-max(obj.r),max(obj.r)]*10; % in mm
                z_plot = ([min(obj.z),max(obj.z)]-obj.simulation_window)/100; % in m
                
                
                % begin the plot
                
                
                fig_double = figure(1);
                if n == 1
                    fig_double.Units = 'normalized';
                    fig_double.OuterPosition = [0.1 0.3 0.8 0.5]; %[0.1 0.3 0.8 0.5]
                end
                
                if ismember(obj.property_plot,{'wakefields','both'})
                    
                    if obj.useAvg
                        skip_axis = 3;
                    else
                        skip_axis = 6;
                    end
                    
                    field_color_max = max(field_half(skip_axis:end,:),[],'all');
                    if field_color_max <= 0; field_color_max = 1; end
                    
                    % set the limits to the axes
                    ax_field = axes('parent',fig_double,'NextPlot','add');
                    ax_field.XLim = z_plot;
                    ax_field.YLim = r_plot;
                    
                    imagesc(ax_field,'XData',z_plot,'YData',r_plot,'CData',field_plot,[-field_color_max field_color_max]);
                    
                    c_field = colorbar('location','eastoutside');
                    colorbar_string = [obj.wakefields_direction,'. wakefields (MV/m)'];
                    c_field.Label.String = colorbar_string;
                    colormap(ax_field,bluewhitered);
                end
                
                if ismember(obj.property_plot,{'density','both'})
                    
                    % stablish the opaque index
                    % the standard deviation is used as a measure the avoid
                    % noisy peaks that sets a wrong scale for the opaqueness
                    % get 3 times the std deviation with no weights
                    meanstd_density = 3*std(density_plot,[],'all');
                    max_opaqueness = 1;
                    ind_opaque = max_opaqueness*density_plot;
                    ind_opaque(density_plot > meanstd_density) = max_opaqueness*meanstd_density;
                    ind_opaque = ind_opaque/max(ind_opaque,[],'all');
                    
                    ax_density = axes('parent',fig_double,'NextPlot','add','color','none');
                    
                    % set limits to the axis
                    ax_density.XLim = z_plot;
                    ax_density.YLim = r_plot;
                    
                    imagesc(ax_density,'XData',z_plot,'YData',r_plot,'CData',double(density_plot>0),'alphadata',ind_opaque);
                    grad = colorGradient([1 1 1],[0 0 0],2);
                    colormap(ax_density,grad);
                end
                
                if strcmp(obj.property_plot,'both')
                    ax_field.Position = ax_density.Position;
                    linkaxes([ax_field,ax_density],'xy')
                end
                
                if obj.title_flag
                    title(['propagation dist. = ',num2str(obj.propagation_distance/100,2),' m',''])
                end
                ylabel('transverse distance (mm)')
                xlabel('propagation distance (m)');
                
                line_flag = false;
                if line_flag
                   hold on
%                    plot(max(z_plot)-[0.003 0.003],r_plot,'LineWidth',2,'k')
                   xline(max(z_plot)-0.0025,'k','LineWidth',2)
                   hold off
                end
                
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
                    struct_movie(f+1) = getframe(gcf);
                    f = f + 1;
                end % create movie
                
                if n < length(obj.dump_list)
                    clf
                end % clear figure
                obj.progress_dump('Plotting 2D',n,length(obj.dump_list))
            end % length dump list
            
            if obj.create_movie
                writeVideo(video,struct_movie);
                save(['save_files/field_density/struct',obj.datadir,...
                    '_',obj.wakefields_direction,'.mat'],'struct_movie')
                close(video);
            end % if create movie
            
        end % field_density_plot
        
        
        function obj = lineout_plot(obj)
            
            if obj.create_movie
                movie_dir = ['movies/lineout/'];
                struct_dir = ['save_files/lineout/'];
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
                f = 0;
                struct_movie(length(obj.dump_list)) = struct('cdata',[],'colormap',[]);
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
                            lineout_plot = obj.ndataOut;
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
                            lineout_plot = obj.ndataOut;
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
                            lineout_plot = obj.([obj.wakefields_direction,'field']);
                            density_plot = obj.(obj.species);
                        end % denorm flag
                        
                end % switch property_plot

                % begin the plot
                
                
                fig_lineout = figure(1);
                if n == 1
                    fig_lineout.Units = 'normalized';
                    fig_lineout.OuterPosition = [0.1 0.3 0.8 0.5]; %[0.1 0.3 0.8 0.5]
                end
                
                if ismember(obj.property_plot,{'wakefields','both'})
                    
                    plot((obj.z-obj.simulation_window)/100,lineout_plot,'LineWidth',1);
                    xlim(([min(obj.z),max(obj.z)]-obj.simulation_window)/100);
                    
                    % set the limits to the axes
%                     ax_field = axes('parent',fig_double,'NextPlot','add');
%                     ax_field.XLim = z_plot;
                                        
                end
                
                if ismember(obj.property_plot,{'density','both'})
                    
                    % stablish the opaque index
                    % the standard deviation is used as a measure the avoid
                    % noisy peaks that sets a wrong scale for the opaqueness
                    % get 3 times the std deviation with no weights
                    meanstd_density = 0.5*std(density_plot,[],'all');
                    max_opaqueness = 1;
                    ind_opaque = max_opaqueness*density_plot;
                    ind_opaque(density_plot > meanstd_density) = max_opaqueness*meanstd_density;
                    ind_opaque = ind_opaque/max(ind_opaque,[],'all');
                    
                    ax_density = axes('parent',fig_lineout,'NextPlot','add','color','none');
                    
                    % set limits to the axis
                    ax_density.XLim = z_plot;
                    ax_density.YLim = r_plot;
                    
                    imagesc(ax_density,'XData',z_plot,'YData',r_plot,'CData',double(density_plot>0),'alphadata',ind_opaque);
                    grad = colorGradient([1 1 1],[0 0 0],2);
                    colormap(ax_density,grad);
                end
                
                if strcmp(obj.property_plot,'both')
                    ax_field.Position = ax_density.Position;
                    linkaxes([ax_field,ax_density],'xy')
                end
                
                if obj.title_flag
                    title(['propagation dist. = ',num2str(obj.propagation_distance/100,2),' m',' (on axis)']) %(r = 1/kp)
                end
                ylabel([obj.wakefields_direction,'. fields (MV/m)'])
                xlabel('propagation distance (m)');
                
                line_flag = false;
                if line_flag
                   hold on
%                    plot(max(z_plot)-[0.003 0.003],r_plot,'LineWidth',2,'k')
                   xline(max(z_plot)-0.0025,'k','LineWidth',2)
                   hold off
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
                    struct_movie(f+1) = getframe(gcf);
                    f = f + 1;
                end % create movie
                
                if n < length(obj.dump_list)
                    clf
                end % clear figure
                obj.progress_dump('Plotting 2D',n,length(obj.dump_list))
            end % length dump list
            
            if obj.create_movie
                writeVideo(video,struct_movie);
                save(['save_files/field_density/struct',obj.datadir,...
                    '_',obj.wakefields_direction,'.mat'],'struct_movie')
                close(video);
            end % if create movie
            
        end % field_density_plot

        
    end % ordinary methods
    
    
    methods(Static)
        
        
    end % static methods
    
end % classdef