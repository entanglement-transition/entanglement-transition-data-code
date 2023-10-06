classdef app_surface_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ui_figure              matlab.ui.Figure
        render_time            matlab.ui.control.Label
        scene_time             matlab.ui.control.Label
        panel_7                matlab.ui.container.Panel
        surf_colormap          matlab.ui.control.DropDown
        ColormapDropDownLabel  matlab.ui.control.Label
        surf_size              matlab.ui.control.Spinner
        SizeSpinnerLabel       matlab.ui.control.Label
        surf_texture_even      matlab.ui.control.DropDown
        EvenDropDownLabel      matlab.ui.control.Label
        surf_texture_odd       matlab.ui.control.DropDown
        OddDropDownLabel       matlab.ui.control.Label
        surf_smooth            matlab.ui.control.CheckBox
        panel_6                matlab.ui.container.Panel
        bt_render              matlab.ui.control.Button
        cam_responsive         matlab.ui.control.CheckBox
        panel_5                matlab.ui.container.Panel
        light_4                matlab.ui.control.CheckBox
        light_3                matlab.ui.control.CheckBox
        light_2                matlab.ui.control.CheckBox
        light_1                matlab.ui.control.CheckBox
        panel_4                matlab.ui.container.Panel
        plane_xy_enable        matlab.ui.control.CheckBox
        plane_yz_enable        matlab.ui.control.CheckBox
        plane_xz_enable        matlab.ui.control.CheckBox
        grid_xy_enable         matlab.ui.control.CheckBox
        grid_yz_enable         matlab.ui.control.CheckBox
        grid_xz_enable         matlab.ui.control.CheckBox
        axis_enable            matlab.ui.control.CheckBox
        image_panel            matlab.ui.container.Panel
        image                  matlab.ui.control.Image
        panel_3                matlab.ui.container.Panel
        ResetButton            matlab.ui.control.Button
        cam_type               matlab.ui.control.DropDown
        cam_angle              matlab.ui.control.Spinner
        AngleSpinnerLabel      matlab.ui.control.Label
        panel_2                matlab.ui.container.Panel
        cam_look_at_z          matlab.ui.control.Spinner
        xLabel_3               matlab.ui.control.Label
        cam_look_at_y          matlab.ui.control.Spinner
        xLabel_2               matlab.ui.control.Label
        cam_look_at_x          matlab.ui.control.Spinner
        xLabel                 matlab.ui.control.Label
        panel                  matlab.ui.container.Panel
        cam_loc_z              matlab.ui.control.Spinner
        xLabel_6               matlab.ui.control.Label
        cam_loc_y              matlab.ui.control.Spinner
        xLabel_5               matlab.ui.control.Label
        cam_loc_x              matlab.ui.control.Spinner
        xLabel_4               matlab.ui.control.Label
    end

    properties (Access = public)
        povlab_dir = "";
        pl povlab = 0;

        cam_location_default = [-7 -7 3];
        cam_location = 0;

        cam_look_at_default = [0 0 0];
        cam_look_at = 0;

        need_create_scene = false;
    end
    
    methods (Access = private)
        % Init pl renderer interface
        function init_pov(app)
            app.povlab_dir = get_root_dir();
            povray_out_dir = app.povlab_dir + "/apps/out";
            povray_version = "3.7";
            
            if isunix
                povray_path = '"/usr/local/bin/povray +A -W1920 -H1080 -L'+ povray_out_dir + '"';
            elseif ispc
                % Windows version of Povray is not supporting image rendering parametrs in command line, so please
                % set up Povray from menu 'Render/Edit Settings' before using Povlab
                povray_path = "C:/Program Files/POV-Ray/v3.7/bin/pvengine64.exe";
            else
                disp('Platform not supported');
            end
            
            app.pl = povlab( povray_version,...
                             povray_path, ...
                             povray_out_dir);

            % Copy 'povlab.inc' to output directory
            povlab_include_file = get_root_dir() + "/include/povlab.inc";
            copyfile(povlab_include_file, povray_out_dir);
            
        end

        % Create camera
        function create_camera(app)
            app.pl.include_begin('camera');
                app.pl.camera('angle', app.cam_angle.Value, 'location', app.cam_location, 'look_at', app.cam_look_at, 'type', app.cam_type.Value);
            app.pl.include_end();
        end
        
        % Create light
        function create_lights(app)
            app.pl.include_begin('lights');
                if(app.light_1.Value)
                    app.pl.light('location', [-10 -17 7], 'color', [1 1 1], 'shadowless', true);
                end
                if(app.light_2.Value)
                    app.pl.light('location', [-10 10 30],  'color', [0.8 0.8 0.8], 'shadowless', true);
                end
                if(app.light_3.Value)
                    app.pl.light('location', [100 200 300], 'color', [0.4 0.4 0.4], 'shadowless', true);
                end
                if(app.light_4.Value)
                    app.pl.light('location', [100 -200 -300], 'color', [1 1 1], 'shadowless', true);
                end
            app.pl.include_end();
        end

        % Create helpers
        function create_helpers(app)
            app.pl.include_begin('helpers');
                % Axis
                % TODO: Move textures to materials.inc
                % tex_axis_gray = app.pl.declare("tex_axis_gray", app.pl.texture('pigment', [0.5 0.5 0.5], 'finish', "phong 1 reflection {0.10 metallic 0.4}"));
                % tex_axis_yellow = app.pl.declare("tex_axis_yellow", app.pl.texture('pigment', [1.0 1.0 0.0], 'finish', "phong 1 reflection {0.10 metallic 0.4}"));
                if(app.axis_enable.Value)
                    app.pl.axis('length', [5 5 5]);
                end

                % TODO: Move textures to materials.inc
                % tex_grid_gray  = app.pl.declare("tex_grid_gray",  app.pl.texture('pigment', [0.5 0.5 0.5], 'finish', "phong 1 reflection {0.10 metallic 0.4}"));
                % tex_grid_red   = app.pl.declare("tex_grid_red",   app.pl.texture('pigment', [1.0 0.0 0.0], 'finish', "phong 1 reflection {0.10 metallic 0.4}"));
                % tex_grid_green = app.pl.declare("tex_grid_green", app.pl.texture('pigment', [0.0 1.0 0.0], 'finish', "phong 1 reflection {0.10 metallic 0.4}"));
                % tex_grid_blue  = app.pl.declare("tex_grid_blue",  app.pl.texture('pigment', [0.0 0.0 1.0], 'finish', "phong 1 reflection {0.10 metallic 0.4}"));
                % pl.grid('width', 10, 'height', 10);pl.grid_2D('width', 10, 'height', 10);pl.grid_2D('width', 10, 'height', 10);
                % Grids
                if(app.grid_xz_enable.Value)
                    app.pl.grid('width', 30, 'height', 30, 'radius', 0.03, 'scale', [0.2, 0.2 0.2]);
                end
                if(app.grid_yz_enable.Value)
                    app.pl.grid('width', 30, 'height', 30, 'radius', 0.03, 'rotate', [0 0 90], 'scale', [0.2, 0.2 0.2]);
                end
                if(app.grid_xy_enable.Value)
                    app.pl.grid('width', 30, 'height', 30, 'radius', 0.03, 'rotate', [90 0 0], 'scale', [0.2, 0.2 0.2]);
                end

                % Planes
                if(app.plane_xz_enable.Value)
                    app.pl.plane('normal', [0 1 0]);
                end
                if(app.plane_yz_enable.Value)

                    app.pl.plane('normal', [1 0 0]);
                end
                if(app.plane_xy_enable.Value)
                    app.pl.plane('normal', [0 0 1]);
                end
                
            app.pl.include_end();
        end
        
        % Create scene
        function create_scene(app)
            % TODO: 1. Show notification
            %       2. Prepare scen before rendering on demand (flag of changes) 
            app.scene_time.Text = sprintf("Preparing ...");
            drawnow();
            tic % Time measure            
            app.pl.scene_begin('scene_file', 'surface.pl', 'image_file', 'surface.png');
                app.pl.include("camera");
                app.pl.include("lights");
                app.pl.include("helpers");
                app.pl.include("textures");
                
                app.pl.global_settings("assumed_gamma 1");
               
                [X,Y,Z] = peaks(app.surf_size.Value);
                % Hide figure
                f = figure('Visible', 'off');
                ax = axes(f);
                s = surf(ax,X,Y,Z);
                
                if(app.surf_texture_odd.Value == "None")
                    texture_odd = "";
                else
                    texture_odd = app.surf_texture_odd.Value;
                end

                if(app.surf_texture_even.Value == "None")
                    texture_even = "";
                else
                    texture_even = app.surf_texture_even.Value;
                end

                app.pl.surface('surface', s, 'smooth', app.surf_smooth.Value, 'colormap', app.surf_colormap.Value, 'scale', [1, 1, 3/10],...
                                'texture_odd', texture_odd, 'texture_even', texture_even);
            app.pl.scene_end();
            app.scene_time.Text = sprintf("Prepared in %0.2f sec", toc);
        end

        % Render
        function render(app)
            app.render_time.Text = sprintf("Rendering ...");
            drawnow();
            tic % Measure rendering            
            % Render
            img_file = app.pl.render();
            img = imread(img_file);
            app.image.ImageSource = img;
            drawnow();
            app.render_time.Text = sprintf("Rendered in: %0.2f sec", toc);
        end
    end   

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startup(app)
            app.cam_location = app.cam_location_default;
            app.cam_loc_x.Value = app.cam_location(1);
            app.cam_loc_y.Value = app.cam_location(2);
            app.cam_loc_z.Value = app.cam_location(3);

            app.cam_look_at = app.cam_look_at_default;
            app.cam_look_at_x.Value = app.cam_look_at(1);
            app.cam_look_at_y.Value = app.cam_look_at(2);
            app.cam_look_at_z.Value = app.cam_look_at(3);

            init_pov(app);
            create_camera(app);
            create_lights(app);
            create_helpers(app);
            create_scene(app);
            render(app);
        end

        % Button pushed function: bt_render
        function on_bt_render(app, event)
%             app.cam_location = [app.cam_loc_x.Value app.cam_loc_y.Value app.cam_loc_z.Value];
%             app.cam_look_at = [app.cam_look_at_x.Value app.cam_look_at_y.Value app.cam_look_at_z.Value];
%             create_camera(app);
%             create_lights(app);
%             create_helpers(app);
            if(app.need_create_scene)
                create_scene(app);
                app.need_create_scene = false;
            end
            render(app);
        end

        % Value changed function: cam_loc_x
        function on_cam_loc_x(app, event)
            app.cam_location(1) = app.cam_loc_x.Value;
            create_camera(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Value changed function: cam_loc_y
        function on_cam_loc_y(app, event)
            app.cam_location(2) = app.cam_loc_y.Value;
            create_camera(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Value changed function: cam_loc_z
        function on_cam_loc_z(app, event)
            app.cam_location(3) = app.cam_loc_z.Value;
            create_camera(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Value changed function: cam_look_at_x
        function on_cam_look_at_x(app, event)
            app.cam_look_at(1) = app.cam_look_at_x.Value;
            create_camera(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Value changed function: cam_look_at_y
        function on_cam_look_at_y(app, event)
            app.cam_look_at(2) = app.cam_look_at_y.Value;
            create_camera(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Value changed function: cam_look_at_z
        function on_cam_look_at_z(app, event)
            app.cam_look_at(3) = app.cam_look_at_z.Value;
            create_camera(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Value changed function: cam_type
        function on_cam_type(app, event)
            create_camera(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Value changed function: cam_angle
        function on_cam_angle(app, event)
            create_camera(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Value changed function: axis_enable, grid_xy_enable, 
        % grid_xz_enable, grid_yz_enable, plane_xy_enable, 
        % plane_xz_enable, plane_yz_enable
        function on_helper_changed(app, event)
            create_helpers(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end

        % Button pushed function: ResetButton
        function on_cam_reset(app, event)
            app.cam_location = app.cam_location_default;
            app.cam_look_at  = app.cam_look_at_default;
            
            app.cam_look_at_x.Value = app.cam_look_at(1);
            app.cam_look_at_y.Value = app.cam_look_at(2);
            app.cam_look_at_z.Value = app.cam_look_at(3);

            app.cam_loc_x.Value = app.cam_location(1);
            app.cam_loc_y.Value = app.cam_location(2);
            app.cam_loc_z.Value = app.cam_location(3);
            app.cam_type.Value  = "perspective";
            app.cam_angle.Value = 45;

            create_camera(app);
            render(app);
        end

        % Value changed function: surf_colormap, surf_size, surf_smooth, 
        % surf_texture_even, surf_texture_odd
        function on_surf_changed(app, event)
            if(app.cam_responsive.Value)
                create_scene(app);                
                render(app);
            else
                app.need_create_scene = true;
            end
        end

        % Value changed function: light_1, light_2, light_3, light_4
        function on_light_changed(app, event)
            create_lights(app);
            if(app.cam_responsive.Value)
                render(app);
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ui_figure and hide until all components are created
            app.ui_figure = uifigure('Visible', 'off');
            app.ui_figure.AutoResizeChildren = 'off';
            app.ui_figure.Color = [0.651 0.651 0.651];
            app.ui_figure.Position = [100 100 1137 598];
            app.ui_figure.Name = 'POV-Lab / Surface demo';
            app.ui_figure.Resize = 'off';

            % Create panel
            app.panel = uipanel(app.ui_figure);
            app.panel.AutoResizeChildren = 'off';
            app.panel.TitlePosition = 'centertop';
            app.panel.Title = 'Location';
            app.panel.BackgroundColor = [0.651 0.651 0.651];
            app.panel.Position = [1003 321 121 117];

            % Create xLabel_4
            app.xLabel_4 = uilabel(app.panel);
            app.xLabel_4.HorizontalAlignment = 'right';
            app.xLabel_4.Position = [15 61 25 22];
            app.xLabel_4.Text = 'X';

            % Create cam_loc_x
            app.cam_loc_x = uispinner(app.panel);
            app.cam_loc_x.ValueChangedFcn = createCallbackFcn(app, @on_cam_loc_x, true);
            app.cam_loc_x.HorizontalAlignment = 'center';
            app.cam_loc_x.Position = [55 61 50 22];

            % Create xLabel_5
            app.xLabel_5 = uilabel(app.panel);
            app.xLabel_5.HorizontalAlignment = 'right';
            app.xLabel_5.Position = [15 40 25 22];
            app.xLabel_5.Text = 'Y';

            % Create cam_loc_y
            app.cam_loc_y = uispinner(app.panel);
            app.cam_loc_y.ValueChangedFcn = createCallbackFcn(app, @on_cam_loc_y, true);
            app.cam_loc_y.HorizontalAlignment = 'center';
            app.cam_loc_y.Position = [55 40 50 22];

            % Create xLabel_6
            app.xLabel_6 = uilabel(app.panel);
            app.xLabel_6.HorizontalAlignment = 'right';
            app.xLabel_6.Position = [15 19 25 22];
            app.xLabel_6.Text = 'Z';

            % Create cam_loc_z
            app.cam_loc_z = uispinner(app.panel);
            app.cam_loc_z.ValueChangedFcn = createCallbackFcn(app, @on_cam_loc_z, true);
            app.cam_loc_z.HorizontalAlignment = 'center';
            app.cam_loc_z.Position = [55 19 50 22];

            % Create panel_2
            app.panel_2 = uipanel(app.ui_figure);
            app.panel_2.AutoResizeChildren = 'off';
            app.panel_2.TitlePosition = 'centertop';
            app.panel_2.Title = 'Look at';
            app.panel_2.BackgroundColor = [0.651 0.651 0.651];
            app.panel_2.Position = [1004 195 119 112];

            % Create xLabel
            app.xLabel = uilabel(app.panel_2);
            app.xLabel.HorizontalAlignment = 'right';
            app.xLabel.Position = [15 57 25 22];
            app.xLabel.Text = 'X';

            % Create cam_look_at_x
            app.cam_look_at_x = uispinner(app.panel_2);
            app.cam_look_at_x.ValueChangedFcn = createCallbackFcn(app, @on_cam_look_at_x, true);
            app.cam_look_at_x.HorizontalAlignment = 'center';
            app.cam_look_at_x.Position = [55 57 50 22];

            % Create xLabel_2
            app.xLabel_2 = uilabel(app.panel_2);
            app.xLabel_2.HorizontalAlignment = 'right';
            app.xLabel_2.Position = [15 36 25 22];
            app.xLabel_2.Text = 'Y';

            % Create cam_look_at_y
            app.cam_look_at_y = uispinner(app.panel_2);
            app.cam_look_at_y.ValueChangedFcn = createCallbackFcn(app, @on_cam_look_at_y, true);
            app.cam_look_at_y.HorizontalAlignment = 'center';
            app.cam_look_at_y.Position = [55 36 50 22];

            % Create xLabel_3
            app.xLabel_3 = uilabel(app.panel_2);
            app.xLabel_3.HorizontalAlignment = 'right';
            app.xLabel_3.Position = [15 15 25 22];
            app.xLabel_3.Text = 'Z';

            % Create cam_look_at_z
            app.cam_look_at_z = uispinner(app.panel_2);
            app.cam_look_at_z.ValueChangedFcn = createCallbackFcn(app, @on_cam_look_at_z, true);
            app.cam_look_at_z.HorizontalAlignment = 'center';
            app.cam_look_at_z.Position = [55 15 50 22];

            % Create panel_3
            app.panel_3 = uipanel(app.ui_figure);
            app.panel_3.AutoResizeChildren = 'off';
            app.panel_3.TitlePosition = 'centertop';
            app.panel_3.Title = 'Camera';
            app.panel_3.BackgroundColor = [0.651 0.651 0.651];
            app.panel_3.Position = [1002 448 122 137];

            % Create AngleSpinnerLabel
            app.AngleSpinnerLabel = uilabel(app.panel_3);
            app.AngleSpinnerLabel.HorizontalAlignment = 'right';
            app.AngleSpinnerLabel.Position = [8 47 46 22];
            app.AngleSpinnerLabel.Text = 'Angle';

            % Create cam_angle
            app.cam_angle = uispinner(app.panel_3);
            app.cam_angle.ValueChangedFcn = createCallbackFcn(app, @on_cam_angle, true);
            app.cam_angle.Position = [66 47 45 22];
            app.cam_angle.Value = 48;

            % Create cam_type
            app.cam_type = uidropdown(app.panel_3);
            app.cam_type.Items = {'orthographic', 'perspective', 'fisheye', 'ultra_wide_angle', 'spherical', 'omnimax', 'panoramic', 'cylinder 1', 'cylinder 2', 'cylinder 3', 'cylinder 4'};
            app.cam_type.ValueChangedFcn = createCallbackFcn(app, @on_cam_type, true);
            app.cam_type.Position = [5 9 114 22];
            app.cam_type.Value = 'perspective';

            % Create ResetButton
            app.ResetButton = uibutton(app.panel_3, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @on_cam_reset, true);
            app.ResetButton.Position = [12 81 100 22];
            app.ResetButton.Text = 'Reset';

            % Create image_panel
            app.image_panel = uipanel(app.ui_figure);
            app.image_panel.AutoResizeChildren = 'off';
            app.image_panel.BackgroundColor = [0.651 0.651 0.651];
            app.image_panel.Position = [198 12 792 573];

            % Create image
            app.image = uiimage(app.image_panel);
            app.image.Position = [10 10 771 553];

            % Create panel_4
            app.panel_4 = uipanel(app.ui_figure);
            app.panel_4.AutoResizeChildren = 'off';
            app.panel_4.TitlePosition = 'centertop';
            app.panel_4.Title = 'Helpers';
            app.panel_4.BackgroundColor = [0.651 0.651 0.651];
            app.panel_4.Position = [15 448 174 137];

            % Create axis_enable
            app.axis_enable = uicheckbox(app.panel_4);
            app.axis_enable.ValueChangedFcn = createCallbackFcn(app, @on_helper_changed, true);
            app.axis_enable.Text = 'Axis';
            app.axis_enable.Position = [65 81 45 22];
            app.axis_enable.Value = true;

            % Create grid_xz_enable
            app.grid_xz_enable = uicheckbox(app.panel_4);
            app.grid_xz_enable.ValueChangedFcn = createCallbackFcn(app, @on_helper_changed, true);
            app.grid_xz_enable.Text = 'Grid XZ';
            app.grid_xz_enable.Position = [9 57 63 22];

            % Create grid_yz_enable
            app.grid_yz_enable = uicheckbox(app.panel_4);
            app.grid_yz_enable.ValueChangedFcn = createCallbackFcn(app, @on_helper_changed, true);
            app.grid_yz_enable.Text = 'Grid YZ';
            app.grid_yz_enable.Position = [10 36 63 22];

            % Create grid_xy_enable
            app.grid_xy_enable = uicheckbox(app.panel_4);
            app.grid_xy_enable.ValueChangedFcn = createCallbackFcn(app, @on_helper_changed, true);
            app.grid_xy_enable.Text = 'Grid XY';
            app.grid_xy_enable.Position = [11 15 64 22];

            % Create plane_xz_enable
            app.plane_xz_enable = uicheckbox(app.panel_4);
            app.plane_xz_enable.ValueChangedFcn = createCallbackFcn(app, @on_helper_changed, true);
            app.plane_xz_enable.Text = 'Plane XZ';
            app.plane_xz_enable.Position = [84 57 71 22];

            % Create plane_yz_enable
            app.plane_yz_enable = uicheckbox(app.panel_4);
            app.plane_yz_enable.ValueChangedFcn = createCallbackFcn(app, @on_helper_changed, true);
            app.plane_yz_enable.Text = 'Plane YZ';
            app.plane_yz_enable.Position = [85 36 71 22];

            % Create plane_xy_enable
            app.plane_xy_enable = uicheckbox(app.panel_4);
            app.plane_xy_enable.ValueChangedFcn = createCallbackFcn(app, @on_helper_changed, true);
            app.plane_xy_enable.Text = 'Plane XY';
            app.plane_xy_enable.Position = [86 15 72 22];

            % Create panel_5
            app.panel_5 = uipanel(app.ui_figure);
            app.panel_5.AutoResizeChildren = 'off';
            app.panel_5.TitlePosition = 'centertop';
            app.panel_5.Title = 'Lights';
            app.panel_5.BackgroundColor = [0.651 0.651 0.651];
            app.panel_5.Position = [16 340 173 100];

            % Create light_1
            app.light_1 = uicheckbox(app.panel_5);
            app.light_1.ValueChangedFcn = createCallbackFcn(app, @on_light_changed, true);
            app.light_1.Text = 'LIght 1';
            app.light_1.Position = [11 42 59 22];
            app.light_1.Value = true;

            % Create light_2
            app.light_2 = uicheckbox(app.panel_5);
            app.light_2.ValueChangedFcn = createCallbackFcn(app, @on_light_changed, true);
            app.light_2.Text = 'LIght 2';
            app.light_2.Position = [89 42 59 22];

            % Create light_3
            app.light_3 = uicheckbox(app.panel_5);
            app.light_3.ValueChangedFcn = createCallbackFcn(app, @on_light_changed, true);
            app.light_3.Text = 'LIght 3';
            app.light_3.Position = [13 8 59 22];

            % Create light_4
            app.light_4 = uicheckbox(app.panel_5);
            app.light_4.ValueChangedFcn = createCallbackFcn(app, @on_light_changed, true);
            app.light_4.Text = 'Light 4';
            app.light_4.Position = [90 8 58 22];

            % Create panel_6
            app.panel_6 = uipanel(app.ui_figure);
            app.panel_6.AutoResizeChildren = 'off';
            app.panel_6.TitlePosition = 'centertop';
            app.panel_6.Title = 'Rendering';
            app.panel_6.BackgroundColor = [0.651 0.651 0.651];
            app.panel_6.Position = [1002 45 122 137];

            % Create cam_responsive
            app.cam_responsive = uicheckbox(app.panel_6);
            app.cam_responsive.Text = 'Responsive';
            app.cam_responsive.Position = [19 82 85 22];
            app.cam_responsive.Value = true;

            % Create bt_render
            app.bt_render = uibutton(app.panel_6, 'push');
            app.bt_render.ButtonPushedFcn = createCallbackFcn(app, @on_bt_render, true);
            app.bt_render.Position = [12 19 100 55];
            app.bt_render.Text = 'Render';

            % Create panel_7
            app.panel_7 = uipanel(app.ui_figure);
            app.panel_7.AutoResizeChildren = 'off';
            app.panel_7.TitlePosition = 'centertop';
            app.panel_7.Title = 'Surface';
            app.panel_7.BackgroundColor = [0.651 0.651 0.651];
            app.panel_7.Position = [19 45 171 286];

            % Create surf_smooth
            app.surf_smooth = uicheckbox(app.panel_7);
            app.surf_smooth.ValueChangedFcn = createCallbackFcn(app, @on_surf_changed, true);
            app.surf_smooth.Text = 'Smooth';
            app.surf_smooth.Position = [25 61 63 22];
            app.surf_smooth.Value = true;

            % Create OddDropDownLabel
            app.OddDropDownLabel = uilabel(app.panel_7);
            app.OddDropDownLabel.HorizontalAlignment = 'right';
            app.OddDropDownLabel.Position = [20 219 23 22];
            app.OddDropDownLabel.Text = 'Odd';

            % Create surf_texture_odd
            app.surf_texture_odd = uidropdown(app.panel_7);
            app.surf_texture_odd.Items = {'None', 'DMFWood6', 'NBglass', 'NBoldglass', 'NBwinebottle', 'NBbeerbottle', 'Ruby_Glass', 'Dark_Green_Glass', 'Yellow_Glass', 'Orange_Glass', 'Vicks_Bottle_Glass', 'Soft_Silver', 'New_Penny', 'Tinny_Brass', 'Gold_Nugget', 'Aluminum', 'Bright_Bronze'};
            app.surf_texture_odd.ValueChangedFcn = createCallbackFcn(app, @on_surf_changed, true);
            app.surf_texture_odd.Position = [52 219 107 22];
            app.surf_texture_odd.Value = 'None';

            % Create EvenDropDownLabel
            app.EvenDropDownLabel = uilabel(app.panel_7);
            app.EvenDropDownLabel.HorizontalAlignment = 'right';
            app.EvenDropDownLabel.Position = [10 184 33 22];
            app.EvenDropDownLabel.Text = 'Even';

            % Create surf_texture_even
            app.surf_texture_even = uidropdown(app.panel_7);
            app.surf_texture_even.Items = {'None', 'DMFWood6', 'NBglass', 'NBoldglass', 'NBwinebottle', 'NBbeerbottle', 'Ruby_Glass', 'Dark_Green_Glass', 'Yellow_Glass', 'Orange_Glass', 'Vicks_Bottle_Glass', 'Soft_Silver', 'New_Penny', 'Tinny_Brass', 'Gold_Nugget', 'Aluminum', 'Bright_Bronze'};
            app.surf_texture_even.ValueChangedFcn = createCallbackFcn(app, @on_surf_changed, true);
            app.surf_texture_even.Position = [52 184 107 22];
            app.surf_texture_even.Value = 'None';

            % Create SizeSpinnerLabel
            app.SizeSpinnerLabel = uilabel(app.panel_7);
            app.SizeSpinnerLabel.HorizontalAlignment = 'right';
            app.SizeSpinnerLabel.Position = [13 150 29 22];
            app.SizeSpinnerLabel.Text = 'Size';

            % Create surf_size
            app.surf_size = uispinner(app.panel_7);
            app.surf_size.Limits = [2 Inf];
            app.surf_size.ValueChangedFcn = createCallbackFcn(app, @on_surf_changed, true);
            app.surf_size.Position = [52 150 107 22];
            app.surf_size.Value = 20;

            % Create ColormapDropDownLabel
            app.ColormapDropDownLabel = uilabel(app.panel_7);
            app.ColormapDropDownLabel.HorizontalAlignment = 'right';
            app.ColormapDropDownLabel.Position = [-4 115 65 22];
            app.ColormapDropDownLabel.Text = 'Colormap';

            % Create surf_colormap
            app.surf_colormap = uidropdown(app.panel_7);
            app.surf_colormap.Items = {'parula', 'turbo', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'jet', 'lines', 'colorcube', 'prism', 'flag', 'white'};
            app.surf_colormap.ValueChangedFcn = createCallbackFcn(app, @on_surf_changed, true);
            app.surf_colormap.Position = [64 115 94 22];
            app.surf_colormap.Value = 'parula';

            % Create scene_time
            app.scene_time = uilabel(app.ui_figure);
            app.scene_time.Position = [17 12 172 23];
            app.scene_time.Text = 'Scene created in';

            % Create render_time
            app.render_time = uilabel(app.ui_figure);
            app.render_time.Position = [1002 12 123 22];
            app.render_time.Text = 'Rendered in';

            % Show the figure after all components are created
            app.ui_figure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app_surface_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ui_figure)

            % Execute the startup function
            runStartupFcn(app, @startup)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ui_figure)
        end
    end
end