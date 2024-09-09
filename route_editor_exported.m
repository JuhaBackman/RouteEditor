%
%  This file is part of ISOBUS Guidance Path editor
%
%  Copyright 2024 Juha Backman / Natural Resources Institute Finland
%
%  ISOBUS Guidance Path editor is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as
%  published by the Free Software Foundation, either version 3 of
%  the License, or (at your option) any later version.
%
%  ISOBUS Guidance Path editor is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public
%  License along with EFDIClientFMIS.
%  If not, see <http://www.gnu.org/licenses/>.
%


classdef route_editor_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ISOBUSGuidancePatheditorUIFigure  matlab.ui.Figure
        FileMenu               matlab.ui.container.Menu
        LoadMenu               matlab.ui.container.Menu
        SaveMenu               matlab.ui.container.Menu
        GridLayout             matlab.ui.container.GridLayout
        TabGroup               matlab.ui.container.TabGroup
        MAPTab                 matlab.ui.container.Tab
        GridLayout4            matlab.ui.container.GridLayout
        GridLayout2            matlab.ui.container.GridLayout
        ExtendButton           matlab.ui.control.Button
        CopyButton             matlab.ui.control.Button
        InsertButton           matlab.ui.control.Button
        ModifypathButton       matlab.ui.control.Button
        PointpropertyButton    matlab.ui.control.Button
        EDITButtonGroup        matlab.ui.container.ButtonGroup
        ObstacleButton         matlab.ui.control.RadioButton
        GuidancePatternButton  matlab.ui.control.RadioButton
        ConnectButton          matlab.ui.control.Button
        DeleteButton           matlab.ui.control.Button
        UIAxes                 matlab.ui.control.UIAxes
        GraphTab               matlab.ui.container.Tab
        GridLayout5            matlab.ui.container.GridLayout
        UIAxes4                matlab.ui.control.UIAxes
        UIAxes3                matlab.ui.control.UIAxes
        UIAxes2                matlab.ui.control.UIAxes
        GridLayout3            matlab.ui.container.GridLayout
        MessageLabel           matlab.ui.control.Label
        READYButton            matlab.ui.control.Button
    end


    properties (Access = private)
        mode

        selectIdxObstacle
        selectGuidancePattern
        selectPartfieldBoundary

        newObstacle
        newGuidancePattern

        constants
        defaults
    end

    properties (Access = public)
        Origin
        Obstacle
        GuidancePattern
        PartfieldBoundary
    end

    methods (Access = private)

        function UpdateFigures(app, resetView)

            if ~exist('resetView','var')
                resetView = false;
            end

            % Update MAP
            cla(app.UIAxes);
            hold(app.UIAxes, 'on');

            for i=1:numel(app.Obstacle)
                plot(app.UIAxes, app.Obstacle(i).x, app.Obstacle(i).y,'ko','ButtonDownFcn',@mouse_click);

                if(app.selectIdxObstacle(i))
                    plot(app.UIAxes, app.Obstacle(i).x, app.Obstacle(i).y,'rx','ButtonDownFcn',@mouse_click);
                end
            end

            for i=1:numel(app.newObstacle)
                plot(app.UIAxes, app.newObstacle(i).x, app.newObstacle(i).y,'ro');
            end

            for i=1:numel(app.GuidancePattern)
                plot(app.UIAxes, app.GuidancePattern(i).x, app.GuidancePattern(i).y,'b','ButtonDownFcn',@mouse_click)

                % TODO: piirtää väärin, jos työalue ei ole yhtenäinen
                plot(app.UIAxes, app.GuidancePattern(i).x(app.GuidancePattern(i).implement(:) == 1), app.GuidancePattern(i).y(app.GuidancePattern(i).implement(:) == 1),'g','ButtonDownFcn',@mouse_click)
                plot(app.UIAxes, app.GuidancePattern(i).x(end), app.GuidancePattern(i).y(end),'bx','ButtonDownFcn',@mouse_click)

                plot(app.UIAxes, app.GuidancePattern(i).x(app.selectGuidancePattern(i).Idx), app.GuidancePattern(i).y(app.selectGuidancePattern(i).Idx),'rx','ButtonDownFcn',@mouse_click)
            end

            plot(app.UIAxes, app.newGuidancePattern.x, app.newGuidancePattern.y,'r');

            for i=1:numel(app.PartfieldBoundary)
                plot(app.UIAxes, app.PartfieldBoundary(i).x, app.PartfieldBoundary(i).y,'k','ButtonDownFcn',@mouse_click)
                plot(app.UIAxes, app.PartfieldBoundary(i).x(app.selectPartfieldBoundary(i).Idx), app.PartfieldBoundary(i).y(app.selectPartfieldBoundary(i).Idx),'rx','ButtonDownFcn',@mouse_click)
            end

            if(resetView)
                axis(app.UIAxes, 'equal');
            end

            function mouse_click(~,eventData)
                UIAxesButtonDown(app,eventData);
            end


            cla(app.UIAxes2);
            hold(app.UIAxes2, 'on');

            cla(app.UIAxes3);
            hold(app.UIAxes3, 'on');

            cla(app.UIAxes4);
            hold(app.UIAxes4, 'on');

            if numel(app.GuidancePattern) > 0
                % Generate time
                dx = diff( [app.GuidancePattern(:).x] );
                dy = diff( [app.GuidancePattern(:).y] );
                v = [app.GuidancePattern(:).speed];
                v = (v(1:end-1) + v(2:end)) / 2;
                dx(v == 0) = 0;
                dy(v == 0) = 0;
                v(v == 0) = 1;
                t = cumsum([0 sqrt(dx.^2 + dy.^2) ./ abs(v)]);

                time = {};
                k = 1;
                for i=1:numel(app.GuidancePattern)
                    time{i} = t(k:(k+numel(app.GuidancePattern(i).speed)-1));
                    k = k + numel(app.GuidancePattern(i).speed);
                end

                % Plot speed
                for i=1:numel(app.GuidancePattern)
                    plot(app.UIAxes2, time{i}, app.GuidancePattern(i).speed,'b');
                    plot(app.UIAxes2, [time{i}(1) time{i}(1)], [0 2],'r');
                end
                plot(app.UIAxes2, [time{end}(end) time{end}(end)], [0 2],'r');

                % Plot heading
                for i=1:numel(app.GuidancePattern)
                    plot(app.UIAxes3, time{i}, app.GuidancePattern(i).yaw,'b');
                    plot(app.UIAxes3, [time{i}(1) time{i}(1)], [0 2],'r');
                end
                plot(app.UIAxes3, [time{end}(end) time{end}(end)], [0 2],'r');

                % Plot state
                for i=1:numel(app.GuidancePattern)
                    plot(app.UIAxes4, time{i}, app.GuidancePattern(i).implement,'b');
                    plot(app.UIAxes4, [time{i}(1) time{i}(1)], [0 2],'r');
                end
                plot(app.UIAxes4, [time{end}(end) time{end}(end)], [0 2],'r');

                linkaxes([app.UIAxes2 app.UIAxes3 app.UIAxes4],'x');
                axis(app.UIAxes2, "auto");
                axis(app.UIAxes3, "auto");
                axis(app.UIAxes4, "auto");
            end
        end

        function message(app, text)

            if(isstring(text) == false)
                text = matlab.unittest.diagnostics.ConstraintDiagnostic.getDisplayableString(text);
            end

            if(strlength(text) == 0)
                app.MessageLabel.Text = "";
            else
                app.MessageLabel.Text = app.MessageLabel.Text + compose("\n") + text;
                disp(text);
            end
        end

        function resetSelection(app)
            app.selectIdxObstacle = false(1,numel(app.Obstacle));
            app.selectGuidancePattern = struct;
            for i = 1:numel(app.GuidancePattern)
                app.selectGuidancePattern(i).Idx = false(1,numel(app.GuidancePattern(i).x));
            end

            app.selectPartfieldBoundary = struct;
            for i = 1:numel(app.PartfieldBoundary)
                app.selectPartfieldBoundary(i).Idx = false(1,numel(app.PartfieldBoundary(i).x));
            end

            app.newObstacle = [];
            app.newGuidancePattern = struct;
            app.newGuidancePattern.x = [];
            app.newGuidancePattern.y = [];
            app.newGuidancePattern.yaw = [];
            app.newGuidancePattern.speed = [];
            app.newGuidancePattern.implement = [];
        end

        function setMode(app, mode)
            app.mode = mode;

            if(strlength(mode) > 0)
                app.READYButton.Enable = "on";
                app.ObstacleButton.Enable = "off";
                app.GuidancePatternButton.Enable = "off";
            else
                app.READYButton.Enable = "off";
                app.ObstacleButton.Enable = "on";
                app.GuidancePatternButton.Enable = "on";
            end
        end

        function modifyPathPoints(app,sideshift,speed,workstate)
            % set new speed
            if exist('speed','var') && ~isnan(speed)
                for i=1:numel(app.GuidancePattern)
                    app.GuidancePattern(i).speed(app.selectGuidancePattern(i).Idx(:)) = speed;
                end
            end

            % set new workstate
            if exist('workstate','var') && ~isnan(workstate)
                for i=1:numel(app.GuidancePattern)
                    app.GuidancePattern(i).implement(app.selectGuidancePattern(i).Idx(:)) = workstate;
                end
            end

            % move path sideways
            for i=1:numel(app.GuidancePattern)
                sideshift2 = sideshift;
                while(abs(sideshift2) > 0 && any(app.selectGuidancePattern(i).Idx))
                    dx = app.GuidancePattern(i).x(2:end) - app.GuidancePattern(i).x(1:(end-1));
                    dy = app.GuidancePattern(i).y(2:end) - app.GuidancePattern(i).y(1:(end-1));

                    length = sqrt(dx.*dx + dy.*dy);

                    dx = dx ./ length;
                    dy = dy ./ length;

                    dx(end+1) = dx(end);
                    dy(end+1) = dy(end);

                    dx(2:end-1) = (dx(1:end-2) + dx(2:end-1)) / 2;
                    dy(2:end-1) = (dy(1:end-2) + dy(2:end-1)) / 2;

                    length = min(length);
                    step = sign(sideshift2)*max(abs(sideshift2), 0.5*length);
                    sideshift2 = sideshift2 - step;

                    app.GuidancePattern(i).x(app.selectGuidancePattern(i).Idx(:)) = app.GuidancePattern(i).x(app.selectGuidancePattern(i).Idx(:)) + step * dy(app.selectGuidancePattern(i).Idx(:));
                    app.GuidancePattern(i).y(app.selectGuidancePattern(i).Idx(:)) = app.GuidancePattern(i).y(app.selectGuidancePattern(i).Idx(:)) - step * dx(app.selectGuidancePattern(i).Idx(:));

                    dx = app.GuidancePattern(i).x(2:end) - app.GuidancePattern(i).x(1:(end-1));
                    dy = app.GuidancePattern(i).y(2:end) - app.GuidancePattern(i).y(1:(end-1));

                    length = sqrt(dx.*dx + dy.*dy);
                    newIdx = true(1,numel(app.GuidancePattern(i).x));

                    for j=1:numel(length)
                        if(length(j) < 0.5)
                            app.GuidancePattern(i).x(j) = (app.GuidancePattern(i).x(j) + app.GuidancePattern(i).x(j+1)) / 2;
                            app.GuidancePattern(i).y(j) = (app.GuidancePattern(i).y(j) + app.GuidancePattern(i).y(j+1)) / 2;
                            app.GuidancePattern(i).speed(j) = (app.GuidancePattern(i).speed(j) + app.GuidancePattern(i).speed(j+1)) / 2;
                            app.GuidancePattern(i).implement(j) = app.GuidancePattern(i).implement(j) == 1 || app.GuidancePattern(i).implement(j+1) == 1;

                            app.selectGuidancePattern(i).Idx(j) = true;

                            newIdx(j+1) = false;
                        end
                    end

                    for j=3:numel(newIdx)
                        if(all(newIdx(j-2:j) == false))
                            newIdx(j-1) = true;
                        end
                    end

                    app.GuidancePattern(i).x = app.GuidancePattern(i).x(newIdx);
                    app.GuidancePattern(i).y = app.GuidancePattern(i).y(newIdx);
                    app.GuidancePattern(i).yaw = app.GuidancePattern(i).yaw(newIdx);
                    app.GuidancePattern(i).speed = app.GuidancePattern(i).speed(newIdx);
                    app.GuidancePattern(i).implement = app.GuidancePattern(i).implement(newIdx);

                    app.selectGuidancePattern(i).Idx = app.selectGuidancePattern(i).Idx(newIdx);
                end
            end

            % Recalculate yaw angle
            for i=1:numel(app.GuidancePattern)
                dx = app.GuidancePattern(i).x(2:end) - app.GuidancePattern(i).x(1:(end-1));
                dy = app.GuidancePattern(i).y(2:end) - app.GuidancePattern(i).y(1:(end-1));

                length = sqrt(dx.*dx + dy.*dy);

                dx = dx ./ length;
                dy = dy ./ length;

                dx(end+1) = dx(end);
                dy(end+1) = dy(end);

                dx = dx .* sign(app.GuidancePattern(i).speed);
                dy = dy .* sign(app.GuidancePattern(i).speed);

                dx(2:end-1) = (dx(1:end-2) + dx(2:end-1)) / 2;
                dy(2:end-1) = (dy(1:end-2) + dy(2:end-1)) / 2;

                yaw = atan2(dy,dx);

                for j=1:numel(yaw)
                    if(isnan(yaw(j)))
                        yaw(j) = yaw(j) - 1;
                    end
                end

                app.GuidancePattern(i).yaw(app.selectGuidancePattern(i).Idx(:)) = yaw(app.selectGuidancePattern(i).Idx(:));
            end
        end


        function connectPaths(app, point0path, point0index, point1path, point1index)

            message(app,"Give path properties. Possible turting types are: 'FASTEST','FASTEST_FORWARDS','FASTEST_BACKWARDS','LRL', 'RLR', 'LSL', 'RSR', 'LSR', 'RSL', 'L|R|L', 'R|L|R'");

            DefaultValues = {};
            VariableNames = {};

            DefaultValues{1} = app.defaults.turningType;
            VariableNames{1} = 'turning type';

            DefaultValues{2} = num2str(app.constants.dt);
            VariableNames{2} = 'calculation resolution';

            DefaultValues{3} = num2str(1/app.constants.K_max);
            VariableNames{3} = 'minimum turning radius';

            DefaultValues{4} = num2str((2*app.constants.K_max)/app.constants.dK_max);
            VariableNames{4} = 'time from side to side';

            DefaultValues{5} = num2str(app.constants.dv_max);
            VariableNames{5} = 'max acceleration';

            DefaultValues{6} = num2str(app.constants.dv_min);
            VariableNames{6} = 'max deceleration';

            DefaultValues{7} = num2str(app.constants.v_start_sequence);
            VariableNames{7} = 'turning speed profile at start';

            DefaultValues{8} = num2str(app.constants.v_end_sequence);
            VariableNames{8} = 'turning speed profile at end';

            DefaultValues{9} = num2str(app.constants.v_center_reverse);
            VariableNames{9} = 'reverse speed';

            answer = inputdlg(VariableNames,"Give turning information",1,DefaultValues);
            if(numel(answer) == 0)
                message(app,"User cancelled path calculation");
                return;
            end

            app.defaults.turningType = answer{1};

            app.constants.dt = str2double(answer{2});

            R_min = str2double(answer{3});
            Turn_time = str2double(answer{4});

            app.constants.K_max = 1/R_min;
            app.constants.K_min = -1/R_min;

            app.constants.dK_max = (2*app.constants.K_max)/Turn_time;
            app.constants.dK_min = (2*app.constants.K_min)/Turn_time;

            app.constants.dv_max = str2double(answer{5});
            app.constants.dv_min = str2double(answer{6});

            app.constants.v_start_sequence = str2double(split(answer{7}));
            app.constants.v_end_sequence = str2double(split(answer{8}));

            app.constants.v_center_reverse = str2double(answer{9});


            x0 = app.GuidancePattern(point0path).x(point0index);
            y0 = app.GuidancePattern(point0path).y(point0index);
            yaw0 = app.GuidancePattern(point0path).yaw(point0index);
            v0 = app.GuidancePattern(point0path).speed(point0index);
            k0 = 0;

            x1 = app.GuidancePattern(point1path).x(point1index);
            y1 = app.GuidancePattern(point1path).y(point1index);
            yaw1 = app.GuidancePattern(point1path).yaw(point1index);
            v1 = app.GuidancePattern(point1path).speed(point1index);
            k1 = 0;

            path = headland(x0, y0, yaw0, v0, k0, x1, y1, yaw1, v1, k1, app.defaults.turningType, app.constants);

            app.newGuidancePattern.x = path.x;
            app.newGuidancePattern.y = path.y;
            app.newGuidancePattern.yaw = path.theta;
            app.newGuidancePattern.yaw(app.newGuidancePattern.yaw > pi) = app.newGuidancePattern.yaw(app.newGuidancePattern.yaw > pi) - 2*pi;
            app.newGuidancePattern.yaw(app.newGuidancePattern.yaw < -pi) = app.newGuidancePattern.yaw(app.newGuidancePattern.yaw < -pi) + 2*pi;
            app.newGuidancePattern.speed = path.v;
            app.newGuidancePattern.implement = ones(1,numel(path.x)) * 2;

            % TODO: What if connection is made to the end of path 1 ??

            if(point0index < numel(app.GuidancePattern(point0path).x))
                message(app,"Remove " + num2str(numel(app.GuidancePattern(point0path).x) - point0index) + " points from the end of path " + num2str(point0path));

                app.GuidancePattern(point0path).x = app.GuidancePattern(point0path).x(1:point0index);
                app.GuidancePattern(point0path).y = app.GuidancePattern(point0path).y(1:point0index);
                app.GuidancePattern(point0path).yaw = app.GuidancePattern(point0path).yaw(1:point0index);
                app.GuidancePattern(point0path).speed = app.GuidancePattern(point0path).speed(1:point0index);
                app.GuidancePattern(point0path).implement = app.GuidancePattern(point0path).implement(1:point0index);
            end

            if(point1index > 1)
                message(app,"Remove " + num2str(point1index-1) + " points from the start of path " + num2str(point1path));

                app.GuidancePattern(point1path).x = app.GuidancePattern(point1path).x(point1index:end);
                app.GuidancePattern(point1path).y = app.GuidancePattern(point1path).y(point1index:end);
                app.GuidancePattern(point1path).yaw = app.GuidancePattern(point1path).yaw(point1index:end);
                app.GuidancePattern(point1path).speed = app.GuidancePattern(point1path).speed(point1index:end);
                app.GuidancePattern(point1path).implement = app.GuidancePattern(point1path).implement(point1index:end);
            end


            if(point1path > point0path+1)
                app.GuidancePattern = [app.GuidancePattern(1:point0path) app.newGuidancePattern app.GuidancePattern(point1path:end) app.GuidancePattern(point0path+1:point1path-1)];
            else
                app.GuidancePattern = [app.GuidancePattern(1:point0path) app.newGuidancePattern app.GuidancePattern(point1path:end)];
            end
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            message(app,"");
            message(app,"ISOBUS Guidance Path editor v0.1");
            message(app,"Copyright (c) 2024 Luke / Juha Backman");

            %% initialize empty structs with given fields
            app.GuidancePattern = struct;
            app.GuidancePattern.x = [];
            app.GuidancePattern.y = [];
            app.GuidancePattern.yaw = [];
            app.GuidancePattern.speed = [];
            app.GuidancePattern.implement = [];
            app.GuidancePattern(1) = [];

            app.Obstacle = struct;
            app.Obstacle.x = [];
            app.Obstacle.y = [];
            app.Obstacle.P = [];
            app.Obstacle(1) = [];

            app.PartfieldBoundary = struct;
            app.PartfieldBoundary.x = [];
            app.PartfieldBoundary.y = [];
            app.PartfieldBoundary(1) = [];

            resetSelection(app);


            %% Initialize vehicle parameteres
            app.constants = struct;

            app.constants.dt = 0.1;                                         % calculation resolution

            %% Tractor
            % R_min = 6;                                                      % minimum turning radius
            % Turn_time = 3;                                                  % time from side to side
            % 
            % app.constants.K_max = 1/R_min;                                  % max left curvature
            % app.constants.K_min = -1/R_min;                                 % max right curvature
            % 
            % app.constants.dK_max = (2*app.constants.K_max)/Turn_time;       % from right to left
            % app.constants.dK_min = (2*app.constants.K_min)/Turn_time;       % from left to right
            % 
            % app.constants.dv_max = 0.5;                                     % max acceleration
            % app.constants.dv_min = -0.5;                                    % max deceleration
            % 
            % 
            % app.constants.v_start_sequence = [1.5];                         % turning speed profile
            % app.constants.v_end_sequence = [1.5];
            % 
            % app.constants.v_center_reverse = -1.0;                          % reverse speed

            %% AKI
            R_min = 2 + 0.1;                                                % minimum turning radius + margin 0.1
            Turn_time = 10;                                                 % time from side to side (measured 6)

            app.constants.K_max = 1/R_min;                                  % max left curvature
            app.constants.K_min = -1/R_min;                                 % max right curvature

            app.constants.dK_max = (2*app.constants.K_max)/Turn_time;       % from right to left
            app.constants.dK_min = (2*app.constants.K_min)/Turn_time;       % from left to right

            app.constants.dv_max = 1;                                       % max acceleration (measured 0.25)
            app.constants.dv_min = -1;                                      % max deceleration


            app.constants.v_start_sequence = [0.5];                         % turning speed profile
            app.constants.v_end_sequence = [0.5];

            app.constants.v_center_reverse = -0.5;                          % reverse speed

            %% Initialize default path parameters
            app.defaults = struct;

            app.defaults.turningType = 'FASTEST';

            app.defaults.direction = '1';
            app.defaults.sideshift = '0.5';
            app.defaults.speed = '0.5';
            app.defaults.workstate = '0';

            app.defaults.minPointDistance = '0.5';

        end

        % Menu selected function: LoadMenu
        function LoadButtonPushed(app, event)

            message(app,"");

            [file,path] = uigetfile('*.XML');
            if isequal(file,0)
                message(app,"User selected Cancel");
            else
                message(app,"User selected " + fullfile(path, file));
                task = readstruct(fullfile(path, file),'AttributeSuffix','_att');

                if(~isfield(task,'PFD'))
                    message(app,"File is not proper partfield XML file");
                    return;
                end

                [indx,tf] = listdlg('ListString',{'Obstacles','Guidance Pattern','Field boundary'});
                if tf == true
                    
                    % Define origin for map projection
                    if(numel(app.Origin) == 0)

                        if( isfield(task.PFD,'GGP') && ...
                                isfield(task.PFD.GGP,'GPN') && ...
                                isfield(task.PFD.GGP.GPN(1),'LSG') && ...
                                isfield(task.PFD.GGP.GPN(1).LSG,'PNT'))

                            app.Origin = [task.PFD.GGP.GPN(1).LSG.PNT(1).C_att, task.PFD.GGP.GPN(1).LSG.PNT(1).D_att, 0];
                        elseif( isfield(task.PFD,'PNT') )

                            app.Origin = [task.PFD.PNT(1).C_att, task.PFD.PNT(1).D_att, 0];
                        elseif( isfield(task.PFD,'PLN') && ...
                                isfield(task.PFD.PLN(1),'LSG') && ...
                                isfield(task.PFD.PLN(1).LSG(1),'PNT') )

                            app.Origin = [task.PFD.PLN(1).LSG(1).PNT(1).C_att, task.PFD.PLN(1).LSG(1).PNT(1).D_att, 0];
                        end

                        message(app,"Set new origin: ")
                        message(app,app.Origin)
                    end

                    % read Obstacles
                    if( any(indx == 1) )
                        app.Obstacle = struct;
                        app.Obstacle.x = [];
                        app.Obstacle.y = [];
                        app.Obstacle.P = [];
                        app.Obstacle(1) = [];

                        if(isfield(task.PFD,'PNT'))
                            [obstacle_x, obstacle_y] = latlon2local([task.PFD.PNT(:).C_att],[task.PFD.PNT(:).D_att], app.Origin(3),app.Origin);

                            for i=1:numel(obstacle_x)
                                app.Obstacle(i).x = obstacle_x(i);
                                app.Obstacle(i).y = obstacle_y(i);
                                app.Obstacle(i).P = [task.PFD.PNT(i).H_att, task.PFD.PNT(i).I_att];
                            end

                        end

                        message(app,"Load " + num2str(numel(app.Obstacle)) + " Obstacles")
                    end

                    % read GuidancePattern
                    if( any(indx == 2) )
                        app.GuidancePattern = struct;
                        app.GuidancePattern.x = [];
                        app.GuidancePattern.y = [];
                        app.GuidancePattern.yaw = [];
                        app.GuidancePattern.speed = [];
                        app.GuidancePattern.implement = [];
                        app.GuidancePattern(1) = [];

                        
                        if( isfield(task.PFD,'GGP') && ...
                            isfield(task.PFD.GGP(1),'GPN') )

                            for i = 1:numel(task.PFD.GGP(1).GPN)

                                [x, y] = latlon2local([task.PFD.GGP(1).GPN(i).LSG.PNT(:).C_att],[task.PFD.GGP(1).GPN(i).LSG.PNT(:).D_att], app.Origin(3),app.Origin);

                                app.GuidancePattern(i).x = x;
                                app.GuidancePattern(i).y = y;

                                app.GuidancePattern(i).yaw = NaN(1,numel(app.GuidancePattern(i).x));
                                app.GuidancePattern(i).speed = NaN(1,numel(app.GuidancePattern(i).x));
                                app.GuidancePattern(i).implement = NaN(1,numel(app.GuidancePattern(i).x));

                                for j = 1:numel(app.GuidancePattern(i).x)
                                    C = sscanf(task.PFD.GGP(1).GPN(i).LSG.PNT(j).B_att,"yaw[%f] speed[%f] implement[%f]");

                                    app.GuidancePattern(i).yaw(j) = C(1);
                                    app.GuidancePattern(i).speed(j) = C(2);
                                    app.GuidancePattern(i).implement(j) = C(3);
                                end
                            end
                        end

                        message(app,"Load " + num2str(numel(app.GuidancePattern)) + " Guidance Patterns");
                    end

                    % read Field boundary
                    if( any(indx == 3) )
                        app.PartfieldBoundary = struct;
                        app.PartfieldBoundary.x = [];
                        app.PartfieldBoundary.y = [];
                        app.PartfieldBoundary(1) = [];

                        if( isfield(task.PFD,'PLN') )                      % Element PLN: Polygon
                            for i = 1:numel(task.PFD.PLN)
                                if( isfield(task.PFD.PLN(i),'A_att') && ...    % Attribute A: PolygonType: 
                                    task.PFD.PLN(i).A_att == 1 )               %   1 = Partfield Boyndary    

                                    if( isfield(task.PFD.PLN(i),'LSG'))
                                        for j = 1:numel(task.PFD.PLN(i).LSG)        % TODO: tarkista onko näitä vain yksi standardissa

                                            [x, y] = latlon2local([task.PFD.PLN(i).LSG(j).PNT(:).C_att],[task.PFD.PLN(i).LSG(j).PNT(:).D_att], app.Origin(3),app.Origin);

                                            app.PartfieldBoundary(end+1).x = x;
                                            app.PartfieldBoundary(end).y = y;
                                        end
                                    end
                                end
                            end
                        end

                        message(app,"Load " + num2str(numel(app.PartfieldBoundary)) + " partfield boundaries");
                    end

                    resetSelection(app);
                    UpdateFigures(app,true);
                end
            end

        end

        % Menu selected function: SaveMenu
        function SaveButtonPushed(app, event)

            message(app,"");

            [file,path] = uiputfile('*.XML');
            if isequal(file,0)
                message(app,"User selected Cancel")
            else
                message(app,"User selected " + fullfile(path, file))

                answer = inputdlg({"Partfield name"},"Give task information");

                taskdata.PFD = struct;
                taskdata.PFD.A_att = "PFD0";
                taskdata.PFD.C_att = answer{1};
                taskdata.PFD.D_att = 0;

                % Write Obstacles
                [lat,lon] = local2latlon([app.Obstacle(:).x],[app.Obstacle(:).y],0,app.Origin);
                for i=1:numel(lat)
                    taskdata.PFD.PNT(i).A_att = 5;
                    taskdata.PFD.PNT(i).C_att = sprintf("%.14f",lat(i));
                    taskdata.PFD.PNT(i).D_att = sprintf("%.14f",lon(i));
                    taskdata.PFD.PNT(i).H_att = app.Obstacle(i).P(1);
                    taskdata.PFD.PNT(i).I_att = app.Obstacle(i).P(2);
                end
                message(app,"Converted " + num2str(numel(app.Obstacle)) + " Obstacles to XML")

                % Write Guidance Patterns
                taskdata.PFD.GGP = struct;
                taskdata.PFD.GGP.A_att = "GGP0";
                taskdata.PFD.GGP.GPN(1) = struct;

                for i = 1:numel(app.GuidancePattern)
                    taskdata.PFD.GGP.GPN(i).A_att = "GPN" + (i-1);
                    taskdata.PFD.GGP.GPN(i).B_att = "state[1]";
                    taskdata.PFD.GGP.GPN(i).C_att = 3;
                    taskdata.PFD.GGP.GPN(i).LSG = struct;
                    taskdata.PFD.GGP.GPN(i).LSG.A_att = 5;
                    taskdata.PFD.GGP.GPN(i).LSG.PNT(1) = struct;

                    [lat,lon] = local2latlon(app.GuidancePattern(i).x,app.GuidancePattern(i).y,0,app.Origin);
                    for j = 1:numel(lat)
                        taskdata.PFD.GGP.GPN(i).LSG.PNT(j).A_att = 9;
                        taskdata.PFD.GGP.GPN(i).LSG.PNT(j).B_att = sprintf("yaw[%f] speed[%f] implement[%i]", app.GuidancePattern(i).yaw(j), app.GuidancePattern(i).speed(j), app.GuidancePattern(i).implement(j));
                        taskdata.PFD.GGP.GPN(i).LSG.PNT(j).C_att = sprintf("%.14f",lat(j));
                        taskdata.PFD.GGP.GPN(i).LSG.PNT(j).D_att = sprintf("%.14f",lon(j));
                    end
                    taskdata.PFD.GGP.GPN(i).LSG.PNT(1).A_att = 6;
                    taskdata.PFD.GGP.GPN(i).LSG.PNT(numel(app.GuidancePattern(i).x)).A_att = 7;
                end
                message(app,"Converted " + num2str(numel(app.GuidancePattern)) + " Guidance Patterns to XML")

                writestruct(taskdata, fullfile(path, file),'AttributeSuffix','_att');
                message(app,"XML written")
            end

        end

        % Button down function: UIAxes
        function UIAxesButtonDown(app, event)

            if(strcmp(app.mode,'Insert') && event.Button == 1)
                if(app.GuidancePatternButton.Value)
                    app.newGuidancePattern.x = [app.newGuidancePattern.x event.IntersectionPoint(1)];
                    app.newGuidancePattern.y = [app.newGuidancePattern.y event.IntersectionPoint(2)];
                    app.newGuidancePattern.yaw = [app.newGuidancePattern.yaw 0];
                    app.newGuidancePattern.speed = [app.newGuidancePattern.speed 0];
                    app.newGuidancePattern.implement = [app.newGuidancePattern.implement 0];
                end
                if(app.ObstacleButton.Value)
                    i = numel(app.newObstacle) + 1;
                    app.newObstacle(i).x = event.IntersectionPoint(1);
                    app.newObstacle(i).y = event.IntersectionPoint(2);
                    app.newObstacle(i).P = [0, 0];
                end
            end

            if(strcmp(app.mode,'Copy') || strcmp(app.mode,'Delete') || strcmp(app.mode,'Extend') || strcmp(app.mode,'Connect') || strcmp(app.mode,'Property') || strcmp(app.mode,'Modify'))
                minDistObstacle = Inf;
                minIdxObstacle = 0;
                obstacleSelected = false;
                minDistGuidancePattern = Inf;
                minIdxGuidancePattern = 0;
                guidanceSelected = false;
                minDistPartfieldBoundary = Inf;
                minIdxPartfieldBoundary = 0;
                boundarySelected = false;

                if(numel(app.Obstacle) > 0)
                    coordinates = [app.Obstacle(:).x; app.Obstacle(:).y];
                    coordinates = [coordinates(1,:) - event.IntersectionPoint(1); coordinates(2,:) - event.IntersectionPoint(2)];
                    dist = coordinates(1,:) .* coordinates(1,:) + coordinates(2,:) .* coordinates(2,:);
                    [minDistObstacle, minIdxObstacle] = min(dist);

                    obstacleSelected = any(app.selectIdxObstacle) == true;
                end

                if(numel(app.GuidancePattern) > 0)
                    coordinates = [app.GuidancePattern(:).x; app.GuidancePattern(:).y];
                    coordinates = [coordinates(1,:) - event.IntersectionPoint(1); coordinates(2,:) - event.IntersectionPoint(2)];
                    dist = coordinates(1,:) .* coordinates(1,:) + coordinates(2,:) .* coordinates(2,:);
                    [minDistGuidancePattern, minIdxGuidancePattern] = min(dist);

                    guidanceSelected = any([app.selectGuidancePattern(:).Idx]) == true;
                end

                if(numel(app.PartfieldBoundary) > 0)
                    coordinates = [app.PartfieldBoundary(:).x; app.PartfieldBoundary(:).y];
                    coordinates = [coordinates(1,:) - event.IntersectionPoint(1); coordinates(2,:) - event.IntersectionPoint(2)];
                    dist = coordinates(1,:) .* coordinates(1,:) + coordinates(2,:) .* coordinates(2,:);
                    [minDistPartfieldBoundary, minIdxPartfieldBoundary] = min(dist);

                    boundarySelected = any([app.selectPartfieldBoundary(:).Idx]) == true;
                end

                if(strcmp(app.mode,'Connect') || strcmp(app.mode,'Modify'))
                    guidanceSelected = true;
                end

                if(~obstacleSelected && ~guidanceSelected && ~boundarySelected)
                    [~, select] = min([minDistObstacle minDistGuidancePattern minDistPartfieldBoundary]);
                    switch select
                        case 1
                            obstacleSelected = true;
                        case 2
                            guidanceSelected = true;
                        case 3
                            boundarySelected = true;
                    end
                end

                if(obstacleSelected)

                    % left mouse button
                    if(event.Button == 1)
                        app.selectIdxObstacle(minIdxObstacle) = app.selectIdxObstacle(minIdxObstacle) == false;
                    end

                    if(strcmp(app.mode,'Copy') && event.Button == 1)
                        if(app.GuidancePatternButton.Value)
                            app.newGuidancePattern.x = [app.newGuidancePattern.x app.Obstacle(minIdxObstacle).x];
                            app.newGuidancePattern.y = [app.newGuidancePattern.y app.Obstacle(minIdxObstacle).y];
                            app.newGuidancePattern.yaw = [app.newGuidancePattern.yaw 0];
                            app.newGuidancePattern.speed = [app.newGuidancePattern.speed 0];
                            app.newGuidancePattern.implement = [app.newGuidancePattern.implement 0];

                        end
                        if(app.ObstacleButton.Value)
                            i = numel(app.newObstacle) + 1;
                            app.newObstacle(i).x = app.Obstacle(minIdxObstacle).x;
                            app.newObstacle(i).y = app.Obstacle(minIdxObstacle).y;
                            app.newObstacle(i).P = [0, 0];
                        end
                    end


                    if(strcmp(app.mode,'Copy') && event.Button == 3)
                        P(:,1) = [app.Obstacle(:).x];
                        P(:,2) = [app.Obstacle(:).y];
                        P(:,3) = 0;

                        if(app.GuidancePatternButton.Value && numel(app.newGuidancePattern.y) > 0)
                            [dist, pos] = point_to_line(P, [app.newGuidancePattern.x(end) app.newGuidancePattern.y(end) 0], [app.Obstacle(minIdxObstacle).x app.Obstacle(minIdxObstacle).y 0]);
                            accept = abs(dist) < 0.5 & pos  > -1;

                            [~, idx] = sort(pos);

                            for i=1:numel(idx)
                                if(accept(idx(i)))
                                    if( ((app.newGuidancePattern.x(end) - app.Obstacle(idx(i)).x)^2 + (app.newGuidancePattern.y(end) - app.Obstacle(idx(i)).y)^2) <= 4*4 )
                                        app.newGuidancePattern.x = [app.newGuidancePattern.x app.Obstacle(idx(i)).x];
                                        app.newGuidancePattern.y = [app.newGuidancePattern.y app.Obstacle(idx(i)).y];
                                        app.newGuidancePattern.yaw = [app.newGuidancePattern.yaw 0];
                                        app.newGuidancePattern.speed = [app.newGuidancePattern.speed 0];
                                        app.newGuidancePattern.implement = [app.newGuidancePattern.implement 0];
                                    else
                                        accept(:) = false;
                                    end
                                end
                            end
                        end
                        if(app.ObstacleButton.Value && numel(app.newObstacle) > 0)

                            [dist, pos] = point_to_line(P, [app.newObstacle(end).x app.newObstacle(end).y 0], [app.Obstacle(minIdxObstacle).x app.Obstacle(minIdxObstacle).y 0]);
                            accept = abs(dist) < 0.5 & pos > -1;

                            [~, idx] = sort(pos);

                            for i=1:numel(app.Obstacle)
                                if(accept(idx(i)))
                                    if( ((app.newObstacle(end).x - app.Obstacle(idx(i)).x)^2 + (app.newObstacle(end).y - app.Obstacle(idx(i)).y)^2) <= 3*3 )
                                        j = numel(app.newObstacle) + 1;
                                        app.newObstacle(j).x = app.Obstacle(idx(i)).x;
                                        app.newObstacle(j).y = app.Obstacle(idx(i)).y;
                                        app.newObstacle(j).P = [0, 0];
                                    else
                                        accept(:) = false;
                                    end
                                end
                            end

                        end
                    end
                end

                if(guidanceSelected)

                    i = 1;
                    while minIdxGuidancePattern > numel(app.selectGuidancePattern(i).Idx)
                        minIdxGuidancePattern = minIdxGuidancePattern - numel(app.selectGuidancePattern(i).Idx);
                        i = i + 1;
                    end

                    % left mouse button
                    if(event.Button == 1)
                        app.selectGuidancePattern(i).Idx(minIdxGuidancePattern) =  app.selectGuidancePattern(i).Idx(minIdxGuidancePattern) == false;
                        % right mouse button
                    elseif(event.Button == 3)
                        app.selectGuidancePattern(i).Idx(:) =  app.selectGuidancePattern(i).Idx(minIdxGuidancePattern) == false;
                    end

                    if(strcmp(app.mode,'Copy') && event.Button == 1)
                        if(app.GuidancePatternButton.Value)
                            app.newGuidancePattern.x = [app.newGuidancePattern.x app.GuidancePattern(i).x(minIdxGuidancePattern)];
                            app.newGuidancePattern.y = [app.newGuidancePattern.y app.GuidancePattern(i).y(minIdxGuidancePattern)];
                            app.newGuidancePattern.yaw = [app.newGuidancePattern.yaw app.GuidancePattern(i).yaw(minIdxGuidancePattern)];
                            app.newGuidancePattern.speed = [app.newGuidancePattern.speed app.GuidancePattern(i).speed(minIdxGuidancePattern)];
                            app.newGuidancePattern.implement = [app.newGuidancePattern.implement app.GuidancePattern(i).implement(minIdxGuidancePattern)];
                        end
                        if(app.ObstacleButton.Value)
                            k = numel(app.newObstacle) + 1;
                            app.newObstacle(k).x = app.GuidancePattern(i).x(minIdxGuidancePattern);
                            app.newObstacle(k).y = app.GuidancePattern(i).y(minIdxGuidancePattern);
                            app.newObstacle(k).P = [0, 0];
                        end
                    end

                    if(strcmp(app.mode,'Copy') && event.Button == 3)
                        if(app.GuidancePatternButton.Value)
                            app.newGuidancePattern.x = [app.newGuidancePattern.x app.GuidancePattern(i).x(1,:)];
                            app.newGuidancePattern.y = [app.newGuidancePattern.y app.GuidancePattern(i).y(1,:)];
                            app.newGuidancePattern.yaw = [app.newGuidancePattern.yaw app.GuidancePattern(i).yaw(1,:)];
                            app.newGuidancePattern.speed = [app.newGuidancePattern.speed app.GuidancePattern(i).speed(1,:)];
                            app.newGuidancePattern.implement = [app.newGuidancePattern.implement app.GuidancePattern(i).implement(1,:)];
                        end
                        if(app.ObstacleButton.Value)
                            k = numel(app.newObstacle);
                            for(j=1:numel(app.GuidancePattern(i).x))
                                app.newObstacle(k+j).x = app.GuidancePattern(i).x(j);
                                app.newObstacle(k+j).y = app.GuidancePattern(i).y(j);
                                app.newObstacle(k+j).P = [0, 0];
                            end
                        end
                    end
                end


                if(boundarySelected)

                    i = 1;
                    while minIdxPartfieldBoundary > numel(app.selectPartfieldBoundary(i).Idx)
                        minIdxPartfieldBoundary = minIdxPartfieldBoundary - numel(app.selectPartfieldBoundary(i).Idx);
                        i = i + 1;
                    end

                    % left mouse button
                    if(event.Button == 1)
                        app.selectPartfieldBoundary(i).Idx(minIdxPartfieldBoundary) = app.selectPartfieldBoundary(i).Idx(minIdxPartfieldBoundary) == false;
                    end

                    if(strcmp(app.mode,'Copy') && event.Button == 1)
                        if(app.GuidancePatternButton.Value)
                            app.newGuidancePattern.x = [app.newGuidancePattern.x app.PartfieldBoundary(i).x(minIdxPartfieldBoundary)];
                            app.newGuidancePattern.y = [app.newGuidancePattern.y app.PartfieldBoundary(i).y(minIdxPartfieldBoundary)];
                            app.newGuidancePattern.yaw = [app.newGuidancePattern.yaw 0];
                            app.newGuidancePattern.speed = [app.newGuidancePattern.speed 0];
                            app.newGuidancePattern.implement = [app.newGuidancePattern.implement 0];

                        end
                        if(app.ObstacleButton.Value)
                            i = numel(app.newObstacle) + 1;
                            app.newObstacle(i).x = app.PartfieldBoundary(i).x(minIdxPartfieldBoundary);
                            app.newObstacle(i).y = app.PartfieldBoundary(i).y(minIdxPartfieldBoundary);
                            app.newObstacle(i).P = [0, 0];
                        end
                    end


                    if(strcmp(app.mode,'Copy') && event.Button == 3)
                        P(:,1) = app.PartfieldBoundary(i).x;
                        P(:,2) = app.PartfieldBoundary(i).y;
                        P(:,3) = 0;

                        if(app.GuidancePatternButton.Value && numel(app.newGuidancePattern.y) > 0)
                            [dist, pos] = point_to_line(P, [app.newGuidancePattern.x(end) app.newGuidancePattern.y(end) 0], [app.PartfieldBoundary(i).x(minIdxPartfieldBoundary) app.PartfieldBoundary(i).y(minIdxPartfieldBoundary) 0]);
                            accept = abs(dist) < 3 & pos > -1;

                            [~, idx] = sort(pos);
                            last_idx = idx(find(accept(idx(:)),1));

                            for j=1:numel(idx)
                                if(accept(idx(j)))
                                    if( abs(last_idx-idx(j)) <= 1 ) 
                                        last_idx = idx(j);
                                        app.newGuidancePattern.x = [app.newGuidancePattern.x app.PartfieldBoundary(i).x(idx(j))];
                                        app.newGuidancePattern.y = [app.newGuidancePattern.y app.PartfieldBoundary(i).y(idx(j))];
                                        app.newGuidancePattern.yaw = [app.newGuidancePattern.yaw 0];
                                        app.newGuidancePattern.speed = [app.newGuidancePattern.speed 0];
                                        app.newGuidancePattern.implement = [app.newGuidancePattern.implement 0];

                                        app.selectPartfieldBoundary(i).Idx(idx(j)) = true;
                                    else
                                        accept(:) = false;
                                    end
                                end
                            end
                        end
                        if(app.ObstacleButton.Value && numel(app.newObstacle) > 0)

                            [dist, pos] = point_to_line(P, [app.newObstacle(end).x app.newObstacle(end).y 0],  [app.PartfieldBoundary(i).x(minIdxPartfieldBoundary) app.PartfieldBoundary(i).y(minIdxPartfieldBoundary) 0]);
                            accept = abs(dist) < 3 & pos  > -1;

                            [~, idx] = sort(pos);

                            for j=1:numel(app.Obstacle)
                                if(accept(idx(j)))
                                    if( abs(last_idx-idx(j)) <= 1 ) 
                                        last_idx = idx(j);

                                        k = numel(app.newObstacle) + 1;
                                        app.newObstacle(k).x = app.PartfieldBoundary(i).x(idx(j));
                                        app.newObstacle(k).y = app.PartfieldBoundary(i).y(idx(j));
                                        app.newObstacle(k).P = [0, 0];

                                        app.selectPartfieldBoundary(i).Idx(idx(j)) = true;
                                    else
                                        accept(:) = false;
                                    end
                                end
                            end

                        end
                    end
                end
            end
            

            UpdateFigures(app);


            function [dist, pos] = point_to_line(P, B, A)
                A = repmat(A,size(P,1),1);
                B = repmat(B,size(P,1),1);
                AB = A - B;
                AP = P - A;

                normsq = sum(AB.^2,2);
                dist = sum(cross(AB,AP,2),2) ./ sqrt(normsq);
                pos = dot(AB,AP,2) ./ normsq;
            end
        end

        % Button pushed function: InsertButton
        function InsertButtonPushed(app, event)

            message(app,"");

            if(app.GuidancePatternButton.Value)
                message(app,"Click on the map for new path points");
                message(app,"Create new guidance pattern by pressing ""Ready"" Button");
                message(app,"Cancel by selecting some other function");
            end
            if(app.ObstacleButton.Value)
                message(app,"Click on the map for inserting new obstacle");
                message(app,"Create new obstacles by pressing ""Ready"" Button");
                message(app,"Cancel by selecting some other function");
            end

            setMode(app,'Insert');
            resetSelection(app);
            UpdateFigures(app);

        end

        % Button pushed function: CopyButton
        function CopyButtonPushed(app, event)

            message(app,"");

            if(app.GuidancePatternButton.Value)
                message(app,"Select points to be copied: left button for single position, right button for whole path");
                message(app,"Create new guidance pattern by pressing ""Ready"" Button");
                message(app,"Cancel by selecting some other function");
            end
            if(app.ObstacleButton.Value)
                message(app,"Select points to be copied: left button for single position, right button for whole path");
                message(app,"Create new obstacles by pressing ""Ready"" Button");
                message(app,"Cancel by selecting some other function");
            end

            setMode(app,'Copy');
            resetSelection(app);
            UpdateFigures(app);

        end

        % Button pushed function: DeleteButton
        function DeleteButtonPushed(app, event)
            message(app,"");
            message(app,"Select points to be deleted: left button for single position, right button for whole path");
            message(app,"Remove points by pressing ""Ready"" Button");
            message(app,"Cancel by selecting some other function");

            setMode(app,'Delete');
            resetSelection(app);
            UpdateFigures(app);
        end

        % Button pushed function: ExtendButton
        function ExtendButtonPushed(app, event)
            message(app,"");
            message(app,"Select path end/start point to extend");
            message(app,"Remove points by pressing ""Ready"" Button");
            message(app,"Cancel by selecting some other function");

            setMode(app,'Extend');
            resetSelection(app);
            UpdateFigures(app);
        end

        % Button pushed function: ConnectButton
        function ConnectButtonPushed(app, event)

            message(app,"");
            message(app,"Select start and end point of connection path");
            message(app,"Generate path by pressing ""Ready"" Button");
            message(app,"Cancel by selecting some other function");
            setMode(app,'Connect');
            resetSelection(app);
            UpdateFigures(app);

        end

        % Button pushed function: PointpropertyButton
        function PointpropertyButtonPushed(app, event)

            message(app,"");
            message(app,"Select points to be edited");
            message(app,"Edit point properties by pressing ""Ready"" Button");
            message(app,"Cancel by selecting some other function");
            setMode(app,'Property');
            resetSelection(app);
            UpdateFigures(app);

        end

        % Button pushed function: ModifypathButton
        function ModifypathButtonPushed(app, event)

            message(app,"");
            message(app,"Select paths to be modified");
            message(app,"Modify paths by pressing ""Ready"" Button");
            message(app,"Cancel by selecting some other function");
            setMode(app,'Modify');
            resetSelection(app);
            UpdateFigures(app);
        end

        % Button pushed function: READYButton
        function READYButtonPushed(app, event)

            message(app,"");

            if(strcmp(app.mode,'Insert'))
                if(numel(app.newObstacle) > 0)
                    startIdx = numel(app.Obstacle)+1;
                    app.Obstacle(startIdx:(startIdx+numel(app.newObstacle)-1)) = app.newObstacle(:);

                    message(app,"Created " + numel(app.newObstacle) + " new obstacles");
                end
                if(numel(app.newGuidancePattern.x) > 0)
                    startIdx = numel(app.GuidancePattern)+1;
                    app.GuidancePattern(startIdx) = app.newGuidancePattern;

                    answer = inputdlg({"Speed","Work state"},"Give path information",1,{app.defaults.speed,app.defaults.workstate});

                    if(numel(answer) > 0)
                        app.defaults.speed = answer{1};
                        app.defaults.workstate = answer{2};

                        app.selectGuidancePattern(startIdx).Idx = true(1,numel(app.newGuidancePattern.x));
                        modifyPathPoints(app,0,str2double(answer{1}),str2double(answer{2}));

                        message(app,"Created new guidance pattern with " + numel(app.newGuidancePattern.x) + " new points");

                        if(startIdx > 1)
                            selection = uiconfirm(app.ISOBUSGuidancePatheditorUIFigure, "Connect to previous path?", "Confirm auto connect","Icon","question");

                            if(strcmp(selection,"OK"))
                                connectPaths(app, startIdx-1, numel(app.GuidancePattern(startIdx-1).x), startIdx, 1);
                                message(app,"Connected to previous path");
                            end
                        end
                    end
                end
            end

            if(strcmp(app.mode,'Delete'))
                app.Obstacle = app.Obstacle(app.selectIdxObstacle == false);

                GuidancePatternNotEmpty = true(1,numel(app.GuidancePattern));

                for i = 1:numel(app.GuidancePattern)
                    app.GuidancePattern(i).x = app.GuidancePattern(i).x(app.selectGuidancePattern(i).Idx == false);
                    app.GuidancePattern(i).y = app.GuidancePattern(i).y(app.selectGuidancePattern(i).Idx == false);
                    app.GuidancePattern(i).yaw = app.GuidancePattern(i).yaw(app.selectGuidancePattern(i).Idx == false);
                    app.GuidancePattern(i).speed = app.GuidancePattern(i).speed(app.selectGuidancePattern(i).Idx == false);
                    app.GuidancePattern(i).implement = app.GuidancePattern(i).implement(app.selectGuidancePattern(i).Idx == false);

                    GuidancePatternNotEmpty(i) = numel(app.GuidancePattern(i).x) > 0;
                end
                app.GuidancePattern = app.GuidancePattern(GuidancePatternNotEmpty);


                if(any(app.selectIdxObstacle))
                    message(app,"Removed " + sum(app.selectIdxObstacle) + " obstacles");
                end

                if(any([app.selectGuidancePattern(:).Idx]))
                    message(app,"Removed " + sum([app.selectGuidancePattern(:).Idx]) + " guidance pattern points");
                end

                if(any(GuidancePatternNotEmpty == false))
                    message(app,"Removed " + sum(GuidancePatternNotEmpty == false) + " guidance patterns");
                end
            end

            if(strcmp(app.mode,'Copy'))
                if(numel(app.newObstacle) > 0)
                    message(app,"Give parameters to be modify");
                    message(app,"Leave input field empty to not modify");
                    answer = inputdlg({"Move x","Move y"},"Modified parameter");

                    if(numel(answer) > 0)
                        if(~isnan(str2double(answer{1})))
                            for(i=1:numel(app.newObstacle))
                                app.newObstacle(i).x = app.newObstacle(i).x + str2double(answer{1});
                            end
                        end
                        if(~isnan(str2double(answer{2})))
                            for(i=1:numel(app.newObstacle))
                                app.newObstacle(i).y = app.newObstacle(i).y + str2double(answer{2});
                            end
                        end

                        startIdx = numel(app.Obstacle)+1;
                        app.Obstacle(startIdx:(startIdx+numel(app.newObstacle)-1)) = app.newObstacle(:);

                        message(app,"Copied " + numel(app.newObstacle) + " new obstacles");
                    else
                        message(app,"User cancelled");
                    end
                end
                if(numel(app.newGuidancePattern.x) > 0)
                    message(app,"Give parameters to be modify");
                    message(app,"Leave input field empty to not modify");
                    answer = inputdlg({"Direction", "Sideshift","Speed","Work state"},"Give path information",1,{app.defaults.direction, app.defaults.sideshift,app.defaults.speed,app.defaults.workstate});

                    if(numel(answer) > 0)
                        app.defaults.direction = answer{1};
                        app.defaults.sideshift = answer{2};
                        app.defaults.speed = answer{3};
                        app.defaults.workstate = answer{4};

                        if(strcmp(app.defaults.direction,"-1") || strcmp(app.defaults.direction,"BACKWARDS"))
                             app.newGuidancePattern.x = app.newGuidancePattern.x(end:-1:1);
                             app.newGuidancePattern.y = app.newGuidancePattern.y(end:-1:1);
                             app.newGuidancePattern.yaw = app.newGuidancePattern.yaw(end:-1:1);
                             app.newGuidancePattern.speed = app.newGuidancePattern.speed(end:-1:1);
                             app.newGuidancePattern.implement = app.newGuidancePattern.implement(end:-1:1);
                        end

                        startIdx = numel(app.GuidancePattern)+1;
                        app.GuidancePattern(startIdx) = app.newGuidancePattern;

                        resetSelection(app);
                        app.selectGuidancePattern(startIdx).Idx = true(1,numel(app.GuidancePattern(startIdx).x));

                        modifyPathPoints(app,str2double(answer{2}),str2double(answer{3}),str2double(answer{4}));

                        message(app,"Copied new guidance pattern with " + numel(app.GuidancePattern(startIdx).x) + " new points");

                        if(startIdx > 1)

                            selection = uiconfirm(app.ISOBUSGuidancePatheditorUIFigure, "Connect to previous path?", "Confirm auto connect","Icon","question");

                            if(strcmp(selection,"OK"))
                                connectPaths(app, startIdx-1, numel(app.GuidancePattern(startIdx-1).x), startIdx, 1);
                                message(app,"Connected to previous path");
                            end
                        end
                    else
                        message(app,"User cancelled");
                    end
                end
            end

            if(strcmp(app.mode,'Extend'))
                points_selected = 0;
                for i = 1:numel(app.GuidancePattern)
                    points_selected = points_selected + sum(app.selectGuidancePattern(i).Idx);
                end

                if(points_selected == 1)
                    for i = 1:numel(app.GuidancePattern)
                        if(any(app.selectGuidancePattern(i).Idx))
                            if(app.selectGuidancePattern(i).Idx(1))
                                % first path point selected

                                message(app,"Give extend length");
                                answer = inputdlg({"Extend"},"Parameter");

                                extend = NaN;
                                if(numel(answer) > 0)
                                    extend = str2double(answer{1});
                                end

                                if(~isnan(extend))

                                    dx = app.GuidancePattern(i).x(2) - app.GuidancePattern(i).x(1);
                                    dy = app.GuidancePattern(i).y(2) - app.GuidancePattern(i).y(1);

                                    length = sqrt(dx*dx + dy*dy);

                                    dx = dx / length * extend;
                                    dy = dy / length * extend;

                                    app.GuidancePattern(i).x = [app.GuidancePattern(i).x(1)-dx,  app.GuidancePattern(i).x];
                                    app.GuidancePattern(i).y = [app.GuidancePattern(i).y(1)-dy,  app.GuidancePattern(i).y];
                                    app.GuidancePattern(i).yaw = [app.GuidancePattern(i).yaw(1),  app.GuidancePattern(i).yaw];
                                    app.GuidancePattern(i).speed = [app.GuidancePattern(i).speed(1),  app.GuidancePattern(i).speed];
                                    app.GuidancePattern(i).implement = [app.GuidancePattern(i).implement(1),  app.GuidancePattern(i).implement];
                                else
                                    message(app,"Not valid parameter!"); 
                                end
                            
                            elseif(app.selectGuidancePattern(i).Idx(end))
                                % last path point selected

                                message(app,"Give extend length");
                                answer = inputdlg({"Extend"},"Parameter");

                                extend = NaN;
                                if(numel(answer) > 0)
                                    extend = str2double(answer{1});
                                end

                                if(~isnan(extend))

                                    dx = app.GuidancePattern(i).x(end) - app.GuidancePattern(i).x(end-1);
                                    dy = app.GuidancePattern(i).y(end) - app.GuidancePattern(i).y(end-1);

                                    length = sqrt(dx*dx + dy*dy);

                                    dx = dx / length * extend;
                                    dy = dy / length * extend;

                                    app.GuidancePattern(i).x = [app.GuidancePattern(i).x, app.GuidancePattern(i).x(end)+dx];
                                    app.GuidancePattern(i).y = [app.GuidancePattern(i).y, app.GuidancePattern(i).y(end)+dy];
                                    app.GuidancePattern(i).yaw = [app.GuidancePattern(i).yaw, app.GuidancePattern(i).yaw(end)];
                                    app.GuidancePattern(i).speed = [app.GuidancePattern(i).speed, app.GuidancePattern(i).speed(end)];
                                    app.GuidancePattern(i).implement = [app.GuidancePattern(i).implement, app.GuidancePattern(i).implement(end)];
                                else
                                    message(app,"Not valid parameter!"); 
                                end

                            else
                                message(app,"Only first or last point can be selected!"); 
                            end
                        end
                    end
                else
                   message(app,"Only one point can be selected!"); 
                end
            end

            if(strcmp(app.mode,'Connect'))
                fail = false;

                point0found = false;
                point0path = NaN;
                point0index = NaN;

                point1found = false;
                point1path = NaN;
                point1index = NaN;

                for i = 1:numel(app.GuidancePattern)
                    if(~fail && any(app.selectGuidancePattern(i).Idx))

                        if(sum(app.selectGuidancePattern(i).Idx) == 1)
                            if(~point0found)
                                point0found = true;
                                point0path = i;
                                point0index = find(app.selectGuidancePattern(i).Idx);
                            elseif(~point1found)
                                point1found = true;
                                point1path = i;
                                point1index = find(app.selectGuidancePattern(i).Idx);
                            else
                                message(app,"Only two points or paths can be selected");
                                fail = true;
                            end
                        else
                            first = find(app.selectGuidancePattern(i).Idx,1,'first');
                            last = find(app.selectGuidancePattern(i).Idx,1,'last');

                            if(any(app.selectGuidancePattern(i).Idx(first:last) == false))
                                message(app,"A uniform area must be chosen from the path");
                                fail = true;
                            end

                            if(~point0found)
                                point0found = true;
                                point0path = i;
                                point0index = [first, last];

                            elseif(~point1found)
                                point1found = true;
                                point1path = i;
                                point1index = [first, last];
                            else
                                message(app,"Only two points or paths can be selected");
                                fail = true;
                            end
                        end
                    end
                end

                if(~fail && ~point1found)
                    message(app,"Two points or paths must be selected");
                    fail = true;
                end

                if(~fail)
                    if(numel(point0index) == 2 || numel(point1index) == 2)
                        if(numel(point0index) == 1)
                            point0index = [point0index point0index];
                        elseif(numel(point1index) == 1)
                            point1index = [point1index point1index];
                        end

                        x0 = [app.GuidancePattern(point0path).x(point0index(1)) app.GuidancePattern(point0path).x(point0index(2))];
                        y0 = [app.GuidancePattern(point0path).y(point0index(1)) app.GuidancePattern(point0path).y(point0index(2))];

                        x1 = [app.GuidancePattern(point1path).x(point1index(1)) app.GuidancePattern(point1path).x(point1index(2))];
                        y1 = [app.GuidancePattern(point1path).y(point1index(1)) app.GuidancePattern(point1path).y(point1index(2))];

                        d11 = (x0(1) - x1(1))^2 + (y0(1) - y1(1))^2;
                        d12 = (x0(1) - x1(2))^2 + (y0(1) - y1(2))^2;
                        d21 = (x0(2) - x1(1))^2 + (y0(2) - y1(1))^2;
                        d22 = (x0(2) - x1(2))^2 + (y0(2) - y1(2))^2;

                        [~, idx] = min([d11 d12 d21 d22]);
                        if(idx==1 || idx==2)
                            point0index = point0index(1);
                        else
                            point0index = point0index(2);
                        end

                        if(idx==1 || idx==3)
                            point1index = point1index(1);
                        else
                            point1index = point1index(2);
                        end
                    end

                    if(point0index == 1 || point1index == numel(app.GuidancePattern(point1path).x))
                        % points are propably in opposite order
                        tmp = point0path;
                        point0path = point1path;
                        point1path = tmp;

                        tmp = point0index;
                        point0index = point1index;
                        point1index = tmp;
                    end

                    message(app,"Connect paths id1: " + num2str(point0path) + " pos1: " + num2str(point0index) + " id2: " + num2str(point1path) + " pos2: " + num2str(point1index));

                    connectPaths(app, point0path, point0index, point1path, point1index);
                end

                message(app,"Connected ");
            end

            if(strcmp(app.mode,'Property'))
                if(any(app.selectIdxObstacle))

                    message(app,"Give parameters to be modify");
                    message(app,"Leave input field empty to not modify");
                    answer = inputdlg({"Move x","Move y"},"Modified parameter");

                    if(numel(answer) > 0 && ~isnan(str2double(answer{1})))
                        for i=1:numel(app.selectIdxObstacle)
                            if(app.selectIdxObstacle(i))
                                app.Obstacle(i).x = app.Obstacle(i).x + str2double(answer{1});
                            end
                        end
                    end
                    if(numel(answer) > 0 && ~isnan(str2double(answer{2})))
                        for i=1:numel(app.selectIdxObstacle)
                            if(app.selectIdxObstacle(i))
                                app.Obstacle(i).y = app.Obstacle(i).y + str2double(answer{2});
                            end
                        end
                    end

                    if(numel(answer) > 0)
                        message(app,"Modified " + sum(app.selectIdxObstacle) + " obstacles");
                    else
                        message(app,"User cancelled");
                    end

                end

                if(any([app.selectGuidancePattern(:).Idx]))

                    message(app,"Give parameters to be modify");
                    message(app,"Leave input field empty to not modify");
                    answer = inputdlg({"Sideshift","Speed","Work state"},"Give path information");

                    if(numel(answer) > 0)
                        modifyPathPoints(app,str2double(answer{1}),str2double(answer{2}),str2double(answer{3}));
                        message(app,"Modified " + sum([app.selectGuidancePattern(:).Idx]) + " guidance pattern points");
                    else
                        message(app,"User cancelled");
                    end
                end
            end

            if(strcmp(app.mode,'Modify'))
                selectedPaths = [];
                PathNames = {};

                for(i=1:numel(app.selectGuidancePattern))
                    if(any(app.selectGuidancePattern(i).Idx))
                        selectedPaths(end+1) = i;
                    end
                    PathNames{i} = "Path" + num2str(i);
                end

                [selectedPaths,tf] = listdlg('ListString',PathNames, 'InitialValue', selectedPaths );

                if(tf)
                    resetSelection(app);

                    for(i=selectedPaths)
                        app.selectGuidancePattern(i).Idx(:) = true;
                    end

                    UpdateFigures(app);

                    [selectedOperations,tf] = listdlg('ListString',{'Merge paths', 'Insert zero to direction change point', 'Cut paths when speed is zero', 'Normalize point distances'});

                    if(tf)
                        % Merge paths
                        if( any(selectedOperations == 1) )
                            if(any(diff(selectedPaths) > 1))
                                message(app,"Only sequential paths can be merged");
                                message(app,"All operations cancelled");
                                selectedOperations = [];
                            else
                                app.newGuidancePattern.x = [app.GuidancePattern(selectedPaths).x];
                                app.newGuidancePattern.y = [app.GuidancePattern(selectedPaths).y];
                                app.newGuidancePattern.yaw = [app.GuidancePattern(selectedPaths).yaw];
                                app.newGuidancePattern.speed = [app.GuidancePattern(selectedPaths).speed];
                                app.newGuidancePattern.implement = [app.GuidancePattern(selectedPaths).implement];

                                if(selectedPaths(1) > 1 && selectedPaths(end) < numel(app.GuidancePattern))
                                    app.GuidancePattern = [app.GuidancePattern(1:(selectedPaths(1)-1)) app.newGuidancePattern app.GuidancePattern((selectedPaths(end)+1):end)];
                                elseif(selectedPaths(1) > 1)
                                    app.GuidancePattern = [app.GuidancePattern(1:(selectedPaths(1)-1)) app.newGuidancePattern];
                                elseif(selectedPaths(end) < numel(app.GuidancePattern))
                                    app.GuidancePattern = [app.newGuidancePattern app.GuidancePattern((selectedPaths(end)+1):end)];
                                else
                                    app.GuidancePattern = [app.newGuidancePattern];
                                end

                                message(app,"Merged " + num2str(numel(selectedPaths)) + " paths");
                                selectedPaths = selectedPaths(1);
                            end
                        end

                        % Insert zero to direction change point
                        if( any(selectedOperations == 2) )

                            for i = 1:numel(app.GuidancePattern)
                                if any(selectedPaths == i)
                                    startPos = 1;

                                    while startPos < numel(app.GuidancePattern(i).x)

                                        if(app.GuidancePattern(i).speed(startPos) > 0)
                                            changePos = find(app.GuidancePattern(i).speed(startPos+1:end) < 0, 1) + startPos;
    
                                            if(isempty(changePos))
                                                break;
                                            end
                                        else
                                            changePos = find(app.GuidancePattern(i).speed(startPos+1:end) > 0, 1) + startPos;
    
                                            if(isempty(changePos))
                                                break;
                                            end
                                        end
    
                                        changePos = changePos-1;
                                        app.GuidancePattern(i).x = [app.GuidancePattern(i).x(1:changePos) app.GuidancePattern(i).x(changePos:end)];
                                        app.GuidancePattern(i).y = [app.GuidancePattern(i).y(1:changePos) app.GuidancePattern(i).y(changePos:end)];
                                        app.GuidancePattern(i).yaw = [app.GuidancePattern(i).yaw(1:changePos) app.GuidancePattern(i).yaw(changePos:end)];
                                        app.GuidancePattern(i).speed = [app.GuidancePattern(i).speed(1:(changePos-1)), 0, app.GuidancePattern(i).speed(changePos+1), app.GuidancePattern(i).speed((changePos+1):end)];
                                        app.GuidancePattern(i).implement = [app.GuidancePattern(i).implement(1:changePos) app.GuidancePattern(i).implement(changePos:end)];

                                        startPos = changePos+1;
                                    end
                                end
                            end
                        end

                        % Cut paths when speed is zero
                        if( any(selectedOperations == 3) )

                            app.newGuidancePattern(1) = [];
                            newSelectedPaths = [];

                            for i = 1:numel(app.GuidancePattern)
                                if any(selectedPaths == i)
                                    endPos = 0;
                                    while endPos < numel(app.GuidancePattern(i).x)
                                        startPos = find(app.GuidancePattern(i).speed(endPos+1:end) ~= 0, 1) + endPos;
                                        if isempty(startPos)
                                            break;
                                        end

                                        endPos = find(app.GuidancePattern(i).speed(startPos+1:end) == 0, 1) + startPos;
                                        if isempty(endPos)
                                            endPos = numel(app.GuidancePattern(i).x);
                                        end

                                        app.newGuidancePattern(end+1).x = app.GuidancePattern(i).x(startPos:endPos);
                                        app.newGuidancePattern(end).y = app.GuidancePattern(i).y(startPos:endPos);
                                        app.newGuidancePattern(end).yaw = app.GuidancePattern(i).yaw(startPos:endPos);
                                        app.newGuidancePattern(end).speed = app.GuidancePattern(i).speed(startPos:endPos);
                                        app.newGuidancePattern(end).implement = app.GuidancePattern(i).implement(startPos:endPos);

                                        newSelectedPaths(end+1) = numel(app.newGuidancePattern);
                                    end
                                else
                                    app.newGuidancePattern(end+1) = app.GuidancePattern(i);
                                end
                            end

                            message(app,"Cut " + num2str(numel(selectedPaths)) + " paths. Created " + num2str(numel(app.newGuidancePattern)-numel(app.GuidancePattern)) + " new paths.");
                            app.GuidancePattern =  app.newGuidancePattern;
                            selectedPaths = newSelectedPaths;
                        end

                        % Normalize point distances
                        if( any(selectedOperations == 4) )
                            answer = inputdlg({"Minimum distance between path points"},"Give path information",1,{app.defaults.minPointDistance});
                            if(numel(answer) > 0)
                                minDistance = str2double(answer{1});
                                if isnan(minDistance)
                                    minDistance = str2double(app.defaults.minPointDistance);
                                    message(app,"Given distance is not proper. Use default: " + app.defaults.minPointDistance);
                                else
                                    app.defaults.minPointDistance = answer{1};
                                end

                                for i = selectedPaths
                                    dx = app.GuidancePattern(i).x(2:end) - app.GuidancePattern(i).x(1:(end-1));
                                    dy = app.GuidancePattern(i).y(2:end) - app.GuidancePattern(i).y(1:(end-1));
                                    length = sqrt(dx.*dx + dy.*dy);

                                    resetSelection(app);
                                    app.newGuidancePattern.x(1) = app.GuidancePattern(i).x(1);
                                    app.newGuidancePattern.y(1) = app.GuidancePattern(i).y(1);
                                    app.newGuidancePattern.yaw(1) = app.GuidancePattern(i).yaw(1);
                                    app.newGuidancePattern.speed(1) = app.GuidancePattern(i).speed(1);
                                    app.newGuidancePattern.implement(1) = app.GuidancePattern(i).implement(1);

                                    distance = 0;
                                    for j=1:numel(length)
                                        distance = distance + length(j);
                                        if(distance >= minDistance || (distance > 0 && app.newGuidancePattern.implement(end) ~= app.GuidancePattern(i).implement(1+j)))
                                            if(distance > minDistance * 2)
                                                steps = floor(distance / minDistance);

                                                dx = (app.GuidancePattern(i).x(1+j) - app.newGuidancePattern.x(end)) / steps;
                                                dy = (app.GuidancePattern(i).y(1+j) - app.newGuidancePattern.y(end)) / steps;
                                                dspeed = (app.GuidancePattern(i).speed(1+j) - app.newGuidancePattern.speed(end)) / steps;

                                                for k=1:(steps-1)
                                                    app.newGuidancePattern.x(end+1) = app.newGuidancePattern.x(end) + dx;
                                                    app.newGuidancePattern.y(end+1) = app.newGuidancePattern.y(end) + dy;
                                                    app.newGuidancePattern.yaw(end+1) = app.GuidancePattern(i).yaw(1+j);
                                                    app.newGuidancePattern.speed(end+1) = app.newGuidancePattern.speed(end) + dspeed;
                                                    app.newGuidancePattern.implement(end+1) = app.GuidancePattern(i).implement(1+j);
                                                end
                                            end

                                            app.newGuidancePattern.x(end+1) = app.GuidancePattern(i).x(1+j);
                                            app.newGuidancePattern.y(end+1) = app.GuidancePattern(i).y(1+j);
                                            app.newGuidancePattern.yaw(end+1) = app.GuidancePattern(i).yaw(1+j);
                                            app.newGuidancePattern.speed(end+1) = app.GuidancePattern(i).speed(1+j);
                                            app.newGuidancePattern.implement(end+1) = app.GuidancePattern(i).implement(1+j);
                                            distance = 0;
                                        end
                                    end
                                    if(distance > 0)
                                        app.newGuidancePattern.x(end+1) = app.GuidancePattern(i).x(end);
                                        app.newGuidancePattern.y(end+1) = app.GuidancePattern(i).y(end);
                                        app.newGuidancePattern.yaw(end+1) = app.GuidancePattern(i).yaw(end);
                                        app.newGuidancePattern.speed(end+1) = app.GuidancePattern(i).speed(end);
                                        app.newGuidancePattern.implement(end+1) = app.GuidancePattern(i).implement(end);

                                        message(app,"Path no. " + num2str(i) + " last distance: " + distance);
                                    else
                                        message(app,"Path no. " + num2str(i) + " all distances greater than " + minDistance);
                                    end

                                    app.GuidancePattern(i) = app.newGuidancePattern;
                                end
                            else
                                message(app,"User cancelled");
                            end
                        end
                    else
                        message(app,"User cancelled");
                    end
                else
                    message(app,"User cancelled");
                end
            end


            setMode(app,'');
            resetSelection(app);
            UpdateFigures(app);
        end

        % Window key press function: ISOBUSGuidancePatheditorUIFigure
        function ISOBUSGuidancePatheditorUIFigureWindowKeyPress(app, event)
            key = event.Key;

            if( key == 'return')
                READYButtonPushed(app);
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ISOBUSGuidancePatheditorUIFigure and hide until all components are created
            app.ISOBUSGuidancePatheditorUIFigure = uifigure('Visible', 'off');
            app.ISOBUSGuidancePatheditorUIFigure.Position = [100 100 920 702];
            app.ISOBUSGuidancePatheditorUIFigure.Name = 'ISOBUS Guidance Path editor';
            app.ISOBUSGuidancePatheditorUIFigure.WindowKeyPressFcn = createCallbackFcn(app, @ISOBUSGuidancePatheditorUIFigureWindowKeyPress, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.ISOBUSGuidancePatheditorUIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadMenu
            app.LoadMenu = uimenu(app.FileMenu);
            app.LoadMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadMenu.Text = 'Load';

            % Create SaveMenu
            app.SaveMenu = uimenu(app.FileMenu);
            app.SaveMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveMenu.Text = 'Save';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.ISOBUSGuidancePatheditorUIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {'1x', 100};

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GridLayout);
            app.GridLayout3.ColumnWidth = {'1x', 200};
            app.GridLayout3.RowHeight = {'1x'};
            app.GridLayout3.Layout.Row = 2;
            app.GridLayout3.Layout.Column = 1;

            % Create READYButton
            app.READYButton = uibutton(app.GridLayout3, 'push');
            app.READYButton.ButtonPushedFcn = createCallbackFcn(app, @READYButtonPushed, true);
            app.READYButton.Enable = 'off';
            app.READYButton.Layout.Row = 1;
            app.READYButton.Layout.Column = 2;
            app.READYButton.Text = 'READY';

            % Create MessageLabel
            app.MessageLabel = uilabel(app.GridLayout3);
            app.MessageLabel.BackgroundColor = [0.902 0.902 0.902];
            app.MessageLabel.VerticalAlignment = 'top';
            app.MessageLabel.Layout.Row = 1;
            app.MessageLabel.Layout.Column = 1;
            app.MessageLabel.Text = 'Message';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GridLayout);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 1;

            % Create MAPTab
            app.MAPTab = uitab(app.TabGroup);
            app.MAPTab.Title = 'MAP';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.MAPTab);
            app.GridLayout4.ColumnWidth = {'1x', 200};
            app.GridLayout4.RowHeight = {'1x'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout4);
            title(app.UIAxes, 'MAP')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Layout.Row = 1;
            app.UIAxes.Layout.Column = 1;
            app.UIAxes.ButtonDownFcn = createCallbackFcn(app, @UIAxesButtonDown, true);

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.GridLayout4);
            app.GridLayout2.ColumnWidth = {'1x'};
            app.GridLayout2.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout2.Layout.Row = 1;
            app.GridLayout2.Layout.Column = 2;

            % Create DeleteButton
            app.DeleteButton = uibutton(app.GridLayout2, 'push');
            app.DeleteButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteButtonPushed, true);
            app.DeleteButton.Layout.Row = 5;
            app.DeleteButton.Layout.Column = 1;
            app.DeleteButton.Text = 'Delete';

            % Create ConnectButton
            app.ConnectButton = uibutton(app.GridLayout2, 'push');
            app.ConnectButton.ButtonPushedFcn = createCallbackFcn(app, @ConnectButtonPushed, true);
            app.ConnectButton.Layout.Row = 7;
            app.ConnectButton.Layout.Column = 1;
            app.ConnectButton.Text = 'Connect';

            % Create EDITButtonGroup
            app.EDITButtonGroup = uibuttongroup(app.GridLayout2);
            app.EDITButtonGroup.Title = 'EDIT';
            app.EDITButtonGroup.Layout.Row = [1 2];
            app.EDITButtonGroup.Layout.Column = 1;

            % Create GuidancePatternButton
            app.GuidancePatternButton = uiradiobutton(app.EDITButtonGroup);
            app.GuidancePatternButton.Text = 'Guidance Pattern';
            app.GuidancePatternButton.Position = [11 59 115 22];
            app.GuidancePatternButton.Value = true;

            % Create ObstacleButton
            app.ObstacleButton = uiradiobutton(app.EDITButtonGroup);
            app.ObstacleButton.Text = 'Obstacle';
            app.ObstacleButton.Position = [11 37 69 22];

            % Create PointpropertyButton
            app.PointpropertyButton = uibutton(app.GridLayout2, 'push');
            app.PointpropertyButton.ButtonPushedFcn = createCallbackFcn(app, @PointpropertyButtonPushed, true);
            app.PointpropertyButton.Layout.Row = 8;
            app.PointpropertyButton.Layout.Column = 1;
            app.PointpropertyButton.Text = 'Point property';

            % Create ModifypathButton
            app.ModifypathButton = uibutton(app.GridLayout2, 'push');
            app.ModifypathButton.ButtonPushedFcn = createCallbackFcn(app, @ModifypathButtonPushed, true);
            app.ModifypathButton.Layout.Row = 9;
            app.ModifypathButton.Layout.Column = 1;
            app.ModifypathButton.Text = 'Modify path';

            % Create InsertButton
            app.InsertButton = uibutton(app.GridLayout2, 'push');
            app.InsertButton.ButtonPushedFcn = createCallbackFcn(app, @InsertButtonPushed, true);
            app.InsertButton.Layout.Row = 3;
            app.InsertButton.Layout.Column = 1;
            app.InsertButton.Text = 'Insert';

            % Create CopyButton
            app.CopyButton = uibutton(app.GridLayout2, 'push');
            app.CopyButton.ButtonPushedFcn = createCallbackFcn(app, @CopyButtonPushed, true);
            app.CopyButton.Layout.Row = 4;
            app.CopyButton.Layout.Column = 1;
            app.CopyButton.Text = 'Copy';

            % Create ExtendButton
            app.ExtendButton = uibutton(app.GridLayout2, 'push');
            app.ExtendButton.ButtonPushedFcn = createCallbackFcn(app, @ExtendButtonPushed, true);
            app.ExtendButton.Layout.Row = 6;
            app.ExtendButton.Layout.Column = 1;
            app.ExtendButton.Text = 'Extend';

            % Create GraphTab
            app.GraphTab = uitab(app.TabGroup);
            app.GraphTab.Title = 'Graph';

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.GraphTab);
            app.GridLayout5.ColumnWidth = {'1x'};
            app.GridLayout5.RowHeight = {'1x', '1x', '1x'};

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.GridLayout5);
            title(app.UIAxes2, 'Speed')
            xlabel(app.UIAxes2, 'time [s]')
            ylabel(app.UIAxes2, 'speed [m/s]')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Layout.Row = 1;
            app.UIAxes2.Layout.Column = 1;

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.GridLayout5);
            title(app.UIAxes3, 'Heading')
            xlabel(app.UIAxes3, 'time [s]')
            ylabel(app.UIAxes3, 'yaw [rad]')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.Layout.Row = 2;
            app.UIAxes3.Layout.Column = 1;

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.GridLayout5);
            title(app.UIAxes4, 'Work state')
            xlabel(app.UIAxes4, 'time [s]')
            ylabel(app.UIAxes4, 'state')
            zlabel(app.UIAxes4, 'Z')
            app.UIAxes4.Layout.Row = 3;
            app.UIAxes4.Layout.Column = 1;

            % Show the figure after all components are created
            app.ISOBUSGuidancePatheditorUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = route_editor_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ISOBUSGuidancePatheditorUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ISOBUSGuidancePatheditorUIFigure)
        end
    end
end