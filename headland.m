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
%  License along with ISOBUS Guidance Path editor.
%  If not, see <http://www.gnu.org/licenses/>.
%



% Generate the fastest path between given positions
function path = headland(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, turn_type, constants)

% check special path types
if(~iscell(turn_type))
    if(strcmp(turn_type, 'FASTEST'))
        turn_type = {'LRL', 'RLR', 'LSL', 'RSR', 'LSR', 'RSL', 'L|R|L', 'R|L|R'};
    elseif(strcmp(turn_type, 'FASTEST_FORWARDS'))
        turn_type = {'LRL', 'RLR', 'LSL', 'RSR', 'LSR', 'RSL'};
    elseif(strcmp(turn_type, 'FASTEST_BACKWARDS'))
        turn_type = {'L|R|L', 'R|L|R'};
    else
        turn_type = {turn_type};
    end
end

paths = cell(0);
for i=1:numel(turn_type)
    if(strcmp(turn_type{i}, 'LRL'))
        paths{i} = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, constants.K_max, constants.K_min, constants.K_max, constants.dK_min, NaN, constants.dK_min,0);
    end
    if(strcmp(turn_type{i}, 'RLR'))
        paths{i} = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, constants.K_min, constants.K_max, constants.K_min, constants.dK_max, NaN, constants.dK_max,0);
    end
    if(strcmp(turn_type{i}, 'LSL'))
        paths{i} = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, constants.K_max, 0, constants.K_max, constants.dK_min, constants.dK_max, constants.dK_min,0);
    end
    if(strcmp(turn_type{i}, 'RSR'))
        paths{i} = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, constants.K_min, 0, constants.K_min, constants.dK_max, constants.dK_min, constants.dK_max,0);
    end
    if(strcmp(turn_type{i}, 'LSR'))
        paths{i} = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, constants.K_max, 0, constants.K_min, constants.dK_min, 0, constants.dK_max,0);
    end
    if(strcmp(turn_type{i}, 'RSL'))
        paths{i} = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, constants.K_min, 0, constants.K_max, constants.dK_max, 0, constants.dK_min,0);
    end
    if(strcmp(turn_type{i}, 'L|R|L'))
        paths{i} = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, constants.K_max, constants.K_min, constants.K_max, constants.dK_min, constants.dK_max, constants.dK_min, constants.v_center_reverse);
    end
    if(strcmp(turn_type{i}, 'R|L|R'))
        paths{i} = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, constants.K_min, constants.K_max, constants.K_min, constants.dK_max, constants.dK_min, constants.dK_max, constants.v_center_reverse);
    end
end

% Find the fastest
disp(' ');
disp('--------------------------------------------------------------------');
disp(' ');

best_t = Inf;
best_idx = NaN;

for i=1:numel(paths)
    if(isfield(paths{i},'t'))
        disp("Path idx: " + num2str(i) + " name: " + turn_type{i} + " t: " + num2str(paths{i}.t(end)));
        if(paths{i}.t(end) < best_t)
            best_idx = i;
            best_t = paths{i}.t(end);
        end
    end
end

if(~isnan(best_idx))
    disp("Fastest path name:  name: " + turn_type{best_idx} + " t: " + num2str(paths{best_idx}.t(end)));

    path = paths{best_idx};
    path.type = turn_type{best_idx};
else
    path = [];
end

disp(' ');
disp('--------------------------------------------------------------------');
disp(' ');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the path between given positions using different turning
% sequences. Used sequence is determined by the parameter sequence.
function path = path_generic_turn(x0, y0, theta0, v0, k0, x1, y1, theta1, v1, k1, constants, k_start, k_center, k_end, dK_start, dK_center, dK_end, v_center)

%% TODO: k_start == 0 or k_end == 0

% Initial null path
path = NaN;

phi_threshold = pi/2;

vC11_idx = [];
vC32_idx = 1;

vC21_idx = [];
vC22_idx = 1;

% Paths:
% S1 -> C1 -> S2 -> C2 -> S3 -> C3 -> S4

fprintf('PATH (dx: %g, dy: %g, dtheta: %g\n',x1-x0,y1-y0,theta1-theta0);
iteration = 1;

recalculate = true;
while(recalculate)

    fprintf('Iteration %d: %g --> %g --> %g --> %g --> %g\n',iteration,k0,k_start,k_center,k_end,k1);
    iteration = iteration+1;

    if(iteration > 200)
        fprintf('Maximum iteration\n');
        return;
    end

    %%
    % Create first connecting spiral (S1)
    %   Start --> S1 --> Left1
    %   k0 --> k_start
    %
    v = generate_speed(v0,[],'start',constants);           % v = v0 --> sequence
    K = generate_curvature(k0,k_start,constants);          % K = k0 --> k_start

    % make sure that the trajectories has equal length
    [v vC11] = reshape_trajectory(v,numel(K),'start');

    % simulation doesn't actually use the last speed and arch first point is exactly
    % the same than spiral last:
    vC11 = [v(end) vC11];

    % fprintf('S1\n');
    % fprintf('numel(K): %g\n', numel(K));
    % fprintf('numel(v): %g\n', numel(v));
    % fprintf('numel(vC11): %g\n', numel(vC11));

    pathS1 = simulate_trajectory(x0,y0,theta0,v,K,constants);

    %plot(pathS1.Ox(end), pathS1.Oy(end),'rx');

    %%
    % Create last connecting spiral (S4)
    %   Left2 --> S4 --> End
    %   k_end --> k1
    %
    v = generate_speed([],v1,'end',constants);             % v = sequence --> v1
    K = generate_curvature(k_end,k1,constants);            % K = k_start --> k1

    % make sure that the trajectories has equal length
    [v vC32] = reshape_trajectory(v,numel(K),'end');

    % fprintf('S4\n');
    % fprintf('numel(K): %g\n', numel(K));
    % fprintf('numel(v): %g\n', numel(v));
    % fprintf('numel(vC32): %g\n', numel(vC32));

    pathS4 = simulate_trajectory(0,0,0,v,K,constants);
    pathS4 = relocate_path(x1,y1,theta1, pathS4, numel(pathS4.t));

    %plot(pathS4.Ox(1), pathS4.Oy(1),'rx');

    %%
    % Create second connecting spiral (S2)
    %   Left1 --> S2 --> Right
    %   k_start --> k_center
    %
    if(isempty(vC11_idx))
        % initial guess: Arch C1 long enought that vC1 is used completely
        vC11_idx = numel(vC11);
    end

    vC12 = [];
    if(v_center < 0)
        %% Center is driven with reverse speed

        K_1 = generate_curvature(k_start,0,constants);          % K = k_start --> 0
        v_1 = generate_speed(vC11(end),0,'',constants);         % v = vC11(end) --> 0

        v_1 = [vC11(vC11_idx:end) v_1];

        %         [v_1 vC12] = reshape_trajectory(v_1,numel(K_1),'end');
        %
        %         if(numel(vC12) > 1 && vC11_idx ~= numel(vC11))
        %             fprintf('Path too short to change speed\n');
        %             return;
        %         end

        if(vC11_idx ~= numel(vC11))
            K_1 = [K_1 zeros(1,numel(v_1)-numel(K_1))];
            vC12 = v_1(1);
        else
            [v_1 vC12] = reshape_trajectory(v_1,numel(K_1),'end');
        end

        K_2 = generate_curvature(0,k_center,constants);         % K = 0 ---> k_center
        v_2 = generate_speed(0,v_center,'',constants);          % v = vC1(end) --> 0

        [v_2 vC21] = reshape_trajectory(v_2,numel(K_2),'start');
        vC21 = [v_2(end) vC21];

        K = [K_1 K_2];
        v = [v_1 v_2];

    else
        K = generate_curvature(k_start,k_center,constants);    % K = k_start --> k_center

        [v vC21] = reshape_trajectory(vC11(vC11_idx:end),numel(K),'start');
        vC21 = [v(end) vC21];
    end

    % fprintf('S2\n');
    % fprintf('numel(K): %g\n', numel(K));
    % fprintf('numel(v): %g\n', numel(v));
    % fprintf('numel(vC21): %g\n', numel(vC21));
    % fprintf('numel(vC12): %g\n', numel(vC12));

    pathS2 = simulate_trajectory(0,0,0,v,K,constants);

    %%
    % Create third connecting spiral (S3)
    %   Right --> S3 --> Left2
    %   k_center --> k_end
    %

    vC31 = [];
    if(v_center < 0)
        %% Center is driven with reverse speed

        K_1 = generate_curvature(k_center,0,constants);        % K = k_center --> 0
        v_1 = generate_speed(v_center,0,'',constants);         % v = vC1(end) --> 0

        [v_1 vC22] = reshape_trajectory(v_1,numel(K_1),'end');

        K_2 = generate_curvature(0,k_end,constants);            % K = 0 ---> k_end
        v_2 = generate_speed(0,vC32(1),'',constants);           % v = vC1(end) --> 0

        v_2 = [v_2 vC32(1:vC32_idx)];

        %         [v_2 vC31] = reshape_trajectory(v_2,numel(K_2),'start');
        %         vC31 = [v_2(end) vC31];
        %
        %         if(numel(vC31) > 1 && vC32_idx ~= 1)
        %             fprintf('path is too short for speed changes\n');
        %             return;
        %         end

        if(vC32_idx ~= 1)
            K_2 = [zeros(1,numel(v_2)-numel(K_2)) K_2];
            vC31 = v_2(end);
        else
            [v_2 vC31] = reshape_trajectory(v_2,numel(K_2),'start');
            vC31 = [v_2(end) vC31];
        end

        K = [K_1 K_2];
        v = [v_1 v_2];

    else
        K = generate_curvature(k_center,k_end,constants);          % K = k_center --> k_end

        % FIXME: WHY THIS IS NOT WORKING? (last value of spiral is not used
        % in simulation
        %[v vC22] = reshape_trajectory([vC32(1:vC32_idx) vC32(vC32_idx)],numel(K),'end');
        [v vC22] = reshape_trajectory(vC32(1:vC32_idx),numel(K),'end');

    end

    % fprintf('S3\n');
    % fprintf('numel(K): %g\n', numel(K));
    % fprintf('numel(v): %g\n', numel(v));
    % fprintf('numel(vC22): %g\n', numel(vC22));
    % fprintf('numel(vC31): %g\n', numel(vC31));

    pathS3 = simulate_trajectory(0,0,0,v,K,constants);

    %%
    if k_center == 0 % ~isnan(dK_center)

        %%
        % Orientation of centre line
        %

        d = norm([pathS1.Ox(end) - pathS4.Ox(1) pathS1.Oy(end) - pathS4.Oy(1)]);

        %      r1:n ja r2:n merkin mukaan ottamalla päästään yksinkertaisemmilla!
        %
        %         Q1 = [pathS2.x(end) pathS2.y(end)];
        %         Q2 = Q1 + [cos(pathS2.theta(end)) sin(pathS2.theta(end))];
        %         P = [pathS2.Ox(1) pathS2.Oy(1)];
        %         r1 = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);  %distance from P to line that goes through Q1 and Q2
        %
        %         Q1 = [pathS3.x(1) pathS3.y(1)];
        %         Q2 = Q1 + [cos(pathS3.theta(1)) sin(pathS3.theta(1))];
        %         P = [pathS3.Ox(end) pathS3.Oy(end)];
        %         r2 = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);
        %
        %         phi_hat = atan2(pathS4.Oy(1) - pathS1.Oy(end), pathS4.Ox(1) - pathS1.Ox(end));
        %
        %         if sign(dK_start) == sign(dK_end)
        %             % LSL or RSR
        %             if( abs(r1 - r2) >  d)
        %                 fprintf('Centre line cannot be created\n');
        %                 return;
        %             end
        %             if(sign(dK_start) < 0)
        %                 % LSL
        %                 phi = phi_hat - asin((r2-r1)/d);
        %             else
        %                 % RSR
        %                 phi = phi_hat + asin((r2-r1)/d);
        %             end
        %         else
        %             % LSR or RSL
        %             if abs(r1 + r2) > d
        %                 fprintf('Centre line cannot be created\n');
        %                 return;
        %             end
        %
        %             if(sign(dK_start) < 0)
        %                 % LSR
        %                 phi = phi_hat + asin((r2+r1)/d);
        %             else
        %                 % RSL
        %                 phi = phi_hat - asin((r2+r1)/d);
        %             end
        %         end

        r1 = cos(pathS2.theta(end))*(pathS2.Oy(1)-pathS2.y(end)) - sin(pathS2.theta(end))*(pathS2.Ox(1)-pathS2.x(end));
        r2 = cos(pathS3.theta(1))*(pathS3.Oy(end)-pathS3.y(1)) - sin(pathS3.theta(1))*(pathS3.Ox(end)-pathS3.x(1));

        phi_hat = atan2(pathS4.Oy(1) - pathS1.Oy(end), pathS4.Ox(1) - pathS1.Ox(end));

        if( abs(r1 - r2) >  d)
            fprintf('Centre line cannot be created\n');
            return;
        end

        phi = phi_hat - asin((r2-r1)/d);

        if isnan(phi)
            %TODO: Quick fix! Find out real reason why phi == NaN when
            %k_start and/or k_end go to zero

            error('phi is nan');
            return;
        end

        if phi < -pi
            phi = phi + 2*pi;
        elseif phi > pi
            phi = phi - 2*pi;
        end

        %%
        % Relocate S2 and S3
        %

        pathS2 = relocate_path3(pathS2, pathS1.Ox(end), pathS1.Oy(end), phi, 0);
        pathS3 = relocate_path3(pathS3, pathS4.Ox(1), pathS4.Oy(1), phi, 1);


        %%
        % Check feasibility
        %


        if ( abs(atan2(pathS3.y(1) - pathS2.y(end), pathS3.x(1) - pathS2.x(end)) - phi) >  phi_threshold )
            k_center = k_center + dK_center*constants.dt;

            if ((k_center - k_start)*sign(dK_start) < 0 || (k_center - k_end)*sign(dK_end) < 0 || dK_center == 0)
                fprintf('Centre arch cannot be reduced more\n');
                return;
            end

            continue;
        end

    else

        %%
        % Position of centre arch
        %

        if(abs(pathS2.K(1)) < eps)
            if(pathS3.K(end) == 0)
                fprintf('first and last arch reduced to zero\n');
                return;
            end

            d3 = sqrt((pathS3.Ox(end)-pathS3.Ox(1))^2 + (pathS3.Oy(end)-pathS3.Oy(1))^2);
            r1 = cos(pathS2.theta(1))*(pathS2.Oy(end)-pathS2.y(1)) - sin(pathS2.theta(1))*(pathS2.Ox(end)-pathS2.x(1));

            d1 = cos(pathS1.theta(end))*(pathS4.Oy(1)-pathS1.y(end)) - sin(pathS1.theta(end))*(pathS4.Ox(1)-pathS1.x(end));

            l1 = d1 - r1;


            if(d3^2 - l1^2 < 0)
                % Error, path cannot be created

                fprintf('centre arch position is not found\n');
                return;
            end

            l2 = sqrt(d3^2 - l1^2);

            C2_Ox = pathS4.Ox(1) + l1*sin(pathS1.theta(end)) + l2*cos(pathS1.theta(end));
            C2_Oy = pathS4.Oy(1) - l1*cos(pathS1.theta(end)) + l2*sin(pathS1.theta(end));

            pathS2 = relocate_path3(pathS2, C2_Ox, C2_Oy, pathS1.theta(end), 1);
            pathS3 = relocate_path2(pathS3, C2_Ox, C2_Oy, pathS4.Ox(1), pathS4.Oy(1));

        elseif(abs(pathS3.K(end)) < eps)

            d2 = sqrt((pathS2.Ox(end)-pathS2.Ox(1))^2 + (pathS2.Oy(end)-pathS2.Oy(1))^2);
            r2 = cos(pathS3.theta(end))*(pathS3.Oy(1)-pathS3.y(end)) - sin(pathS3.theta(end))*(pathS3.Ox(1)-pathS3.x(end));

            d1 = cos(pathS4.theta(1))*(pathS1.Oy(end)-pathS4.y(1)) - sin(pathS4.theta(1))*(pathS1.Ox(end)-pathS4.x(1));

            l1 = d1 - r2;


            if(d2^2 - l1^2 < 0)
                % Error, path cannot be created

                fprintf('centre arch position is not found\n');
                return;
            end

            l2 = sqrt(d2^2 - l1^2);

            C2_Ox = pathS1.Ox(end) + l1*sin(pathS4.theta(1)) - l2*cos(pathS4.theta(1));
            C2_Oy = pathS1.Oy(end) - l1*cos(pathS4.theta(1)) - l2*sin(pathS4.theta(1));

            pathS2 = relocate_path2(pathS2, pathS1.Ox(end), pathS1.Oy(end), C2_Ox, C2_Oy);
            pathS3 = relocate_path3(pathS3, C2_Ox, C2_Oy, pathS4.theta(1), 0);

        else

            d2 = sqrt((pathS2.Ox(end)-pathS2.Ox(1))^2 + (pathS2.Oy(end)-pathS2.Oy(1))^2);
            d3 = sqrt((pathS3.Ox(end)-pathS3.Ox(1))^2 + (pathS3.Oy(end)-pathS3.Oy(1))^2);

            d14 = sqrt((pathS1.Ox(end) - pathS4.Ox(1))^2 + (pathS1.Oy(end) - pathS4.Oy(1))^2);

            l1 = (d14^2 - d3^2 + d2^2) / (2 * d14);
            D = l1/d14;
            x = pathS1.Ox(end) + D*(pathS4.Ox(1) - pathS1.Ox(end));
            y = pathS1.Oy(end) + D*(pathS4.Oy(1) - pathS1.Oy(end));

            if(d2^2 - l1^2 < 0)
                % Error, path cannot be created

                fprintf('centre arch position is not found\n');
                return;
            end

            l2 = sqrt(d2^2 - l1^2);
            D = l2/l1;
            if k_center < 0
                % Centre turn right
                C2_Ox = x + D*(pathS1.Oy(end) - y);
                C2_Oy = y - D*(pathS1.Ox(end) - x);
            else
                % Centre turn left
                C2_Ox = x - D*(pathS1.Oy(end) - y);
                C2_Oy = y + D*(pathS1.Ox(end) - x);
            end

            %%
            % Relocate S2 and S3
            %


            %plot(pathS1.Ox(end), pathS1.Oy(end),'rx');
            %plot(pathS4.Ox(1), pathS4.Oy(1),'rx');
            %plot(C2_Ox,C2_Oy,'rx');


            pathS2 = relocate_path2(pathS2, pathS1.Ox(end), pathS1.Oy(end), C2_Ox, C2_Oy);
            pathS3 = relocate_path2(pathS3, C2_Ox, C2_Oy, pathS4.Ox(1), pathS4.Oy(1));
        end
    end

    %%
    % Create arch between S1 and S2 (C1)
    %

    [pathC1 vC11_idxn vC12_idxn] = create_arch(pathS1.Ox(end),pathS1.Oy(end),pathS1.x(end),pathS1.y(end),pathS2.x(1),pathS2.y(1),sign(k_start),vC11,vC12,constants);

    %%
    % Create arch/line between S2 and S3 (C2)


    if(isempty(vC21_idx))
        % initial guess: Arch C2 long enought that vC21 is used completely
        vC21_idx = numel(vC21);
    end


    if(k_center == 0)
        [pathC2 vC21_idxn vC22_idxn] = create_straight_line(pathS2.x(end),pathS2.y(end),pathS3.x(1),pathS3.y(1),vC21(1:vC21_idx),vC22(vC22_idx:numel(vC22)),constants);
    else
        [pathC2 vC21_idxn vC22_idxn] = create_arch(pathS2.Ox(end),pathS2.Oy(end),pathS2.x(end),pathS2.y(end),pathS3.x(1),pathS3.y(1),sign(k_center),vC21(1:vC21_idx),vC22(vC22_idx:numel(vC22)),constants);
    end

    %%
    % Create arch between S3 and S4 (C3)
    %

    [pathC3 vC31_idxn vC32_idxn] = create_arch(pathS3.Ox(end),pathS3.Oy(end),pathS3.x(end),pathS3.y(end),pathS4.x(1),pathS4.y(1),sign(k_end),vC31,vC32,constants);

    if(isempty(pathC3))
        % This should be thrown as an exeption, but now special return
        % value indicates this error
        fprintf('Cannot change direction in arch\n');
        return;
    end

    % fprintf('vC11_idxn = %g (%g)\n', vC11_idxn, vC11_idx);
    % fprintf('vC12_idxn = %g\n', vC12_idxn);
    % fprintf('vC21_idxn = %g (%g)\n', vC21_idxn, vC21_idx);
    % fprintf('vC22_idxn = %g (%g)\n', vC22_idxn, vC22_idx);
    % fprintf('vC31_idxn = %g\n', vC31_idxn);
    % fprintf('vC32_idxn = %g (%g)\n', vC32_idxn, vC32_idx);

    %%
    % Check feasibility
    %

    recalculate = false;

    % Arch C1 OK? Is there loop in the path (90% of full circle)
    if(abs(max_turning_angle([pathS1.theta pathC1.theta pathS2.theta])) > 1.8*pi)

        k_start = k_start + dK_start*constants.dt;
        if(k_start + dK_start*constants.dt < k0  && sign(dK_start) < 0 )
            % C1 = Right
            k_start = k0;
        elseif(k_start + dK_start*constants.dt > k0 && sign(dK_start) > 0 )
            % C1 = Left
            k_start = k0;
        end
        recalculate = true;
    end

    % Arch C2 OK? Is there loop in the path (90% of full circle) if reverse
    if(v_center < 0 && ~isnan(dK_center) && abs(max_turning_angle([pathS2.theta pathC2.theta pathS3.theta])) > 1.8*pi)

        k_center = k_center + dK_center*constants.dt;
        if(k_center + dK_center*constants.dt < 0  && sign(dK_center) < 0 )
            % C2 = Right
            k_center = 0;
        elseif(k_center + dK_center*constants.dt > 0 && sign(dK_center) > 0 )
            % C2 = Left
            k_center = 0;
        end
        recalculate = true;
    end

    % Arch C3 OK? Is there loop in the path (90% of full circle)
    if(abs(max_turning_angle([pathS3.theta pathC3.theta pathS4.theta])) > 1.8*pi)
        k_end = k_end + dK_end*constants.dt;
        if(k_end + dK_end*constants.dt < k1 && sign(dK_end) < 0 )
            % C3 = Right
            k_end = k1;
        elseif(k_end + dK_end*constants.dt > k1 && sign(dK_end) > 0)
            % C3 = Left
            k_end = k1;
        end
        recalculate = true;
    end

    % vC1 index OK?
    if(vC11_idx ~= vC11_idxn)
        vC11_idx = vC11_idxn;
        recalculate = true;
    end

    % vC3 index OK?
    if(vC32_idx ~= vC32_idxn)
        vC32_idx = vC32_idxn;
        recalculate = true;
    end


    % vC2 index OK?
    if(recalculate == true)
        % Reset center indexes if some other indexes changes
        vC21_idx = [];
        vC22_idx = 1;

    elseif(vC21_idxn ~= vC21_idx || vC22_idxn ~= 1)

        % Path is not long enough for speed changes
        % Either acceleration in reverse motion or deceleration on forward
        % motion --> reduce max backward or increase min forward speed

        if(isfield(pathC2,'angle'))
            path_length = pathC2.angle*pathC2.r;
        else
            if(numel(pathC2.x > 1))
                path_length = norm([pathC2.x(end)-pathC2.x(1) pathC2.y(end)-pathC2.y(1)]);
            else
                path_length = 0;
            end
        end

        dvC22 = 0;
        dvC21 = 0;

        while(abs(sum(vC21(1:vC21_idx)*constants.dt) + sum(vC22(vC22_idx:numel(vC22))*constants.dt)) > path_length)

            if(vC22_idx < numel(vC22) && dvC21 >= dvC22)
                vC22_idx = vC22_idx + 1;
                dvC22 = dvC22 + abs(vC22(vC22_idx) - vC22(vC22_idx-1));
            elseif(vC21_idx > 1 && dvC21 < dvC22)
                vC21_idx = vC21_idx - 1;
                dvC21 = dvC21 + abs(vC21(vC21_idx) - vC21(vC21_idx+1));
            else
                fprintf('Could not match the speed in center section\n');
                return;
            end
        end

        recalculate = true;
    end

    % Final feasibility for reduced turns at start or end
    if(recalculate == false)

        if(isempty(pathC1.t))

            % difference of theta between spirals cannot be greater than
            % difference of thetas inside spirals
            dtheta_max = 0;
            if(numel(pathS1.theta) >  1)
                dtheta_max = abs(pathS1.theta(end) - pathS1.theta(end-1));
            end
            if(numel(pathS2.theta) >  1)
                dtheta_max = max(dtheta_max, abs(pathS2.theta(2) - pathS2.theta(1)));
            end

            if( abs(pathS1.theta(end) - pathS2.theta(1)) > dtheta_max )
                fprintf('C1 theta mismatch %g vs %g \n',pathS1.theta(end),pathS2.theta(1));
                return;
            end
        end

        if(isempty(pathC3.t))

            dtheta_max = 0;
            if(numel(pathS3.theta) >  1)
                dtheta_max = abs(pathS3.theta(end) - pathS3.theta(end-1));
            else
                dtheta_max = abs(pathC2.theta(end) - pathC2.theta(end-1));
            end
            if(numel(pathS4.theta) >  1)
                dtheta_max = max(dtheta_max, abs(pathS4.theta(2) - pathS4.theta(1)));
            end

            if( abs(pathS3.theta(end) - pathS4.theta(1)) > dtheta_max )
                fprintf('C3 theta mismatch %g vs %g (diff: %g, max %g)\n',pathS3.theta(end),pathS4.theta(1),abs(pathS3.theta(end)-pathS4.theta(1)),dtheta_max);
                return;
            end
        end
    end
end

%%
% Connect paths
%

pathS1.ID = ones(1, numel(pathS1.t)) * 1;
pathC1.ID = ones(1, numel(pathC1.t)) * 2;
pathS2.ID = ones(1, numel(pathS2.t)) * 3;
pathC2.ID = ones(1, numel(pathC2.t)) * 4;
pathS3.ID = ones(1, numel(pathS3.t)) * 5;
pathC3.ID = ones(1, numel(pathC3.t)) * 6;
pathS4.ID = ones(1, numel(pathS4.t)) * 7;

path = connect_paths(pathS1,pathC1);
path = connect_paths(path,pathS2);
path = connect_paths(path,pathC2);
path = connect_paths(path,pathS3);
path = connect_paths(path,pathC3);
path = connect_paths(path,pathS4);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate speed trajectory from v0 to v1
function v = generate_speed(v0,v1,params,constants)
    v = [];
    
    if(strcmp(params,'start'))
        % starting sequence from v0
        
        v = [generate_speed(v0,constants.v_start_sequence(1),[],constants), constants.v_start_sequence];
        
    elseif(strcmp(params,'end'))
        % ending sequence to v1
        
        v = [constants.v_end_sequence, generate_speed(constants.v_end_sequence(end),v1,[],constants)];
        
    elseif(v0 <= v1)
        % maximum acceleration from v0 to v1
        
        v(1) = v0;
        while(v(end) < v1)
            v(end+1) = v(end) + constants.dv_max*constants.dt;
        end
        v(end) = v1;
        
    elseif(v0 > v1)
        % maximum deceleration from v0 to v1
        
        v(1) = v0;
        while(v(end) > v1)
            v(end+1) = v(end) + constants.dv_min*constants.dt;
        end
        v(end) = v1;
        
    end
    
    
    % Time to change direction
    if(v0 == 0)
        v = [0 0 0 v];
    end
    
    if(v1 == 0)
        v = [v 0 0 0];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate curvature trajectory from k0 to k1
function K = generate_curvature(k0,k1,constants)
    K = [];
    
    if(k0 <= k1)
        % maximum right to left
        
        K(1) = k0;
        while(K(end) < k1)
            K(end+1) = K(end) + constants.dK_max*constants.dt;
        end
        K(end) = k1;
        
    elseif(k0 > k1)
        % maximum left to right
        
        K(1) = k0;
        while(K(end) > k1)
            K(end+1) = K(end) + constants.dK_min*constants.dt;
        end
        K(end) = k1;
        
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out1 out2] = reshape_trajectory(in, length, options)

if(strcmp(options,'start'))

    if(length >= numel(in))
        out1 = in;
        out1(end+1:length) = in(end);
        out2 = out1(end);
    else
        out1 = in(1:length);
        out2 = in(length+1:end);
    end
    
elseif(strcmp(options,'end'))
    
    if(length >= numel(in))
        start_idx = length-numel(in);
        out1(1:start_idx) = in(1);
        out1(start_idx+1:length) = in(1:end);
        out2 = in(1);
    else
        out1 = in(numel(in)-length+1:end);
        out2 = in(1:numel(in)-length);
    end
   
else
    error('incorrect parameters');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the trajectory from given position using given control values
function path = simulate_trajectory(x0, y0, theta0, v, K, constants)

% control trajectories
path.v = v;
path.K = K;

% initial conditions
path.x = x0;
path.y = y0;
path.theta = theta0;
path.t = 0;

path.Ox = [];
path.Oy = [];

% Simulate trajectory
for i=2:numel(K)

    x = [path.x(i-1), path.y(i-1), path.theta(i-1)];
    u = [path.v(i-1), path.K(i-1)];
    
    % RK4 integrator

    k1 = KinematicEquations(x, u);
    k2 = KinematicEquations(x + 0.5*constants.dt*k1, u);
    k3 = KinematicEquations(x + 0.5*constants.dt*k2, u);
    k4 = KinematicEquations(x + constants.dt*k3, u);
    
    x = x + constants.dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    path.x(i) = x(1);
    path.y(i) = x(2);
    path.theta(i) = x(3);
    
    % time
    path.t(i) = path.t(i-1) + constants.dt;
end

% center positions of momentary turning circles
for i=1:numel(K);
    if(path.K(i) == 0)
        path.Ox(i) = path.x(i);
        path.Oy(i) = path.y(i);
    else
        path.Ox(i) = path.x(i) - (1/path.K(i))*sin(path.theta(i));
        path.Oy(i) = path.y(i) + (1/path.K(i))*cos(path.theta(i));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dstate = KinematicEquations(state,control)
    
    v = control(1);
    K = control(2);
    
    x = state(1);
    y = state(2);
    theta = state(3);

    % Kinematic equation
    dx = v * cos(theta);
    dy = v * sin(theta);
    dtheta = v * K;

    dstate = [dx, dy, dtheta];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move given path such that it goes through given position at given time
function path = relocate_path(x, y, theta, path, index)
    
    x_RP = path.x(index);
    y_RP = path.y(index);

    dtheta = theta - path.theta(index);
    R = [cos(dtheta) -sin(dtheta); sin(dtheta) cos(dtheta)];
    
    P = R*[path.x - x_RP; path.y - y_RP];
    path.x = P(1,:) + x;
    path.y = P(2,:) + y;
     
    P = R*[path.Ox - x_RP; path.Oy - y_RP];
    path.Ox = P(1,:) + x;
    path.Oy = P(2,:) + y;
    
    path.theta = path.theta + dtheta;

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move given path such that starting and ending momentary turning points
% equals to given points
function path = relocate_path2(path, Ox1, Oy1, Oxend, Oyend)
    
    x_RP = path.Ox(1);
    y_RP = path.Oy(1);

    ang1 = atan2(path.Oy(end)-path.Oy(1),path.Ox(end)-path.Ox(1));
    ang2 = atan2(Oyend-Oy1,Oxend-Ox1);
            
    dtheta = ang2 - ang1;
    R = [cos(dtheta) -sin(dtheta); sin(dtheta) cos(dtheta)];

    P = R*[path.x - x_RP; path.y - y_RP];
    path.x = P(1,:) + Ox1;
    path.y = P(2,:) + Oy1;
     
    P = R*[path.Ox - x_RP; path.Oy - y_RP];
    path.Ox = P(1,:) + Ox1;
    path.Oy = P(2,:) + Oy1;
    
    path.theta = path.theta + dtheta;

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move given path such that its momentary turning point in the beginning/end 
% (flag = 0 or !0, respectively) equals (Ox1,Oy1).
% In addition, path is rotated so that its ending/beginning (flag = 0 or !0, respectively) heading equals phi.
function path = relocate_path3(path, Ox1, Oy1, phi, flag)
    
    if flag == 0
        x_RP = path.Ox(1);
        y_RP = path.Oy(1);
        dtheta = phi - path.theta(end);
    else
        x_RP = path.Ox(end);
        y_RP = path.Oy(end);
        dtheta = phi - path.theta(1);
    end
    
    R = [cos(dtheta) -sin(dtheta); sin(dtheta) cos(dtheta)];

    P = R*[path.x - x_RP; path.y - y_RP];
    path.x = P(1,:) + Ox1;
    path.y = P(2,:) + Oy1;
     
    P = R*[path.Ox - x_RP; path.Oy - y_RP];
    path.Ox = P(1,:) + Ox1;
    path.Oy = P(2,:) + Oy1;
    
    path.theta = path.theta + dtheta;  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create arch between (x1,y1) and (x2,y2) using (Ox,Oy) as a center point
% and using v as a speed profile 
%
% angle_step = 1 --> counterclockwise
% angle_step = -1 --> clockwise
function [path v1_idx v2_idx] = create_arch(Ox, Oy, x1, y1, x2, y2, angle_step, v1, v2, constants)

v1_idx = 1;
v2_idx = numel(v2);

v_direction = sum([sign(v1) sign(v2)]);

if(abs(v_direction) ~= (numel(v1)+numel(v2)) )
    path = [];
    return;
end

angle_step = angle_step * sign(v_direction);

r = sqrt((x1-Ox)^2+(y1-Oy)^2);
r2 = sqrt((x2-Ox)^2+(y2-Oy)^2);
ang1 = atan2(y1-Oy,x1-Ox);
ang2 = atan2(y2-Oy,x2-Ox);

if(r == 0 || r2 == 0)
    % TODO: check that (x1,y1) == (x2,y2)
    
    path.t = [];
    path.theta = [];
    
    path.angle = 0;
    path.r = 0;
    return;
end
    
if(angle_step > 0 && ang2 < ang1)
    ang2 = ang2 + 2*pi;
end
if(angle_step < 0 && ang2 > ang1)
    ang2 = ang2 - 2*pi;
end
    
angle = angle_step*(ang2 - ang1);

K = sign(v_direction)*angle_step/r;

path1.v = [];
path1.K = [];

path1.x = [];
path1.y = [];
path1.theta = [];
path1.t = [];

path1.Ox = [];
path1.Oy = [];

ang = ang1;
if(numel(v1) > 0)
    % Simulate forwards
    
    while(1)

        w = abs(v1(v1_idx) / r);          % direction is chosen by anglestep
        ang = ang + angle_step*w*constants.dt;

        v1_idx = v1_idx + 1;

        if(v1_idx > numel(v1))
            v1_idx = numel(v1);
            
            if(numel(v2) > 0)
                % continue with second speed vector if exists
                break;
            end
        end

        if(angle_step*(ang2 - ang) <= 0)
            break;
        end

        path1.v(end+1) = v1(v1_idx);
        path1.K(end+1) = K;

        path1.x(end+1) = Ox + r * cos(ang);
        path1.y(end+1) = Oy + r * sin(ang);

        path1.theta(end+1) = ang + sign(v_direction)*angle_step*pi/2;

        path1.Ox(end+1) = Ox;
        path1.Oy(end+1) = Oy;
        
    end
end

path1.t = (0:(numel(path1.v)-1))*constants.dt;


path2.v = [];
path2.K = [];

path2.x = [];
path2.y = [];
path2.theta = [];
path2.t = [];

path2.Ox = [];
path2.Oy = [];


ang1 = ang;
ang = ang2;

if(numel(v2) > 0)
    % simulate bakwards

    while(1)

        w = abs(v2(v2_idx) / r);
        ang = ang - angle_step*w*constants.dt;

        if(angle_step*(ang - ang1) <= 0)
            break;
        end
        
        path2.v = [v2(v2_idx) path2.v];
        path2.K = [K path2.K];

        path2.x = [Ox + r * cos(ang) path2.x];
        path2.y = [Oy + r * sin(ang) path2.y];

        path2.theta = [ang + sign(v_direction)*angle_step*pi/2 path2.theta];
        
        path2.Ox = [Ox path2.Ox];
        path2.Oy = [Oy path2.Oy];
        
        v2_idx = v2_idx - 1;

        if(v2_idx < 1)
            v2_idx = 1;
        end

    end

end

path2.t = (0:(numel(path2.v)-1))*constants.dt;

path = connect_paths(path1,path2);
path.angle = angle;
path.r = r;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% connect two paths 
function path = connect_paths(P1,P2)

if(numel(P1.t) == 0)
    path = P2;
    return;
end

if(numel(P2.t) == 0)
    path = P1;
    return;
end

path.x = [P1.x P2.x];
path.y = [P1.y P2.y];
path.theta = [P1.theta P2.theta];
path.Ox = [P1.Ox P2.Ox];
path.Oy = [P1.Oy P2.Oy];
path.v = [P1.v P2.v];
path.K = [P1.K P2.K];

if(isfield(P1,'ID'))
    path.ID = P1.ID;
end

if(isfield(P2,'ID'))
    path.ID = [path.ID P2.ID];
end

if(isfield(P1,'dK') && isfield(P2,'dK'))
    path.dK = [P1.dK P2.dK];
end

dt = abs(sqrt((P1.x(end)-P2.x(1))^2 + (P1.y(end)-P2.y(1))^2) / P1.v(end));

path.t = [P1.t P2.t+P1.t(end)+dt]; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angle = max_turning_angle(theta)

dtheta = diff(theta);

idx = dtheta > pi;
dtheta(idx) = dtheta(idx) - 2*pi;

idx = dtheta < -pi;
dtheta(idx) = dtheta(idx) + 2*pi;

direction = sign(sum(sign(dtheta)));
dtheta = dtheta(sign(dtheta) == direction);

angle = sum(dtheta);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create line between (x1,y1) and (x2,y2) using v1 and v2 as a speed profile 
% NOTE: v1(end) and v2(1) must be same and space between points must be
% long enough for those

function [path v1_idx v2_idx] = create_straight_line(x1, y1, x2, y2, v1, v2, constants)

v1_idx = 1;
v2_idx = numel(v2);

theta = atan2(y2 - y1, x2 - x1);

path1.v = [];
path1.K = [];

path1.x = [];
path1.y = [];
path1.theta = [];
path1.t = [];

path1.Ox = [];
path1.Oy = [];

pos1 = [x1 y1];
if(numel(v1) > 0)
    % Simulate forwards
    
    % v1_idx = 1;
        
    pos_diff = [x2-x1 y2-y1];
    dir_vec = pos_diff/norm(pos_diff);
    
    while(1)

        pos1 = pos1 + dir_vec*v1(v1_idx)*constants.dt;

        v1_idx = v1_idx + 1;

        if(v1_idx > numel(v1))
            v1_idx = numel(v1);
            
            if(numel(v2) > 0)
                % continue with second speed vector if exists
                break;
            end
        end
        
        if norm(pos1-[x1 y1]) >= norm(pos_diff)
            break;
        end
        
        path1.v(end+1) = v1(v1_idx);
        path1.K(end+1) = 0;

        path1.x(end+1) = pos1(1);
        path1.y(end+1) = pos1(2);

        path1.theta(end+1) = theta;

        path1.Ox(end+1) = pos1(1);
        path1.Oy(end+1) = pos1(2);

    end
end

path1.t = (0:(numel(path1.v)-1))*constants.dt;


path2.v = [];
path2.K = [];

path2.x = [];
path2.y = [];
path2.theta = [];
path2.t = [];

path2.Ox = [];
path2.Oy = [];

pos2 = [x2 y2];
if(numel(v2) > 0)
    % Simulate backwards
    
    % v2_idx = numel(v2);
    
    pos_diff = pos2-pos1;
    dir_vec = -pos_diff/norm(pos_diff);
    
    while(1)

        pos2 = pos2 + dir_vec*v2(v2_idx)*constants.dt;
        
        if norm(pos2-[x2 y2]) >= norm(pos_diff)
            break;
        end
        
        path2.v = [v2(v2_idx) path2.v];
        path2.K = [0 path2.K];

        path2.x = [pos2(1) path2.x];
        path2.y = [pos2(2) path2.y];

        path2.theta = [theta path2.theta];

        path2.Ox = [pos2(1) path2.Ox];
        path2.Oy = [pos2(2) path2.Oy];

        v2_idx = v2_idx - 1;

        if(v2_idx < 1)
            v2_idx = 1;
        end

    end
end

path2.t = (0:(numel(path2.v)-1))*constants.dt;

path = connect_paths(path1,path2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

