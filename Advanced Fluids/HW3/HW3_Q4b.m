function HW3_Q4b()

set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultFigureColor', 'white');

    % Define constants
    a = 1.0;  % Constant that I don't understand
    mu = 1.83e-5; % Dynamic viscosity, in Pa-s
    rho = 1000; % Density of the particle, kg/m^3
    g = [0; -9.81]; % Gravity vector (assuming 2D motion)

    % Define parameter variations
    r_values = [1e-5, 5e-5, 10e-5]; % Different particle radii in m
    up0_values = [[10; 0], [20; 0], [30; 0]]; % Different initial velocities in m/s

    % Time span
    tspan = [0, 100]; % Solve from t=0 to t=10

    figure; hold on; % Prepare figure for plotting
    
    % Loop over different values of r and up0
    for i = 1:length(r_values)
        for j = 1:size(up0_values, 2)
            r = r_values(i);
            up0 = up0_values(:, j);
            xp0 = [0; 0];  % Initial position (fixed)

            initial_conditions = [xp0; up0];
            options = odeset('Events', @(t, y) event_func(t, y)); % Set event function

            % Solve ODE (FIXED: now using 'options')
            [t, sol] = ode45(@(t, y) odefunc(t, y, a, mu, r, rho, g), tspan, initial_conditions, options);

            % Extract position
            xp = sol(:, 1:2);

            % Plot trajectory
            plot(xp(:,1), xp(:,2), 'DisplayName', sprintf('$r = %.1e, \\; \\mathbf{u}_{p0} = [%d, %d]$', r, up0(1), up0(2)));
        end
    end

    % Plot formatting
    xlabel('X Position'); ylabel('Y Position');
    title('Particle Trajectories');
    axis equal; grid on;
    legend show;
end

function dydt = odefunc(t, y, a, mu, r, rho, g)
    % Extract variables
    xp = y(1:2); % Position vector
    up = y(3:4); % Velocity vector

    % Compute fluid velocity using the given function
    eta = norm(xp);
    uf = (1 / (1 + a * eta^2)^2) * [1; 0]; % Assuming flow is in the x-direction

    % Compute acceleration
    dupdt = (9 * mu / (2 * r^2 * rho)) * (uf - up) + g;

    % Construct dydt
    dydt = [up; dupdt];
end

function [value, isterminal, direction] = event_func(~, y)
    % Event function to stop integration when y-position reaches -0.5
    value = y(2) + 0.5;  % y reaches -0.5 when value = 0
    isterminal = 1; % Stop the integration
    direction = -1;  % Detect when y is decreasing through -0.5
end
