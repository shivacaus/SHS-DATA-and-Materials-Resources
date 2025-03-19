% Robust SEIR Model with Error Handling and Optimization
format long;

% Load the data
data = readtable('ssss2.xlsx');
time = data.time;
S = data.S;
E = data.E;
I = data.I;
R = data.R;
temperature = data.temperature;
humidity = data.humidity;
windspeed = data.windspeed;

% Smooth the data using a moving average
windowSize = 5; % Configurable window size
S = smoothdata(S, 'movmean', windowSize);
E = smoothdata(E, 'movmean', windowSize);
I = smoothdata(I, 'movmean', windowSize);
R = smoothdata(R, 'movmean', windowSize);
temperature = smoothdata(temperature, 'movmean', windowSize);
humidity = smoothdata(humidity, 'movmean', windowSize);
windspeed = smoothdata(windspeed, 'movmean', windowSize);

% Total population
N = S + E + I + R;

% Normalize the SEIR compartments
s = S ./ N;
e = E ./ N;
i = I ./ N;
r = R ./ N;

% Normalize environmental factors using z-score
temp_norm = zscore(temperature);
humidity_norm = zscore(humidity);
wind_norm = zscore(windspeed);

% Combine all data into a matrix
data_combined = [s, e, i, r, temp_norm, humidity_norm, wind_norm];

% Define parameter bounds
lb = [0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0];
ub = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

% Initial guess for parameters (optional)
% x0 = [0.29, 0.2, 0.19, 0.01, 0.01, -0.01, -0.01, 0.01, -0.01, -0.01, -0.01, 0.01, 0.01, -0.01, 0.1];

% Set up the optimization options
options = optimoptions('particleswarm', ...
    'Display', 'iter', ...
    'MaxIterations', 10, ... % Increased iterations
    'SwarmSize', 100, ...    % Adjusted swarm size
    'FunctionTolerance', 1e-6, ...
    'MaxStallIterations', 50); % Reduced stall iterations

% Perform the optimization using particle swarm
objective = @(x) robust_objective(x, time, data_combined);
[x, fval] = particleswarm(objective, length(lb), lb, ub, options);

% Extract the optimized parameters
beta0 = x(1);
epsilon0 = x(2);
rho0 = x(3);
a1 = x(4); a2 = x(5); % Seasonal parameters
b_t = x(6); b_h = x(7); b_w = x(8); % Temperature, humidity, wind effects on beta
e_t = x(9); e_h = x(10); e_w = x(11); % Effects on epsilon
r_t = x(12); r_h = x(13); r_w = x(14); % Effects on rho
intervention_strength = x(15);

% Run the model with optimized parameters
[s_fit, e_fit, i_fit, r_fit] = seir_improved_model(x, time, data_combined);

% Display the optimized parameters
display_parameters(x, temp_norm, humidity_norm, wind_norm);

% Calculate goodness of fit metrics
[rmse_s, rmse_e, rmse_i, rmse_r, mae_s, mae_e, mae_i, mae_r, r2_s, r2_e, r2_i, r2_r] = calculate_goodness_of_fit(s, e, i, r, s_fit, e_fit, i_fit, r_fit);

% Display goodness of fit metrics
fprintf('Goodness of Fit Metrics:\n');
fprintf('RMSE (S): %.4f, RMSE (E): %.4f, RMSE (I): %.4f, RMSE (R): %.4f\n', rmse_s, rmse_e, rmse_i, rmse_r);
fprintf('MAE (S): %.4f, MAE (E): %.4f, MAE (I): %.4f, MAE (R): %.4f\n', mae_s, mae_e, mae_i, mae_r);
fprintf('R-squared (S): %.4f, R-squared (E): %.4f, R-squared (I): %.4f, R-squared (R): %.4f\n', r2_s, r2_e, r2_i, r2_r);

% Plot the results
plot_results(time, s, e, i, r, s_fit, e_fit, i_fit, r_fit, rmse_s, rmse_e, rmse_i, rmse_r, r2_s, r2_e, r2_i, r2_r);

% Supporting Functions

function obj_value = robust_objective(params, t, data_combined)
    try
        [s_fit, e_fit, i_fit, r_fit] = seir_improved_model(params, t, data_combined);
        residuals = [s_fit - data_combined(:,1); 
                     e_fit - data_combined(:,2); 
                     i_fit - data_combined(:,3); 
                     r_fit - data_combined(:,4)];
        obj_value = sum(residuals.^2);
    catch ME
        fprintf('Error in robust_objective: %s\n', ME.message); % Log error
        obj_value = 1e10; % Return a large value if an error occurs
    end
end

function [s, e, i, r] = seir_improved_model(params, t, data_combined)
    beta0 = params(1);
    epsilon0 = params(2);
    rho0 = params(3);
    a1 = params(4); a2 = params(5);
    b_t = params(6); b_h = params(7); b_w = params(8);
    e_t = params(9); e_h = params(10); e_w = params(11);
    r_t = params(12); r_h = params(13); r_w = params(14);
    intervention_strength = params(15);
    
    omega = 2.13e-5; % Given constant
    alpha = 5.59e-5;

    s0 = data_combined(1,1); e0 = data_combined(1,2); i0 = data_combined(1,3); r0 = data_combined(1,4);
    temp = data_combined(:,5); hum = data_combined(:,6); wind = data_combined(:,7);
    
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1);
    [t_sol, y] = ode45(@(t,y) seir_ode(t, y, beta0, epsilon0, rho0, alpha, omega, a1, a2, b_t, b_h, b_w, e_t, e_h, e_w, r_t, r_h, r_w, intervention_strength, temp, hum, wind), t, [s0; e0; i0; r0], options);
    
    % Interpolate the solution to match the original time points
    s = interp1(t_sol, y(:,1), t, 'pchip');
    e = interp1(t_sol, y(:,2), t, 'pchip');
    i = interp1(t_sol, y(:,3), t, 'pchip');
    r = interp1(t_sol, y(:,4), t, 'pchip');
end

function dydt = seir_ode(t, y, beta0, epsilon0, rho0, alpha, omega, a1, a2, b_t, b_h, b_w, e_t, e_h, e_w, r_t, r_h, r_w, intervention_strength, temp, hum, wind)
    s = max(0, min(1, y(1)));
    e = max(0, min(1, y(2)));
    i = max(0, min(1, y(3)));
    r = max(0, min(1, y(4)));
    
    temp_t = interp1(1:length(temp), temp, t, 'pchip', 'extrap');
    hum_t = interp1(1:length(hum), hum, t, 'pchip', 'extrap');
    wind_t = interp1(1:length(wind), wind, t, 'pchip', 'extrap');
    
    % Time-varying parameters
    beta = beta0 * (1 + a1*cos(2*pi*t/365) + a2*sin(2*pi*t/365)) * ...
           (1 + b_t * temp_t + b_h * hum_t + b_w * wind_t) * ...
           exp(-intervention_strength * t/100);

    epsilon = epsilon0 * (1 + e_t * temp_t + e_h * hum_t + e_w * wind_t) * ...
              exp(-intervention_strength * t/150);

    rho = rho0 * (1 + r_t * temp_t + r_h * hum_t + r_w * wind_t) * ...
           exp(-intervention_strength * t/200);

    dsdt = alpha - beta * s * i - omega * s;
    dedt = beta * s * i - epsilon * e - omega * e;
    didt = epsilon * e - rho * i - omega * i;
    drdt = rho * i - omega * r;
    
    dydt = [dsdt; dedt; didt; drdt];
end

function display_parameters(params, temp, humidity, wind)
    beta0 = params(1);
    epsilon0 = params(2);
    rho0 = params(3);
    b_t = params(6); b_h = params(7); b_w = params(8);
    e_t = params(9); e_h = params(10); e_w = params(11);
    r_t = params(12); r_h = params(13); r_w = params(14);
    intervention_strength = params(15);
    
    fprintf('Optimized Parameters:\n');
    fprintf('Beta0: %.4f\n', beta0);
    fprintf('Epsilon0: %.4f\n', epsilon0);
    fprintf('Rho0: %.4f\n', rho0);
    fprintf('Intervention Strength: %.4f\n\n', intervention_strength);
    
    fprintf('Environmental Effects on Beta:\n');
    fprintf('Temperature: %.4f, Humidity: %.4f, Wind: %.4f\n\n', b_t, b_h, b_w);
    
    fprintf('Environmental Effects on Epsilon:\n');
    fprintf('Temperature: %.4f, Humidity: %.4f, Wind: %.4f\n\n', e_t, e_h, e_w);
    
    fprintf('Environmental Effects on Rho:\n');
    fprintf('Temperature: %.4f, Humidity: %.4f, Wind: %.4f\n\n', r_t, r_h, r_w);
    
    % Plot effects of environmental factors
    plot_environmental_effects(params, temp, humidity, wind);
end

function plot_environmental_effects(params, temp, humidity, wind)
    b_t = params(6); b_h = params(7); b_w = params(8);
    e_t = params(9); e_h = params(10); e_w = params(11);
    r_t = params(12); r_h = params(13); r_w = params(14);
    
    figure;
    subplot(3,1,1);
    plot(temp, 1 + b_t * temp, 'r-', humidity, 1 + b_h * humidity, 'g-', wind, 1 + b_w * wind, 'b-');
    title('Environmental Effects on Beta');
    legend('Temperature', 'Humidity', 'Wind');
    ylabel('Multiplier');
    
    subplot(3,1,2);
    plot(temp, 1 + e_t * temp, 'r-', humidity, 1 + e_h * humidity, 'g-', wind, 1 + e_w * wind, 'b-');
    title('Environmental Effects on Epsilon');
    legend('Temperature', 'Humidity', 'Wind');
    ylabel('Multiplier');
    
    subplot(3,1,3);
    plot(temp, 1 + r_t * temp, 'r-', humidity, 1 + r_h * humidity, 'g-', wind, 1 + r_w * wind, 'b-');
    title('Environmental Effects on Rho');
    legend('Temperature', 'Humidity', 'Wind');
    ylabel('Multiplier');
    xlabel('Normalized Environmental Factor');
end

function [rmse_s, rmse_e, rmse_i, rmse_r, mae_s, mae_e, mae_i, mae_r, r2_s, r2_e, r2_i, r2_r] = calculate_goodness_of_fit(s, e, i, r, s_fit, e_fit, i_fit, r_fit)
    rmse_s = sqrt(mean((s_fit - s).^2));
    rmse_e = sqrt(mean((e_fit - e).^2));
    rmse_i = sqrt(mean((i_fit - i).^2));
    rmse_r = sqrt(mean((r_fit - r).^2));

    mae_s = mean(abs(s_fit - s));
    mae_e = mean(abs(e_fit - e));
    mae_i = mean(abs(i_fit - i));
    mae_r = mean(abs(r_fit - r));

    r2_s = 1 - sum((s_fit - s).^2) / sum((s - mean(s)).^2);
    r2_e = 1 - sum((e_fit - e).^2) / sum((e - mean(e)).^2);
    r2_i = 1 - sum((i_fit - i).^2) / sum((i - mean(i)).^2);
    r2_r = 1 - sum((r_fit - r).^2) / sum((r - mean(r)).^2);
end

function plot_results(time, s, e, i, r, s_fit, e_fit, i_fit, r_fit, rmse_s, rmse_e, rmse_i, rmse_r, r2_s, r2_e, r2_i, r2_r)
    figure;
    subplot(2,2,1);
    plot(time, s, 'b.', time, s_fit, 'r-');
    title(sprintf('Susceptible (RMSE: %.4f, R^2: %.4f)', rmse_s, r2_s));
    legend('Data', 'Fit');

    subplot(2,2,2);
    plot(time, e, 'b.', time, e_fit, 'r-');
    title(sprintf('Exposed (RMSE: %.4f, R^2: %.4f)', rmse_e, r2_e));
    legend('Data', 'Fit');

    subplot(2,2,3);
    plot(time, i, 'b.', time, i_fit, 'r-');
    title(sprintf('Infected (RMSE: %.4f, R^2: %.4f)', rmse_i, r2_i));
    legend('Data', 'Fit');

    subplot(2,2,4);
    plot(time, r, 'b.', time, r_fit, 'r-');
    title(sprintf('Recovered (RMSE: %.4f, R^2: %.4f)', rmse_r, r2_r));
    legend('Data', 'Fit');
end
