clear;

%% Define constants and parameters
capital_share = 0.25;  % Share of capital in output
savings_rate  = 0.2;  % Proportion of output saved
deprec_rate   = 0.15;  % Depreciation of capital
pop_growth    = 0.01; % Growth rate of population
tech_progress = 0.02; % Technological advancement rate

%% Define a function to compute capital stock in the next period
next_capital = @(cap) (1/((1 + pop_growth) * (1 + tech_progress))) * (savings_rate * cap.^capital_share + (1 - deprec_rate) * cap);

%% Set range for capital stock and plot capital dynamics
cap_min = 0.01;
cap_max = 3;
cap_vals = linspace(cap_min, cap_max, 300);  % Create a grid of capital values

next_cap_vals = next_capital(cap_vals);

figure;
plot(cap_vals, cap_vals, 'LineWidth', 1.5); % Plot 45-degree line
hold on;
plot(cap_vals, next_cap_vals, 'LineWidth', 2); % Plot k(t+1) vs k(t)
hold off;
xlabel('$\hat{k}(t)$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\hat{k}(t+1)$', 'Interpreter', 'latex', 'FontSize', 12);
legend('45-degree line', '$\hat{k}(t)$', 'Interpreter', 'latex', 'Location', 'NorthWest', 'FontSize', 14);
title('Capital Dynamics in the Modified Solow Model', 'FontSize', 14);

%% Plot investment versus break-even investment
real_investment = savings_rate * cap_vals.^capital_share;
threshold_investment = (deprec_rate + pop_growth + tech_progress + pop_growth*tech_progress) * cap_vals;
output_vals = cap_vals.^capital_share;

figure;
plot(cap_vals, real_investment, 'LineWidth', 2);  % Real investment
hold on;
plot(cap_vals, threshold_investment, 'LineWidth', 2);  % Threshold investment
plot(cap_vals, output_vals, 'LineWidth', 2);  % Output
hold off;
xlabel('$\hat{k}(t)$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Real investment', 'Break-even investment', 'Output', 'Location', 'NorthWest', 'FontSize', 14);
title('Investment and Output in the Solow Model', 'FontSize', 14);

%% Solve for steady-state capital
options = optimset('Display', 'iter');  % Display iterative process
[steady_cap, fval, exitflag] = fzero(@(k) capital_root(k, pop_growth, tech_progress, savings_rate, capital_share, deprec_rate), [cap_min, cap_max], options);

%% Simulate convergence of capital over time
max_iter = 300;      % Maximum number of iterations
tolerance = 1e-5;    % Convergence tolerance

capital_stock = nan(max_iter, 1);  % Preallocate array for capital stock
capital_stock(1) = 0.2;            % Initial value for capital

t = 1;
while t < max_iter
    capital_stock(t+1) = next_capital(capital_stock(t));  % Calculate next period capital
    diff_cap = abs(capital_stock(t+1) - capital_stock(t));  % Calculate difference
    if diff_cap < tolerance
        capital_stock(t+2:end) = capital_stock(t+1);  % If converged, fill rest with steady-state
        break;
    end
    t = t + 1;
end

% Plot capital convergence path
figure;
plot(1:t, capital_stock(1:t), 'LineWidth', 2);
xlabel('Time (t)', 'FontSize', 12);
ylabel('$\hat{k}(t)$', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\hat{k}(t)$', 'Location', 'NorthWest', 'Interpreter', 'latex', 'FontSize', 14);
title('Convergence Path of Capital Stock', 'Interpreter', 'latex', 'FontSize', 14);

%% Function to compute the capital root (for finding the steady state)
function f = capital_root(cap, pop_rate, tech_growth, save_rate, cap_share, depreciation)
    future_cap = (1 / ((1 + pop_rate) * (1 + tech_growth))) * (save_rate * cap.^cap_share + (1 - depreciation) * cap);
    f = cap - future_cap;
end
