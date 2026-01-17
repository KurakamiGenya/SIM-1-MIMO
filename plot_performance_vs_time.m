%% Plot Performance vs Execution Time
% This script creates a scatter plot showing the trade-off between
% final Sum Rate (performance) and execution time for all algorithms

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('  Loading Performance and Time Data... \n');
fprintf('========================================\n');

%% Load Gradient Descent results
try
    load('Capacity_K_10.mat');
    gd_capacity = Capacity_K_10(10);  % L=10 (Max_L)
    load('GD_execution_time.mat');
    gd_time = algorithm_execution_time;
    fprintf('✓ Loaded Gradient Descent: Capacity=%.4f, Time=%.2fs\n', gd_capacity, gd_time);
    has_GD = true;
catch
    fprintf('✗ Gradient Descent data not found\n');
    has_GD = false;
end

%% Load Standard PSO results
try
    load('Capacity_K_10_PSO.mat');
    pso_capacity = Capacity_K_10_PSO(10);  % L=10 (Max_L)
    load('PSO_execution_time.mat');
    pso_time = algorithm_execution_time;
    fprintf('✓ Loaded Standard PSO: Capacity=%.4f, Time=%.2fs\n', pso_capacity, pso_time);
    has_PSO = true;
catch
    fprintf('✗ Standard PSO data not found\n');
    has_PSO = false;
end

%% Load PSO-TVIW results
try
    load('Capacity_K_10_PSO_TVIW.mat');
    tviw_capacity = Capacity_K_10_PSO_TVIW(10);  % L=10 (Max_L)
    load('PSO_TVIW_execution_time.mat');
    tviw_time = algorithm_execution_time;
    fprintf('✓ Loaded PSO-TVIW: Capacity=%.4f, Time=%.2fs\n', tviw_capacity, tviw_time);
    has_TVIW = true;
catch
    fprintf('✗ PSO-TVIW data not found\n');
    has_TVIW = false;
end

%% Load Optimal Capacity (optional, for reference)
try
    load('Capacity_max.mat');
    optimal_capacity = Capacity_max(10);  % L=10 (Max_L)
    fprintf('✓ Loaded Optimal Capacity: %.4f bits/s/Hz\n', optimal_capacity);
    has_optimal = true;
catch
    fprintf('ℹ Optimal capacity not found (run test_max.m to generate)\n');
    has_optimal = false;
end

fprintf('========================================\n\n');

%% Check if at least one algorithm has data
if ~has_GD && ~has_PSO && ~has_TVIW
    error('No algorithm data found! Please run the algorithms first.');
end

%% Create Performance vs Time Plot
figure('Position', [100, 100, 900, 700]);
hold on;
grid on;

% Plot each algorithm as a scatter point with labels
marker_size = 200;

if has_GD
    scatter(gd_time, gd_capacity, marker_size, 'b', 'o', 'filled', ...
        'DisplayName', 'Gradient Descent', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(gd_time, gd_capacity, '  GD', 'FontSize', 11, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

if has_PSO
    scatter(pso_time, pso_capacity, marker_size, 'r', 's', 'filled', ...
        'DisplayName', 'PSO (w=0.7)', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(pso_time, pso_capacity, '  PSO', 'FontSize', 11, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

if has_TVIW
    scatter(tviw_time, tviw_capacity, marker_size, 'g', '^', 'filled', ...
        'DisplayName', 'PSO-TVIW (w: 0.9→0.4)', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(tviw_time, tviw_capacity, '  PSO-TVIW', 'FontSize', 11, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

%% Add optimal capacity reference line (horizontal)
if has_optimal
    xlim_current = xlim;
    plot(xlim_current, [optimal_capacity optimal_capacity], 'k--', 'LineWidth', 2, ...
        'DisplayName', sprintf('Optimal (%.2f)', optimal_capacity));
    text(mean(xlim_current), optimal_capacity, sprintf('  Optimal: %.2f', optimal_capacity), ...
        'FontSize', 10, 'VerticalAlignment', 'bottom', 'Color', 'k');
end

xlabel('Execution Time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 14, 'FontWeight', 'bold');
title('Performance vs Execution Time Trade-off', 'FontSize', 15, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;

%% Save figure
savefig('performance_vs_time.fig');
saveas(gcf, 'performance_vs_time.jpg');
fprintf('✓ Performance vs Time plot saved to:\n');
fprintf('  - performance_vs_time.fig\n');
fprintf('  - performance_vs_time.jpg\n\n');

%% Print summary statistics
fprintf('========================================\n');
fprintf('  Performance vs Time Summary          \n');
fprintf('========================================\n\n');

if has_GD && has_PSO && has_TVIW
    % Calculate efficiency metrics
    fprintf('Algorithm Comparison:\n');
    fprintf('--------------------\n');
    
    % Best performance
    [best_capacity, best_idx] = max([gd_capacity, pso_capacity, tviw_capacity]);
    alg_names = {'Gradient Descent', 'PSO', 'PSO-TVIW'};
    fprintf('Best Performance: %s (%.4f bits/s/Hz)\n', alg_names{best_idx}, best_capacity);
    
    % Fastest algorithm
    [min_time, fast_idx] = min([gd_time, pso_time, tviw_time]);
    fprintf('Fastest Algorithm: %s (%.2f seconds)\n', alg_names{fast_idx}, min_time);
    
    % Efficiency (Capacity per second)
    gd_efficiency = gd_capacity / gd_time;
    pso_efficiency = pso_capacity / pso_time;
    tviw_efficiency = tviw_capacity / tviw_time;
    
    fprintf('\nEfficiency (Capacity/Time):\n');
    fprintf('  GD:       %.6f bits/s/Hz per second\n', gd_efficiency);
    fprintf('  PSO:      %.6f bits/s/Hz per second\n', pso_efficiency);
    fprintf('  PSO-TVIW: %.6f bits/s/Hz per second\n', tviw_efficiency);
    
    [best_eff, eff_idx] = max([gd_efficiency, pso_efficiency, tviw_efficiency]);
    fprintf('  → Most Efficient: %s\n', alg_names{eff_idx});
    
    if has_optimal
        fprintf('\nGap to Optimal:\n');
        fprintf('  GD:       %.4f bits/s/Hz (%.2f%%)\n', optimal_capacity - gd_capacity, ...
            (optimal_capacity - gd_capacity)/optimal_capacity*100);
        fprintf('  PSO:      %.4f bits/s/Hz (%.2f%%)\n', optimal_capacity - pso_capacity, ...
            (optimal_capacity - pso_capacity)/optimal_capacity*100);
        fprintf('  PSO-TVIW: %.4f bits/s/Hz (%.2f%%)\n', optimal_capacity - tviw_capacity, ...
            (optimal_capacity - tviw_capacity)/optimal_capacity*100);
    end
end

fprintf('\n========================================\n');
