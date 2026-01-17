%% Plot Convergence Comparison: Gradient Descent vs PSO vs PSO-TVIW
% This script loads convergence history from all three algorithms and
% creates a comparison plot showing Sum Rate vs Iterations

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('  Loading Convergence Data...          \n');
fprintf('========================================\n');

%% Load Gradient Descent convergence
try
    data = load('GD_convergence.mat');
    gd_conv = data.convergence_history;
    gd_time = data.convergence_time;
    fprintf('✓ Loaded Gradient Descent convergence\n');
    has_GD = true;
catch
    fprintf('✗ Gradient Descent convergence not found\n');
    has_GD = false;
end

%% Load Standard PSO convergence
try
    data = load('PSO_convergence.mat');
    pso_conv = data.convergence_history;
    pso_time = data.convergence_time;
    fprintf('✓ Loaded Standard PSO convergence\n');
    has_PSO = true;
catch
    fprintf('✗ Standard PSO convergence not found\n');
    has_PSO = false;
end

%% Load PSO-TVIW convergence
try
    data = load('PSO_TVIW_convergence.mat');
    tviw_conv = data.convergence_history;
    tviw_time = data.convergence_time;
    fprintf('✓ Loaded PSO-TVIW convergence\n');
    has_TVIW = true;
catch
    fprintf('✗ PSO-TVIW convergence not found\n');
    has_TVIW = false;
end

fprintf('========================================\n\n');

%% Load Optimal Capacity (if available)
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

%% Check if at least one algorithm has convergence data
if ~has_GD && ~has_PSO && ~has_TVIW
    error('No convergence data found! Please run the algorithms first.');
end

%% Create Convergence Comparison Plot
figure('Position', [100, 100, 900, 600]);
hold on;
grid on;

if has_GD
    plot(1:length(gd_conv), gd_conv, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 6, ...
        'MarkerIndices', 1:max(1,floor(length(gd_conv)/10)):length(gd_conv), ...
        'DisplayName', 'Gradient Descent');
end

if has_PSO
    plot(1:length(pso_conv), pso_conv, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 6, ...
        'MarkerIndices', 1:max(1,floor(length(pso_conv)/10)):length(pso_conv), ...
        'DisplayName', 'PSO (w=0.7)');
end

if has_TVIW
    plot(1:length(tviw_conv), tviw_conv, 'g-^', 'LineWidth', 2.5, 'MarkerSize', 6, ...
        'MarkerIndices', 1:max(1,floor(length(tviw_conv)/10)):length(tviw_conv), ...
        'DisplayName', 'PSO-TVIW (w: 0.9→0.4)');
end

xlabel('Iterations', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 14, 'FontWeight', 'bold');
title('Convergence Comparison: GD vs PSO vs PSO-TVIW', 'FontSize', 15, 'FontWeight', 'bold');

%% Add optimal capacity reference line
if has_optimal
    xlim_current = xlim;
    plot(xlim_current, [optimal_capacity optimal_capacity], 'k--', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Optimal (%.2f)', optimal_capacity));
end

legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;

%% Save figure
savefig('convergence_comparison.fig');
saveas(gcf, 'convergence_comparison.jpg');
fprintf('✓ Convergence plot saved to:\n');
fprintf('  - convergence_comparison.fig\n');
fprintf('  - convergence_comparison.jpg\n\n');

%% Print convergence statistics
fprintf('========================================\n');
fprintf('  Convergence Statistics                \n');
fprintf('========================================\n\n');

if has_GD
    fprintf('Gradient Descent:\n');
    fprintf('  Iterations:      %d\n', length(gd_conv));
    fprintf('  Initial Capacity: %.4f bits/s/Hz\n', gd_conv(1));
    fprintf('  Final Capacity:   %.4f bits/s/Hz\n', gd_conv(end));
    fprintf('  Improvement:      %.4f bits/s/Hz (%.2f%%)\n\n', ...
        gd_conv(end)-gd_conv(1), (gd_conv(end)-gd_conv(1))/gd_conv(1)*100);
end

if has_PSO
    fprintf('Standard PSO (w=0.7):\n');
    fprintf('  Iterations:       %d\n', length(pso_conv));
    fprintf('  Initial Capacity: %.4f bits/s/Hz\n', pso_conv(1));
    fprintf('  Final Capacity:   %.4f bits/s/Hz\n', pso_conv(end));
    fprintf('  Improvement:      %.4f bits/s/Hz (%.2f%%)\n\n', ...
        pso_conv(end)-pso_conv(1), (pso_conv(end)-pso_conv(1))/pso_conv(1)*100);
end

if has_TVIW
    fprintf('PSO-TVIW (w: 0.9→0.4):\n');
    fprintf('  Iterations:       %d\n', length(tviw_conv));
    fprintf('  Initial Capacity: %.4f bits/s/Hz\n', tviw_conv(1));
    fprintf('  Final Capacity:   %.4f bits/s/Hz\n', tviw_conv(end));
    fprintf('  Improvement:      %.4f bits/s/Hz (%.2f%%)\n\n', ...
        tviw_conv(end)-tviw_conv(1), (tviw_conv(end)-tviw_conv(1))/tviw_conv(1)*100);
end

fprintf('========================================\n');
