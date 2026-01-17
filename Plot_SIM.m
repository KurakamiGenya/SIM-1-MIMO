%% Plot Comparison: Gradient Descent vs PSO vs PSO-TVIW
% This script loads and compares results from three optimization algorithms:
% 1. Gradient Descent (GD)
% 2. Standard PSO (constant inertia weight w=0.7)
% 3. PSO with Time-Varying Inertia Weight (PSO-TVIW, w: 0.9→0.4)

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('  Loading Algorithm Results...         \n');
fprintf('========================================\n');

%% Load Gradient Descent results
try
    load('NMSE_K_10.mat');
    load('Capacity_K_10.mat');
    fprintf('✓ Loaded Gradient Descent results\n');
    has_GD = true;
catch
    fprintf('✗ Gradient Descent results not found\n');
    has_GD = false;
end

%% Load Standard PSO results
try
    load('NMSE_K_10_PSO.mat');
    load('Capacity_K_10_PSO.mat');
    fprintf('✓ Loaded Standard PSO results\n');
    has_PSO = true;
catch
    fprintf('✗ Standard PSO results not found\n');
    has_PSO = false;
end

%% Load PSO-TVIW results
try
    load('NMSE_K_10_PSO_TVIW.mat');
    load('Capacity_K_10_PSO_TVIW.mat');
    fprintf('✓ Loaded PSO-TVIW results\n');
    has_TVIW = true;
catch
    fprintf('✗ PSO-TVIW results not found\n');
    has_TVIW = false;
end

fprintf('========================================\n\n');

%% Check if at least one algorithm has results
if ~has_GD && ~has_PSO && ~has_TVIW
    error('No algorithm results found! Please run the algorithms first.');
end

%% Determine Max_L
if has_GD
    Max_L = length(NMSE_K_10);
elseif has_PSO
    Max_L = length(NMSE_K_10_PSO);
else
    Max_L = length(NMSE_K_10_PSO_TVIW);
end

%% Figure 1: NMSE Comparison
figure('Position', [100, 100, 800, 600]);
hold on;
grid on;

if has_GD
    plot(1:Max_L, NMSE_K_10, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Gradient Descent');
end

if has_PSO
    plot(1:Max_L, NMSE_K_10_PSO, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', 'DisplayName', 'PSO (w=0.7)');
end

if has_TVIW
    plot(1:Max_L, NMSE_K_10_PSO_TVIW, 'g-^', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'g', 'DisplayName', 'PSO-TVIW (w: 0.9→0.4)');
end

xlabel('Number of TX-SIM Layers (L)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('NMSE', 'FontSize', 14, 'FontWeight', 'bold');
title('NMSE Comparison: GD vs PSO vs PSO-TVIW', 'FontSize', 15, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;

%% Figure 2: Capacity Comparison
figure('Position', [150, 150, 800, 600]);
hold on;
grid on;

if has_GD
    plot(1:Max_L, Capacity_K_10, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Gradient Descent');
end

if has_PSO
    plot(1:Max_L, Capacity_K_10_PSO, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', 'DisplayName', 'PSO (w=0.7)');
end

if has_TVIW
    plot(1:Max_L, Capacity_K_10_PSO_TVIW, 'g-^', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'g', 'DisplayName', 'PSO-TVIW (w: 0.9→0.4)');
end

xlabel('Number of TX-SIM Layers (L)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 14, 'FontWeight', 'bold');
title('Capacity Comparison: GD vs PSO vs PSO-TVIW', 'FontSize', 15, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;

%% Print numerical comparison table
fprintf('========================================\n');
fprintf('  Numerical Results (L=1 to L=%d)      \n', Max_L);
fprintf('========================================\n\n');

fprintf('CAPACITY (bits/s/Hz):\n');
fprintf('%-5s | %-15s | %-15s | %-15s\n', 'L', 'GD', 'PSO', 'PSO-TVIW');
fprintf('------|-----------------|-----------------|------------------\n');
for i = 1:Max_L
    fprintf('%-5d | ', i);
    if has_GD
        fprintf('%-15.4f | ', Capacity_K_10(i));
    else
        fprintf('%-15s | ', 'N/A');
    end
    if has_PSO
        fprintf('%-15.4f | ', Capacity_K_10_PSO(i));
    else
        fprintf('%-15s | ', 'N/A');
    end
    if has_TVIW
        fprintf('%-15.4f\n', Capacity_K_10_PSO_TVIW(i));
    else
        fprintf('%-15s\n', 'N/A');
    end
end

fprintf('\nNMSE:\n');
fprintf('%-5s | %-15s | %-15s | %-15s\n', 'L', 'GD', 'PSO', 'PSO-TVIW');
fprintf('------|-----------------|-----------------|------------------\n');
for i = 1:Max_L
    fprintf('%-5d | ', i);
    if has_GD
        fprintf('%-15.6f | ', NMSE_K_10(i));
    else
        fprintf('%-15s | ', 'N/A');
    end
    if has_PSO
        fprintf('%-15.6f | ', NMSE_K_10_PSO(i));
    else
        fprintf('%-15s | ', 'N/A');
    end
    if has_TVIW
        fprintf('%-15.6f\n', NMSE_K_10_PSO_TVIW(i));
    else
        fprintf('%-15s\n', 'N/A');
    end
end

%% Performance comparison at Max_L
fprintf('\n========================================\n');
fprintf('  Performance at L=%d                   \n', Max_L);
fprintf('========================================\n');

if has_GD
    fprintf('Gradient Descent:\n');
    fprintf('  Capacity: %.4f bits/s/Hz\n', Capacity_K_10(Max_L));
    fprintf('  NMSE:     %.6f\n\n', NMSE_K_10(Max_L));
end

if has_PSO
    fprintf('Standard PSO (w=0.7):\n');
    fprintf('  Capacity: %.4f bits/s/Hz\n', Capacity_K_10_PSO(Max_L));
    fprintf('  NMSE:     %.6f\n\n', NMSE_K_10_PSO(Max_L));
end

if has_TVIW
    fprintf('PSO-TVIW (w: 0.9→0.4):\n');
    fprintf('  Capacity: %.4f bits/s/Hz\n', Capacity_K_10_PSO_TVIW(Max_L));
    fprintf('  NMSE:     %.6f\n\n', NMSE_K_10_PSO_TVIW(Max_L));
end

%% Calculate improvements
if has_GD && has_PSO
    cap_imp_pso = ((Capacity_K_10_PSO(Max_L) - Capacity_K_10(Max_L)) / Capacity_K_10(Max_L)) * 100;
    nmse_imp_pso = ((NMSE_K_10(Max_L) - NMSE_K_10_PSO(Max_L)) / NMSE_K_10(Max_L)) * 100;
    fprintf('PSO vs GD:\n');
    fprintf('  Capacity improvement: %+.2f%%\n', cap_imp_pso);
    fprintf('  NMSE improvement:     %+.2f%%\n\n', nmse_imp_pso);
end

if has_GD && has_TVIW
    cap_imp_tviw = ((Capacity_K_10_PSO_TVIW(Max_L) - Capacity_K_10(Max_L)) / Capacity_K_10(Max_L)) * 100;
    nmse_imp_tviw = ((NMSE_K_10(Max_L) - NMSE_K_10_PSO_TVIW(Max_L)) / NMSE_K_10(Max_L)) * 100;
    fprintf('PSO-TVIW vs GD:\n');
    fprintf('  Capacity improvement: %+.2f%%\n', cap_imp_tviw);
    fprintf('  NMSE improvement:     %+.2f%%\n\n', nmse_imp_tviw);
end

if has_PSO && has_TVIW
    cap_imp_tviw_pso = ((Capacity_K_10_PSO_TVIW(Max_L) - Capacity_K_10_PSO(Max_L)) / Capacity_K_10_PSO(Max_L)) * 100;
    nmse_imp_tviw_pso = ((NMSE_K_10_PSO(Max_L) - NMSE_K_10_PSO_TVIW(Max_L)) / NMSE_K_10_PSO(Max_L)) * 100;
    fprintf('PSO-TVIW vs PSO:\n');
    fprintf('  Capacity improvement: %+.2f%%\n', cap_imp_tviw_pso);
    fprintf('  NMSE improvement:     %+.2f%%\n\n', nmse_imp_tviw_pso);
end

fprintf('========================================\n');
fprintf('Note: Positive improvement = better performance\n');
fprintf('========================================\n');
