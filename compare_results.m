% % Compare Results : Gradient Descent vs PSO vs PSO -
    TVIW %
        This script loads and
            compares the results from three optimization algorithms
    : %
      1. Gradient Descent(GD) -
    from test_SIM.m % 2. Standard PSO - from test_SIM_PSO.m % 3. PSO with Time
    - Varying Inertia Weight(PSO - TVIW) -
    from test_SIM_PSO_TVIW.m

    clc;
clear;
close all;

fprintf('========================================\n');
fprintf('  Comparing Optimization Algorithms     \n');
fprintf('========================================\n\n');

% % Load Results %
        Check and load Gradient Descent results if exist ('NMSE_K_10.mat',
                                                          'file') &&
    exist('Capacity_K_10.mat', 'file') load('NMSE_K_10.mat');
load('Capacity_K_10.mat');
fprintf('[✓] Loaded Gradient Descent results\n');
has_GD = true;
else fprintf('[✗] Gradient Descent results not found. Run test_SIM.m first.\n');
has_GD = false;
end

        % Check and load Standard PSO results if exist ('NMSE_K_10_PSO.mat',
                                                        'file') &&
    exist('Capacity_K_10_PSO.mat', 'file') load('NMSE_K_10_PSO.mat');
load('Capacity_K_10_PSO.mat');
fprintf('[✓] Loaded Standard PSO results\n');
has_PSO = true;
else fprintf('[✗] Standard PSO results not found. Run test_SIM_PSO.m first.\n');
has_PSO = false;
end

            % Check and load PSO -
        TVIW results if exist ('NMSE_K_10_PSO_TVIW.mat', 'file') &&
    exist('Capacity_K_10_PSO_TVIW.mat', 'file') load('NMSE_K_10_PSO_TVIW.mat');
load('Capacity_K_10_PSO_TVIW.mat');
fprintf('[✓] Loaded PSO-TVIW results\n');
has_TVIW = true;
else fprintf(
    '[✗] PSO-TVIW results not found. Run test_SIM_PSO_TVIW.m first.\n');
has_TVIW = false;
end

    fprintf('\n');

if
  ~has_GD && ~has_PSO &&
      ~has_TVIW error(
          'No results found! Please run at least one optimization script first.');
end

    % % Determine the number of layers Max_L = 10;
if has_GD
  Max_L = length(NMSE_K_10);
elseif has_PSO Max_L = length(NMSE_K_10_PSO);
elseif has_TVIW Max_L = length(NMSE_K_10_PSO_TVIW);
end

    % % Create Comparison Plots figure('Position', [ 100, 100, 1400, 600 ]);

% % Plot 1 : Capacity Comparison subplot(1, 2, 1);
hold on;
grid on;

if has_GD
  plot(1 : Max_L, Capacity_K_10, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8,
       'MarkerFaceColor', 'b', 'DisplayName', 'Gradient Descent');
end if has_PSO plot(1 : Max_L, Capacity_K_10_PSO, 'r-s', 'LineWidth', 2.5,
                    'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName',
                    'PSO (Constant w)');
end if has_TVIW plot(1 : Max_L, Capacity_K_10_PSO_TVIW, 'g-^', 'LineWidth', 2.5,
                     'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName',
                     'PSO-TVIW');
end

    xlabel('Number of TX-SIM Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Spectral Efficiency Comparison', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 11);
hold off;

% % Plot 2 : NMSE Comparison subplot(1, 2, 2);
hold on;
grid on;

if has_GD
  semilogy(1 : Max_L, NMSE_K_10, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8,
           'MarkerFaceColor', 'b', 'DisplayName', 'Gradient Descent');
end if has_PSO semilogy(1 : Max_L, NMSE_K_10_PSO, 'r-s', 'LineWidth', 2.5,
                        'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName',
                        'PSO (Constant w)');
end if has_TVIW semilogy(1 : Max_L, NMSE_K_10_PSO_TVIW, 'g-^', 'LineWidth', 2.5,
                         'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName',
                         'PSO-TVIW');
end

    xlabel('Number of TX-SIM Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE', 'FontSize', 12, 'FontWeight', 'bold');
title('NMSE Comparison', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 11);
hold off;

% %
    Print Numerical Comparison
        fprintf('========================================\n');
fprintf('  Numerical Results Comparison          \n');
fprintf('========================================\n\n');

fprintf('%-5s | %-20s | %-20s | %-20s\n', 'L', 'Capacity (bits/s/Hz)',
        'Capacity (bits/s/Hz)', 'Capacity (bits/s/Hz)');
fprintf('%-5s | %-20s | %-20s | %-20s\n', '', 'GD', 'PSO', 'PSO-TVIW');
fprintf(
    '------+----------------------+----------------------+----------------------\n');

for
  i = 1 : Max_L fprintf('%-5d | ', i);
if has_GD
  fprintf('%-20.4f | ', Capacity_K_10(i));
else
  fprintf('%-20s | ', 'N/A');
end if has_PSO fprintf('%-20.4f | ', Capacity_K_10_PSO(i));
else fprintf('%-20s | ', 'N/A');
end if has_TVIW fprintf('%-20.4f\n', Capacity_K_10_PSO_TVIW(i));
else fprintf('%-20s\n', 'N/A');
end end

    fprintf('\n');
fprintf('%-5s | %-20s | %-20s | %-20s\n', 'L', 'NMSE', 'NMSE', 'NMSE');
fprintf('%-5s | %-20s | %-20s | %-20s\n', '', 'GD', 'PSO', 'PSO-TVIW');
fprintf(
    '------+----------------------+----------------------+----------------------\n');

for
  i = 1 : Max_L fprintf('%-5d | ', i);
if has_GD
  fprintf('%-20.6f | ', NMSE_K_10(i));
else
  fprintf('%-20s | ', 'N/A');
end if has_PSO fprintf('%-20.6f | ', NMSE_K_10_PSO(i));
else fprintf('%-20s | ', 'N/A');
end if has_TVIW fprintf('%-20.6f\n', NMSE_K_10_PSO_TVIW(i));
else fprintf('%-20s\n', 'N/A');
end end

    % %
    Performance Summary fprintf('\n========================================\n');
fprintf('  Performance Summary (L=%d)            \n', Max_L);
fprintf('========================================\n\n');

if has_GD
  fprintf('Gradient Descent:\n');
fprintf('  Capacity: %.4f bits/s/Hz\n', Capacity_K_10(Max_L));
fprintf('  NMSE:     %.6f\n\n', NMSE_K_10(Max_L));
end

    if has_PSO fprintf('Standard PSO (Constant w=0.7):\n');
fprintf('  Capacity: %.4f bits/s/Hz\n', Capacity_K_10_PSO(Max_L));
fprintf('  NMSE:     %.6f\n\n', NMSE_K_10_PSO(Max_L));
end

    if has_TVIW fprintf('PSO-TVIW (w: 0.9→0.4):\n');
fprintf('  Capacity: %.4f bits/s/Hz\n', Capacity_K_10_PSO_TVIW(Max_L));
fprintf('  NMSE:     %.6f\n\n', NMSE_K_10_PSO_TVIW(Max_L));
end

    % % Relative Performance Analysis if has_GD &&has_PSO cap_improvement_PSO =
    ((Capacity_K_10_PSO(Max_L) - Capacity_K_10(Max_L)) / Capacity_K_10(Max_L)) *
    100;
nmse_improvement_PSO =
    ((NMSE_K_10(Max_L) - NMSE_K_10_PSO(Max_L)) / NMSE_K_10(Max_L)) * 100;

fprintf('PSO vs GD:\n');
fprintf('  Capacity improvement: %+.2f%%\n', cap_improvement_PSO);
fprintf('  NMSE improvement:     %+.2f%%\n\n', nmse_improvement_PSO);
end

    if has_GD &&has_TVIW cap_improvement_TVIW =
        ((Capacity_K_10_PSO_TVIW(Max_L) - Capacity_K_10(Max_L)) /
         Capacity_K_10(Max_L)) *
        100;
nmse_improvement_TVIW =
    ((NMSE_K_10(Max_L) - NMSE_K_10_PSO_TVIW(Max_L)) / NMSE_K_10(Max_L)) * 100;

fprintf('PSO-TVIW vs GD:\n');
fprintf('  Capacity improvement: %+.2f%%\n', cap_improvement_TVIW);
fprintf('  NMSE improvement:     %+.2f%%\n\n', nmse_improvement_TVIW);
end

    if has_PSO &&has_TVIW cap_improvement_TVIW_PSO =
        ((Capacity_K_10_PSO_TVIW(Max_L) - Capacity_K_10_PSO(Max_L)) /
         Capacity_K_10_PSO(Max_L)) *
        100;
nmse_improvement_TVIW_PSO =
    ((NMSE_K_10_PSO(Max_L) - NMSE_K_10_PSO_TVIW(Max_L)) /
     NMSE_K_10_PSO(Max_L)) *
    100;

fprintf('PSO-TVIW vs Standard PSO:\n');
fprintf('  Capacity improvement: %+.2f%%\n', cap_improvement_TVIW_PSO);
fprintf('  NMSE improvement:     %+.2f%%\n\n', nmse_improvement_TVIW_PSO);
end

    fprintf('========================================\n');
fprintf('Note: Positive improvement means better performance\n');
fprintf('========================================\n');
