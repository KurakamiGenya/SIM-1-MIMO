%% Master Script: Run All Optimization Algorithms
% This script provides instructions and timing for running all three algorithms.
% Due to MATLAB workspace limitations with 'clearvars' in each script,
% you need to run each algorithm separately and then generate plots.
%
% RECOMMENDED USAGE:
% Option 1: Run this script - it will pause between algorithms
% Option 2: Run each script manually in sequence

clc;
clear;
close all;

fprintf('========================================\n');
fprintf('  SIM-MIMO Algorithm Comparison Suite  \n');
fprintf('========================================\n\n');

fprintf('This script will run all three optimization algorithms:\n');
fprintf('  1. Gradient Descent (test_SIM.m)\n');
fprintf('  2. Standard PSO (test_SIM_PSO.m)\n');
fprintf('  3. PSO-TVIW (test_SIM_PSO_TVIW.m)\n\n');

fprintf('After completion, comparison plots will be generated.\n');
fprintf('========================================\n\n');



%% Algorithm 1: Gradient Descent
test_SIM;

fprintf('\n----------------------------------------\n');
fprintf('✓ Gradient Descent COMPLETED!\n');
fprintf('Results saved to: NMSE_K_10.mat, Capacity_K_10.mat\n');
fprintf('========================================\n\n');

%% Algorithm 2: Standard PSO
test_SIM_PSO;

fprintf('\n----------------------------------------\n');
fprintf('✓ Standard PSO COMPLETED!\n');
fprintf('Results saved to: NMSE_K_10_PSO.mat, Capacity_K_10_PSO.mat\n');
fprintf('========================================\n\n');

%% Algorithm 3: PSO-TVIW
test_SIM_PSO_TVIW;

fprintf('\n----------------------------------------\n');
fprintf('✓ PSO-TVIW COMPLETED!\n');
fprintf('Results saved to: NMSE_K_10_PSO_TVIW.mat, Capacity_K_10_PSO_TVIW.mat\n');
fprintf('========================================\n\n');

%% All algorithms completed
fprintf('╔════════════════════════════════════════╗\n');
fprintf('║   ALL ALGORITHMS COMPLETED!           ║\n');
fprintf('╚════════════════════════════════════════╝\n');
fprintf('\n');

%% Generate convergence comparison plot
fprintf('Generating convergence comparison plot...\n');
plot_convergence;
fprintf('✓ Convergence plot completed!\n\n');

%% Generate performance vs time plot
fprintf('Generating performance vs time plot...\n');
plot_performance_vs_time;
fprintf('✓ Performance vs time plot completed!\n\n');

fprintf('╔════════════════════════════════════════╗\n');
fprintf('║   ALL PLOTS GENERATED SUCCESSFULLY!   ║\n');
fprintf('╚════════════════════════════════════════╝\n');
fprintf('\n');

%% Load execution times from individual algorithm files
try
    load('GD_execution_time.mat');
    gd_time = algorithm_execution_time;
    has_gd_time = true;
catch
    gd_time = 0;
    has_gd_time = false;
end

try
    load('PSO_execution_time.mat');
    pso_time = algorithm_execution_time;
    has_pso_time = true;
catch
    pso_time = 0;
    has_pso_time = false;
end

try
    load('PSO_TVIW_execution_time.mat');
    tviw_time = algorithm_execution_time;
    has_tviw_time = true;
catch
    tviw_time = 0;
    has_tviw_time = false;
end

%% Save combined execution times
if has_gd_time && has_pso_time && has_tviw_time
    execution_times = struct();
    execution_times.GradientDescent = gd_time;
    execution_times.StandardPSO = pso_time;
    execution_times.PSO_TVIW = tviw_time;
    execution_times.TotalTime = gd_time + pso_time + tviw_time;
    
    save('algorithm_execution_times.mat', 'execution_times');
    fprintf('Combined execution times saved to: algorithm_execution_times.mat\n\n');
    
    %% Display execution time summary
    fprintf('========================================\n');
    fprintf('  Execution Time Summary                \n');
    fprintf('========================================\n');
    fprintf('Gradient Descent: %.2f seconds (%.2f minutes)\n', gd_time, gd_time/60);
    fprintf('Standard PSO:     %.2f seconds (%.2f minutes)\n', pso_time, pso_time/60);
    fprintf('PSO-TVIW:         %.2f seconds (%.2f minutes)\n', tviw_time, tviw_time/60);
    fprintf('----------------------------------------\n');
    fprintf('Total Time:       %.2f seconds (%.2f minutes)\n', execution_times.TotalTime, execution_times.TotalTime/60);
    fprintf('========================================\n\n');
else
    fprintf('⚠ Warning: Some execution time files are missing.\n\n');
end

%% Verify result files
fprintf('Verifying result files...\n');
files_ok = true;

required_files = {
    'NMSE_K_10.mat', 'Capacity_K_10.mat', ...
    'NMSE_K_10_PSO.mat', 'Capacity_K_10_PSO.mat', ...
    'NMSE_K_10_PSO_TVIW.mat', 'Capacity_K_10_PSO_TVIW.mat'
};

for i = 1:length(required_files)
    if isfile(required_files{i})
        fprintf('  ✓ %s\n', required_files{i});
    else
        fprintf('  ✗ %s (MISSING)\n', required_files{i});
        files_ok = false;
    end
end

fprintf('\n');

%% Generate comparison plots
if files_ok
    fprintf('========================================\n');
    fprintf('  Generating Comparison Plots...       \n');
    fprintf('========================================\n\n');
    
    Plot_SIM;  % Final results comparison (Sum Rate vs L)
    plot_convergence;  % Convergence comparison (Sum Rate vs Iterations)
    
    fprintf('\n✓ All comparison plots generated successfully!\n');
else
    fprintf('⚠ Warning: Some result files are missing.\n');
    fprintf('   Please check for errors in algorithm execution.\n');
end

fprintf('\n========================================\n');
fprintf('  Workflow Completed!                  \n');
fprintf('========================================\n');
