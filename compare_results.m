%% Script để so sánh kết quả giữa Gradient Descent và PSO
clc;
clear;
close all;

%% Load kết quả từ thuật toán gốc (Gradient Descent)
if exist('Capacity_K_10.mat', 'file')
    load('Capacity_K_10.mat');
    capacity_gd = Capacity_K_10;
else
    warning('Chưa có kết quả từ thuật toán Gradient Descent. Vui lòng chạy test_SIM.m trước.');
    capacity_gd = [];
end

if exist('NMSE_K_10.mat', 'file')
    load('NMSE_K_10.mat');
    nmse_gd = NMSE_K_10;
else
    warning('Chưa có kết quả NMSE từ thuật toán Gradient Descent.');
    nmse_gd = [];
end

%% Load kết quả từ PSO
if exist('Capacity_K_10_PSO.mat', 'file')
    load('Capacity_K_10_PSO.mat');
    capacity_pso = Capacity_K_10_PSO;
else
    warning('Chưa có kết quả từ thuật toán PSO. Vui lòng chạy test_SIM_PSO.m trước.');
    capacity_pso = [];
end

if exist('NMSE_K_10_PSO.mat', 'file')
    load('NMSE_K_10_PSO.mat');
    nmse_pso = NMSE_K_10_PSO;
else
    warning('Chưa có kết quả NMSE từ thuật toán PSO.');
    nmse_pso = [];
end

%% So sánh Spectral Efficiency (Sum Rate)
if ~isempty(capacity_gd) && ~isempty(capacity_pso)
    figure('Position', [100, 100, 800, 600]);
    L_values = 1:length(capacity_gd);
    
    plot(L_values, capacity_gd, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Gradient Descent');
    hold on;
    plot(L_values, capacity_pso, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', 'DisplayName', 'Standard PSO');
    
    xlabel('Number of TX-SIM Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Spectral Efficiency Comparison: Gradient Descent vs PSO', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    grid on;
    set(gca, 'FontSize', 11);
    
    % Tính toán và hiển thị improvement
    fprintf('\n========== SPECTRAL EFFICIENCY COMPARISON ==========\n');
    fprintf('Layer (L) | GD (bits/s/Hz) | PSO (bits/s/Hz) | Improvement (%%)\n');
    fprintf('------------------------------------------------------------\n');
    for i = 1:length(L_values)
        improvement = ((capacity_pso(i) - capacity_gd(i)) / capacity_gd(i)) * 100;
        fprintf('   %2d     |    %8.4f    |    %8.4f     |    %+7.2f%%\n', ...
            L_values(i), capacity_gd(i), capacity_pso(i), improvement);
    end
    fprintf('------------------------------------------------------------\n');
    fprintf('Average   |    %8.4f    |    %8.4f     |    %+7.2f%%\n', ...
        mean(capacity_gd), mean(capacity_pso), ...
        ((mean(capacity_pso) - mean(capacity_gd)) / mean(capacity_gd)) * 100);
    fprintf('====================================================\n\n');
end

%% So sánh NMSE
if ~isempty(nmse_gd) && ~isempty(nmse_pso)
    figure('Position', [150, 150, 800, 600]);
    L_values = 1:length(nmse_gd);
    
    semilogy(L_values, nmse_gd, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Gradient Descent');
    hold on;
    semilogy(L_values, nmse_pso, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', 'DisplayName', 'Standard PSO');
    
    xlabel('Number of TX-SIM Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('NMSE', 'FontSize', 12, 'FontWeight', 'bold');
    title('NMSE Comparison: Gradient Descent vs PSO', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    grid on;
    set(gca, 'FontSize', 11);
    
    % Tính toán và hiển thị improvement
    fprintf('\n================ NMSE COMPARISON ==================\n');
    fprintf('Layer (L) |    GD (NMSE)   |   PSO (NMSE)   | Reduction (%%)\n');
    fprintf('------------------------------------------------------------\n');
    for i = 1:length(L_values)
        reduction = ((nmse_gd(i) - nmse_pso(i)) / nmse_gd(i)) * 100;
        fprintf('   %2d     |    %10.6f    |   %10.6f   |    %+7.2f%%\n', ...
            L_values(i), nmse_gd(i), nmse_pso(i), reduction);
    end
    fprintf('------------------------------------------------------------\n');
    fprintf('Average   |    %10.6f    |   %10.6f   |    %+7.2f%%\n', ...
        mean(nmse_gd), mean(nmse_pso), ...
        ((nmse_gd(1) - nmse_pso(1)) / nmse_gd(1)) * 100);
    fprintf('====================================================\n\n');
end

%% Tạo bảng so sánh tổng hợp
if ~isempty(capacity_gd) && ~isempty(capacity_pso) && ...
   ~isempty(nmse_gd) && ~isempty(nmse_pso)
    
    figure('Position', [200, 200, 1000, 400]);
    
    % Subplot 1: Bar chart so sánh Capacity
    subplot(1,2,1);
    L_values = 1:min(length(capacity_gd), length(capacity_pso));
    bar_data = [capacity_gd(L_values)', capacity_pso(L_values)'];
    b = bar(L_values, bar_data);
    b(1).FaceColor = [0.2, 0.4, 0.8];
    b(2).FaceColor = [0.8, 0.2, 0.2];
    xlabel('Number of TX-SIM Layers (L)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Spectral Efficiency Comparison', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Gradient Descent', 'Standard PSO', 'Location', 'best');
    grid on;
    
    % Subplot 2: Bar chart so sánh NMSE
    subplot(1,2,2);
    bar_data2 = [nmse_gd(L_values)', nmse_pso(L_values)'];
    b2 = bar(L_values, bar_data2);
    b2(1).FaceColor = [0.2, 0.4, 0.8];
    b2(2).FaceColor = [0.8, 0.2, 0.2];
    xlabel('Number of TX-SIM Layers (L)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('NMSE', 'FontSize', 11, 'FontWeight', 'bold');
    title('NMSE Comparison', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Gradient Descent', 'Standard PSO', 'Location', 'best');
    grid on;
    set(gca, 'YScale', 'log');
end

%% Tổng kết
fprintf('\n==================== SUMMARY ======================\n');
if ~isempty(capacity_pso) && ~isempty(capacity_gd)
    if mean(capacity_pso) > mean(capacity_gd)
        fprintf('✓ PSO achieves BETTER spectral efficiency than Gradient Descent\n');
    else
        fprintf('✗ PSO achieves WORSE spectral efficiency than Gradient Descent\n');
    end
end

if ~isempty(nmse_pso) && ~isempty(nmse_gd)
    if mean(nmse_pso) < mean(nmse_gd)
        fprintf('✓ PSO achieves LOWER NMSE than Gradient Descent\n');
    else
        fprintf('✗ PSO achieves HIGHER NMSE than Gradient Descent\n');
    end
end
fprintf('===================================================\n\n');
