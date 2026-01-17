% % Script để so sánh kết quả giữa Gradient Descent,
    Standard PSO và TVIW - PSO clc;
clear;
close all;

% %
    Load kết quả từ thuật toán
    gốc(Gradient Descent) if exist ('Capacity_K_10.mat', 'file')
        load('Capacity_K_10.mat');
capacity_gd = Capacity_K_10;
else warning(
    'Chưa có kết quả từ thuật toán Gradient Descent. Vui lòng chạy test_SIM.m trước.');
capacity_gd = [];
end

    if exist ('NMSE_K_10.mat', 'file') load('NMSE_K_10.mat');
nmse_gd = NMSE_K_10;
else warning('Chưa có kết quả NMSE từ thuật toán Gradient Descent.');
nmse_gd = [];
end

    % %
    Load kết quả từ Standard PSO if exist ('Capacity_K_10_PSO.mat', 'file')
        load('Capacity_K_10_PSO.mat');
capacity_pso = Capacity_K_10_PSO;
else warning(
    'Chưa có kết quả từ thuật toán Standard PSO. Vui lòng chạy test_SIM_PSO.m trước.');
capacity_pso = [];
end

    if exist ('NMSE_K_10_PSO.mat', 'file') load('NMSE_K_10_PSO.mat');
nmse_pso = NMSE_K_10_PSO;
else warning('Chưa có kết quả NMSE từ thuật toán Standard PSO.');
nmse_pso = [];
end

        % % Load kết quả từ TVIW -
    PSO if exist ('Capacity_K_10_PSO_TVIW.mat', 'file')
        load('Capacity_K_10_PSO_TVIW.mat');
capacity_tviw = Capacity_K_10_PSO_TVIW;
else warning(
    'Chưa có kết quả từ thuật toán TVIW-PSO. Vui lòng chạy test_SIM_PSO_TVIW.m trước.');
capacity_tviw = [];
end

    if exist ('NMSE_K_10_PSO_TVIW.mat', 'file') load('NMSE_K_10_PSO_TVIW.mat');
nmse_tviw = NMSE_K_10_PSO_TVIW;
else warning('Chưa có kết quả NMSE từ thuật toán TVIW-PSO.');
nmse_tviw = [];
end

        % % So sánh Spectral Efficiency(Sum Rate) if ~isempty(capacity_gd) &&
    (~isempty(capacity_pso) || ~isempty(capacity_tviw))
        figure('Position', [ 100, 100, 900, 600 ]);
L_values = 1
    : max([ length(capacity_gd), length(capacity_pso), length(capacity_tviw) ]);

plot(1 : length(capacity_gd), capacity_gd, 'b-o', 'LineWidth', 2.5,
     'MarkerSize', 8, ... 'MarkerFaceColor', 'b', 'DisplayName',
     'Gradient Descent');
hold on;

if
  ~isempty(capacity_pso)
      plot(1 : length(capacity_pso), capacity_pso, 'r-s', 'LineWidth', 2.5,
           'MarkerSize', 8, ... 'MarkerFaceColor', 'r', 'DisplayName',
           'Standard PSO');
end

    if ~isempty(capacity_tviw)
        plot(1 : length(capacity_tviw), capacity_tviw, 'g-^', 'LineWidth', 2.5,
             'MarkerSize', 8, ... 'MarkerFaceColor', 'g', 'DisplayName',
             'TVIW-PSO');
end

    xlabel('Number of TX-SIM Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Spectral Efficiency Comparison: GD vs PSO vs TVIW-PSO', 'FontSize', 14,
      'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);

% Tính toán và hiển thị improvement
        fprintf('\n========== SPECTRAL EFFICIENCY COMPARISON ==========\n');
fprintf('L  |    GD    |    PSO   |  TVIW-PSO  | PSO vs GD | TVIW vs GD\n');
fprintf('------------------------------------------------------------------\n');

min_len = min([length(capacity_gd) length(capacity_pso) length(capacity_tviw)]);
if isempty (capacity_pso)
  , capacity_pso = zeros(size(capacity_gd));
end if isempty (capacity_tviw), capacity_tviw = zeros(size(capacity_gd)); end
    
    for i = 1:min_len
        imp_pso = ((capacity_pso(i) - capacity_gd(i)) / capacity_gd(i)) * 100;
imp_tviw = ((capacity_tviw(i) - capacity_gd(i)) / capacity_gd(i)) * 100;

fprintf('%2d | %8.4f | %8.4f | %8.4f   |  %+6.2f%%  |  %+6.2f%%\n', ... i,
        capacity_gd(i), capacity_pso(i), capacity_tviw(i), imp_pso, imp_tviw);
end fprintf(
    '------------------------------------------------------------------\n\n');
end

        % % So sánh NMSE if ~isempty(nmse_gd) &&
    (~isempty(nmse_pso) || ~isempty(nmse_tviw))
        figure('Position', [ 150, 150, 900, 600 ]);

semilogy(1 : length(nmse_gd), nmse_gd, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8,
         ... 'MarkerFaceColor', 'b', 'DisplayName', 'Gradient Descent');
hold on;

if
  ~isempty(nmse_pso)
      semilogy(1 : length(nmse_pso), nmse_pso, 'r-s', 'LineWidth', 2.5,
               'MarkerSize', 8, ... 'MarkerFaceColor', 'r', 'DisplayName',
               'Standard PSO');
end

    if ~isempty(nmse_tviw)
        semilogy(1 : length(nmse_tviw), nmse_tviw, 'g-^', 'LineWidth', 2.5,
                 'MarkerSize', 8, ... 'MarkerFaceColor', 'g', 'DisplayName',
                 'TVIW-PSO');
end

    xlabel('Number of TX-SIM Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE', 'FontSize', 12, 'FontWeight', 'bold');
title('NMSE Comparison: GD vs PSO vs TVIW-PSO', 'FontSize', 14, 'FontWeight',
      'bold');
legend('Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);
end

        % % Tạo bảng so sánh tổng hợp(Bar Chart) if ~isempty(capacity_gd) &&
    ~isempty(capacity_pso) &&
    ~isempty(capacity_tviw)

        figure('Position', [ 200, 200, 1000, 400 ]);
L_values = 1
    : min([ length(capacity_gd), length(capacity_pso), length(capacity_tviw) ]);

% Subplot 1 : Bar chart so sánh Capacity subplot(1, 2, 1);
    bar_data = [capacity_gd(L_values)', capacity_pso(L_values)', capacity_tviw(L_values)'];
    b = bar(L_values, bar_data);
    b(1).FaceColor = [0.2, 0.4, 0.8]; % Blue
    b(2).FaceColor = [0.8, 0.2, 0.2]; % Red
    b(3).FaceColor = [0.2, 0.8, 0.2]; % Green
    xlabel('Number of TX-SIM Layers (L)');
    ylabel('Sum Rate (bits/s/Hz)');
    title('Spectral Efficiency Comparison');
    legend('GD', 'PSO', 'TVIW-PSO', 'Location', 'best');
    grid on;
    
    % Subplot 2: Bar chart so sánh NMSE
    subplot(1,2,2);
    bar_data2 = [nmse_gd(L_values)', nmse_pso(L_values)', nmse_tviw(L_values)'];
    b2 = bar(L_values, bar_data2);
    b2(1).FaceColor = [0.2, 0.4, 0.8];
    b2(2).FaceColor = [0.8, 0.2, 0.2];
    b2(3).FaceColor = [0.2, 0.8, 0.2];
    xlabel('Number of TX-SIM Layers (L)');
    ylabel('NMSE');
    title('NMSE Comparison');
    legend('GD', 'PSO', 'TVIW-PSO', 'Location', 'best');
    grid on;
    set(gca, 'YScale', 'log');
end
