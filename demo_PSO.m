%% Demo PSO cho SIM-MIMO (Fast Version - chỉ test với L=1,2,3)
% Script này chạy nhanh hơn để demo thuật toán PSO
clc;
clearvars;
close all;

fprintf('========================================\n');
fprintf('  PSO Demo for SIM-MIMO Optimization  \n');
fprintf('========================================\n\n');

%% System Parameters
Thickness = 0.05;
Pt = 10^(20/10);
Sigma2 = 10^(-110/10);
c = 3*10^8;
f0 = 28*10^9;
lambda = c/f0;
N_max = 10;
PL = -20*log10(4*pi/lambda)-35*log10(250);
pathloss = 10^(PL/10);
M = 100;
N = 100;
d_element_spacing = lambda/2;
S = 4;
MonteCarlo = 3;  % Giảm từ 10 xuống 3 để demo nhanh
Max_L = 3;  % Chỉ test 3 layer đầu tiên
K = 10;

%% PSO Parameters
num_particles = 20;  % Giảm từ 30 xuống 20
max_iterations = 50;  % Giảm từ 100 xuống 50
w = 0.7;  % CONSTANT (Standard PSO)
c1 = 1.5;
c2 = 1.5;

fprintf('PSO Parameters:\n');
fprintf('  - Number of particles: %d\n', num_particles);
fprintf('  - Max iterations: %d\n', max_iterations);
fprintf('  - Inertia weight (w): %.2f (CONSTANT)\n', w);
fprintf('  - Cognitive param (c1): %.2f\n', c1);
fprintf('  - Social param (c2): %.2f\n\n', c2);

NMSE_demo = zeros(MonteCarlo,1);
Capacity_demo = zeros(MonteCarlo,1);
NMSE_average_demo = zeros(Max_L,1);
Capacity_average_demo = zeros(Max_L,1);

for ii = 1:Max_L
    L = ii;
    fprintf('\n======== Testing with L = %d TX-SIM layers ========\n', L);
    tic
    
    W_T = zeros(M,M);
    Corr_T = zeros(M,M);
    U_R = zeros(N,N);
    C_single_stream = zeros(S,1);
    Corr_R = zeros(N,N);
    d_layer_spacing_transmit = Thickness/L;
    d_layer_spacing_receive = Thickness/K;
    W_T_1 = zeros(M,S);
    U_R_1 = zeros(S,N);
    
    %% Calculate transmission matrices
    for mm1 = 1:M
        m_z = ceil(mm1/N_max);
        m_x = mod(mm1-1,N_max)+1;
        for mm2 = 1:M
            n_z = ceil(mm2/N_max);
            n_x = mod(mm2-1,N_max)+1;
            d_temp  = sqrt((m_x-n_x)^2 + (m_z-n_z)^2)*d_element_spacing;
            d_temp2 = sqrt(d_layer_spacing_transmit^2 + d_temp^2);
            W_T(mm2,mm1) = lambda/4/pi/d_temp2*exp(-1i*2*pi*d_temp2/lambda);
            Corr_T(mm2,mm1) = sinc(2*d_temp/lambda);
        end
    end
    
    for nn1 = 1:N
        m_z = ceil(nn1/N_max);
        m_x = mod(nn1-1,N_max)+1;
        for nn2 = 1:N
            n_z = ceil(nn2/N_max);
            n_x = mod(nn2-1,N_max)+1;
            d_temp  = sqrt((m_x-n_x)^2 + (m_z-n_z)^2)*d_element_spacing;
            d_temp2 = sqrt(d_layer_spacing_receive^2 + d_temp^2);
            U_R(nn2,nn1) = lambda/4/pi/d_temp2*exp(-1i*2*pi*d_temp2/lambda);
            Corr_R(nn2,nn1) = sinc(2*d_temp/lambda);
        end
    end
    
    for mm = 1:M
        m_z = ceil(mm/N_max);
        m_x = mod(mm-1,N_max)+1;
        for nn = 1:S
            d_transmit = sqrt(d_layer_spacing_transmit^2 + ...
                ((m_x-(1+N_max)/2)*d_element_spacing)^2 + ...
                ((m_z-(1+N_max)/2)*d_element_spacing - (nn-(1+S)/2)*lambda/2)^2);
            W_T_1(mm,nn) = lambda/4/pi/d_transmit*exp(-1i*2*pi*d_transmit/lambda);
        end
    end
    
    for mm = 1:N
        m_z = ceil(mm/N_max);
        m_x = mod(mm-1,N_max)+1;
        for nn = 1:S
            d_receive = sqrt(d_layer_spacing_receive^2 +...
                ((m_x-(1+N_max)/2)*d_element_spacing)^2 +...
                ((m_z-(1+N_max)/2)*d_element_spacing - (nn-(1+S)/2)*lambda/2)^2);
            U_R_1(nn,mm) = lambda/4/pi/d_receive*exp(-1i*2*pi*d_receive/lambda);
        end
    end
    
    rng(1)
    for jj = 1:MonteCarlo
        fprintf('\n--- Monte Carlo iteration %d/%d ---\n', jj, MonteCarlo);
        
        G_independent = sqrt(1/2)*(randn(N,M)+1i*randn(N,M));
        G = sqrt(pathloss)*(Corr_R)^(1/2)*G_independent*(Corr_T)^(1/2);
        [G_left, G_svd, G_right] = svd(G);
        H_true = G_svd(1:S,1:S);
        H_true_vec = H_true(:);  % Vectorize matrix (no toolbox needed)
        Norm_H = norm(H_true_vec)^2;
        h_diag = diag(H_true);
        
        if S == 1
            PA_WF = Pt;
        else
            [PA_WF] = WF(Pt, Sigma2, h_diag);
        end
        
        %% PSO Initialization
        dim = M*L + N*K;
        fprintf('Particle dimension: %d (TX: %d, RX: %d)\n', dim, M*L, N*K);
        
        particles = -pi + 2*pi*rand(num_particles, dim);
        velocities = -pi + 2*pi*rand(num_particles, dim);
        pBest = particles;
        pBest_fitness = -inf(num_particles, 1);
        gBest = zeros(1, dim);
        gBest_fitness = -inf;
        
        % Initial evaluation
        fprintf('Evaluating initial population...\n');
        for p = 1:num_particles
            fitness = evaluate_fitness_PSO(particles(p,:), M, L, N, K, W_T, W_T_1, U_R, U_R_1, G, PA_WF, Sigma2, H_true, H_true_vec, Norm_H);
            if fitness > pBest_fitness(p)
                pBest_fitness(p) = fitness;
                pBest(p,:) = particles(p,:);
            end
            if fitness > gBest_fitness
                gBest_fitness = fitness;
                gBest = particles(p,:);
            end
        end
        fprintf('Initial best fitness: %.4f bits/s/Hz\n', gBest_fitness);
        
        %% PSO Main Loop
        for iter = 1:max_iterations
            for p = 1:num_particles
                % Standard PSO update: v = w*v + c1*rand*(pBest - x) + c2*rand*(gBest - x)
                velocities(p,:) = w * velocities(p,:) + ...
                                  c1 * rand(1,dim) .* (pBest(p,:) - particles(p,:)) + ...
                                  c2 * rand(1,dim) .* (gBest - particles(p,:));
                
                max_velocity = pi/2;
                velocities(p,:) = max(min(velocities(p,:), max_velocity), -max_velocity);
                
                particles(p,:) = particles(p,:) + velocities(p,:);
                particles(p,:) = wrapToPi(particles(p,:));
                
                fitness = evaluate_fitness_PSO(particles(p,:), M, L, N, K, W_T, W_T_1, U_R, U_R_1, G, PA_WF, Sigma2, H_true, H_true_vec, Norm_H);
                
                if fitness > pBest_fitness(p)
                    pBest_fitness(p) = fitness;
                    pBest(p,:) = particles(p,:);
                end
                
                if fitness > gBest_fitness
                    gBest_fitness = fitness;
                    gBest = particles(p,:);
                end
            end
            
            if mod(iter, 10) == 0 || iter == max_iterations
                fprintf('  Iteration %3d: Best Capacity = %.4f bits/s/Hz\n', iter, gBest_fitness);
            end
        end
        
        %% Extract and evaluate final solution
        phase_transmit_angles = reshape(gBest(1:M*L), M, L);
        phase_receive_angles = reshape(gBest(M*L+1:end), N, K);
        
        phase_transmit = exp(1i*phase_transmit_angles);
        phase_receive = exp(1i*phase_receive_angles);
        
        P = diag(phase_transmit(:,1))*W_T_1;
        for l=1:L-1
            P = diag(phase_transmit(:,l+1))*W_T*P;
        end
        
        Q = U_R_1*diag(phase_receive(:,1));
        for k = 1:K-1
            Q = Q*U_R*diag(phase_receive(:,k+1));
        end
        
        H_SIM = Q*G*P;
        H_SIM_vec = H_SIM(:);  % Vectorize matrix (no toolbox needed)
        Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec;
        
        Error_new = norm(Factor*H_SIM-H_true)^2/Norm_H;
        NMSE_demo(jj) = Error_new;
        
        for pp = 1:S
            C_single_stream(pp) = log2(1+PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2/ ...
                (Sigma2+(abs(Factor*H_SIM(pp,:)).^2*PA_WF-PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2)));
        end
        Capacity_demo(jj) = sum(C_single_stream);
        
        fprintf('Final: Capacity = %.4f bits/s/Hz, NMSE = %.6f\n', Capacity_demo(jj), NMSE_demo(jj));
    end
    
    NMSE_average_demo(ii) = mean(NMSE_demo);
    Capacity_average_demo(ii) = mean(Capacity_demo);
    
    elapsed = toc;
    fprintf('\nL=%d completed in %.2f seconds\n', L, elapsed);
    fprintf('Average Capacity: %.4f bits/s/Hz\n', Capacity_average_demo(ii));
    fprintf('Average NMSE: %.6f\n', NMSE_average_demo(ii));
end

%% Plot results
figure('Position', [100, 100, 1200, 400]);

subplot(1,2,1);
plot(1:Max_L, Capacity_average_demo, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Number of TX-SIM Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('PSO: Spectral Efficiency vs L', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

subplot(1,2,2);
semilogy(1:Max_L, NMSE_average_demo, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
xlabel('Number of TX-SIM Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE', 'FontSize', 12, 'FontWeight', 'bold');
title('PSO: NMSE vs L', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

%% Summary
fprintf('\n\n========================================\n');
fprintf('         DEMO COMPLETED!               \n');
fprintf('========================================\n');
fprintf('Results Summary:\n');
for i = 1:Max_L
    fprintf('  L=%d: Capacity=%.4f bits/s/Hz, NMSE=%.6f\n', i, Capacity_average_demo(i), NMSE_average_demo(i));
end
fprintf('========================================\n\n');

fprintf('Để chạy full experiment (L=1 to 10, MC=10):\n');
fprintf('  >> test_SIM_PSO\n\n');

fprintf('Để so sánh với thuật toán gốc:\n');
fprintf('  >> compare_results\n\n');

%% Fitness Function
function fitness = evaluate_fitness_PSO(particle, M, L, N, K, W_T, W_T_1, U_R, U_R_1, G, PA_WF, Sigma2, H_true, H_true_vec, Norm_H)
    phase_transmit_angles = reshape(particle(1:M*L), M, L);
    phase_receive_angles = reshape(particle(M*L+1:end), N, K);
    
    phase_transmit = exp(1i*phase_transmit_angles);
    phase_receive = exp(1i*phase_receive_angles);
    
    P = diag(phase_transmit(:,1))*W_T_1;
    for l=1:L-1
        P = diag(phase_transmit(:,l+1))*W_T*P;
    end
    
    Q = U_R_1*diag(phase_receive(:,1));
    for k = 1:K-1
        Q = Q*U_R*diag(phase_receive(:,k+1));
    end
    
    H_SIM = Q*G*P;
    H_SIM_vec = H_SIM(:);
    S = size(H_SIM, 1);
    
    % CORRECT: Calculate Factor to match H_true (like Gradient Descent)
    Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec;
    
    if isnan(Factor) || isinf(Factor)
        fitness = -1e10;
        return;
    end
    
    sum_rate = 0;
    for pp = 1:S
        SINR = PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2 / ...
               (Sigma2 + (abs(Factor*H_SIM(pp,:)).^2*PA_WF - PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2));
        
        if SINR > 0 && ~isnan(SINR) && ~isinf(SINR)
            sum_rate = sum_rate + log2(1 + SINR);
        end
    end
    
    fitness = sum_rate;
    
    if isnan(fitness) || isinf(fitness)
        fitness = -1e10;
    end
end
