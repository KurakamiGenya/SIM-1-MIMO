clc;
clearvars;
close all;

%% Start execution timer
algorithm_start_time = tic;

%% Algorithm Header
fprintf('\n');
fprintf('╔════════════════════════════════════════╗\n');
fprintf('║        PSO-TVIW ALGORITHM              ║\n');
fprintf('╚════════════════════════════════════════╝\n');
fprintf('Status: RUNNING...\n');
fprintf('Method: PSO with Time-Varying Inertia Weight\n');
fprintf('Parameters: w varies from 0.9 to 0.4\n');
fprintf('----------------------------------------\n\n');

%% System Parameters
Thickness = 0.05; %% Thickness of TX-SIM and RX-SIM
Pt = 10^(20/10); %% Transmit power
Sigma2 = 10^(-110/10); %% Average noise power at the receiver
c = 3*10^8; %% Speed of light
f0 = 28*10^9; %% Radio frequency
lambda = c/f0; %% Wavelength
N_max = 10; %% Number of meta-atoms on each row
PL = -20*log10(4*pi/lambda)-35*log10(250); %% Pathloss in dB
pathloss = 10^(PL/10); %% Pathloss
M = 100; %% Number of meta-atoms on each layer of TX-SIM
N = 100; %% Number of meta-atoms on each layer of RX-SIM
d_element_spacing = lambda/2; %% Element spacing
S = 4; %% Number of data streams
MonteCarlo = 10; %% Number of independent experiments
Max_L = 10; %% The maximum number of metasurface layers in TX-SIM
K = 10;  %% The number of metasurface layers in RX-SIM

%% PSO Parameters with Time-Varying Inertia Weight (TVIW)
num_particles = 30; %% Number of particles in swarm
max_iterations = 100; %% Maximum number of PSO iterations
w_max = 0.9; %% Maximum inertia weight (TIME-VARYING)
w_min = 0.4; %% Minimum inertia weight (TIME-VARYING)
c1 = 1.5; %% Cognitive parameter
c2 = 1.5; %% Social parameter

NMSE = zeros(MonteCarlo,1);
Capacity = zeros(MonteCarlo,1);
NMSE_average = zeros(Max_L,1);
Capacity_average = zeros(Max_L,1);
for ii = 1:Max_L
    L = ii; %% The number of metasurface layers in TX-SIM
    fprintf('\n======== L = %d ========\n', L);
    tic
    
    W_T = zeros(M,M);
    Corr_T = zeros(M,M);
    U_R = zeros(N,N);
    C_single_stream = zeros(S,1);
    Corr_R = zeros(N,N);
    d_layer_spacing_transmit = Thickness/L; %% Adjacent layer spacing in TX-SIM
    d_layer_spacing_receive = Thickness/K; %% Adjacent layer spacing in RX-SIM
    W_T_1 = zeros(M,S);
    U_R_1 = zeros(S,N);
    
    %% Calculate inter-layer transmission coefficient matrix W_T and channel correlation matrix Corr_T associated with TX-SIM
    for mm1 = 1:M
        m_z = ceil(mm1/N_max);
        m_x = mod(mm1-1,N_max)+1;
        for mm2 = 1:M
            n_z = ceil(mm2/N_max);
            n_x = mod(mm2-1,N_max)+1;
            d_temp  = sqrt(  (m_x-n_x)^2 +  (m_z-n_z) ^2 )*d_element_spacing;
            d_temp2 = sqrt(d_layer_spacing_transmit^2 + d_temp^2);
            W_T(mm2,mm1) = lambda/4/pi/d_temp2*exp(-1i*2*pi*d_temp2/lambda);
            Corr_T(mm2,mm1) = sinc(2*d_temp/lambda);
        end
    end
    
    %% Calculate inter-layer transmission coefficient matrix U_R and channel correlation matrix Corr_R associated with RX-SIM
    for nn1 = 1:N
        m_z = ceil(nn1/N_max);
        m_x = mod(nn1-1,N_max)+1;
        for nn2 = 1:N
            n_z = ceil(nn2/N_max);
            n_x = mod(nn2-1,N_max)+1;
            d_temp  = sqrt( (m_x-n_x)^2 + (m_z-n_z)^2 )*d_element_spacing;
            d_temp2 = sqrt(d_layer_spacing_receive^2 + d_temp^2);
            U_R(nn2,nn1) = lambda/4/pi/d_temp2*exp(-1i*2*pi*d_temp2/lambda);
            Corr_R(nn2,nn1) = sinc(2*d_temp/lambda);
        end
    end
    
    %% The channel from transmitter to the first layer of TX-SIM
    for mm = 1:M
        m_z = ceil(mm/N_max);
        m_x = mod(mm-1,N_max)+1;
        for nn = 1:S
            d_transmit = sqrt(d_layer_spacing_transmit^2 + ...
                ( (m_x-(1+N_max)/2)*d_element_spacing )^2 + ...
                ( (m_z-(1+N_max)/2)*d_element_spacing - (nn-(1+S)/2)*lambda/2 )^2 );
            W_T_1(mm,nn) = lambda/4/pi/d_transmit*exp(-1i*2*pi*d_transmit/lambda);
        end
    end
    
    %% The channel from the last layer of RX-SIM to the receiver
    for mm = 1:N
        m_z = ceil(mm/N_max);
        m_x = mod(mm-1,N_max)+1;
        for nn = 1:S
            d_receive = sqrt(d_layer_spacing_receive^2 +...
                ( (m_x-(1+N_max)/2)*d_element_spacing  )^2 +...
                ( (m_z-(1+N_max)/2)*d_element_spacing - (nn-(1+S)/2)*lambda/2 )^2 );
            U_R_1(nn,mm) = lambda/4/pi/d_receive*exp(-1i*2*pi*d_receive/lambda);
        end
    end
    
    rng(1)
    for jj = 1:MonteCarlo
        G_independent = sqrt(1/2)*(randn(N,M)+1i*randn(N,M));
        G = sqrt(pathloss)*(Corr_R)^(1/2)*G_independent*(Corr_T)^(1/2); %% HMIMO channel
        [G_left, G_svd, G_right] = svd(G); %% SVD of HMIMO channel
        H_true = G_svd(1:S,1:S); %% Target channel
        H_true_vec = H_true(:);  % Vectorize matrix (no toolbox needed)
        Norm_H = norm(H_true_vec)^2;
        h_diag = diag(H_true);
        
        %% Power allocation using water-filling algorithm
        if S == 1
            PA_WF = Pt;
        else
            [ PA_WF ] = WF( Pt, Sigma2, h_diag );
        end
        
        %% PSO Initialization
        % Dimension: M*L (TX phases) + N*K (RX phases)
        dim = M*L + N*K;
        
        % Initialize particle positions (phase angles in [-pi, pi])
        particles = -pi + 2*pi*rand(num_particles, dim);
        
        % Initialize velocities
        velocities = -pi + 2*pi*rand(num_particles, dim);
        
        % Initialize personal best positions and fitness
        pBest = particles;
        pBest_fitness = -inf(num_particles, 1);
        
        % Initialize global best position and fitness
        gBest = zeros(1, dim);
        gBest_fitness = -inf;
        
        % Evaluate initial fitness for all particles
        for p = 1:num_particles
            fitness = evaluate_fitness_PSO(particles(p,:), M, L, N, K, W_T, W_T_1, U_R, U_R_1, G, PA_WF, Sigma2, H_true, H_true_vec, Norm_H);
            
            % Update personal best
            if fitness > pBest_fitness(p)
                pBest_fitness(p) = fitness;
                pBest(p,:) = particles(p,:);
            end
            
            % Update global best
            if fitness > gBest_fitness
                gBest_fitness = fitness;
                gBest = particles(p,:);
            end
        end
        
        %% PSO Main Loop with Time-Varying Inertia Weight
        % Initialize convergence tracking for L=Max_L
        if L == Max_L && jj == MonteCarlo
            convergence_history = zeros(max_iterations, 1);
            convergence_time = zeros(max_iterations, 1);  % Track elapsed time
            convergence_start_time = tic;  % Start timer for convergence tracking
        end
        
        for iter = 1:max_iterations
            % Calculate time-varying inertia weight (linearly decreasing)
            w = w_max - (w_max - w_min) * iter / max_iterations;
            
            for p = 1:num_particles
                % Update velocity using PSO formula with time-varying w
                % v = w(iter)*v + c1*rand*(pBest - x) + c2*rand*(gBest - x)
                velocities(p,:) = w * velocities(p,:) + ...
                                  c1 * rand(1,dim) .* (pBest(p,:) - particles(p,:)) + ...
                                  c2 * rand(1,dim) .* (gBest - particles(p,:));
                
                % Limit velocity to prevent explosion
                max_velocity = pi/2;
                velocities(p,:) = max(min(velocities(p,:), max_velocity), -max_velocity);
                
                % Update position
                particles(p,:) = particles(p,:) + velocities(p,:);
                
                % Wrap phase angles to [-pi, pi]
                particles(p,:) = wrapToPi(particles(p,:));
                
                % Evaluate fitness
                fitness = evaluate_fitness_PSO(particles(p,:), M, L, N, K, W_T, W_T_1, U_R, U_R_1, G, PA_WF, Sigma2, H_true, H_true_vec, Norm_H);
                
                % Update personal best
                if fitness > pBest_fitness(p)
                    pBest_fitness(p) = fitness;
                    pBest(p,:) = particles(p,:);
                end
                
                % Update global best
                if fitness > gBest_fitness
                    gBest_fitness = fitness;
                    gBest = particles(p,:);
                end
            end
            
            % Track convergence for L=Max_L
            if L == Max_L && jj == MonteCarlo
                convergence_history(iter) = gBest_fitness;
                convergence_time(iter) = toc(convergence_start_time);  % Time since PSO loop started
            end
            
            % Display progress
            if mod(iter, 20) == 0
                fprintf('L=%d, MC=%d, Iter=%d, w=%.3f, Best Capacity=%.4f bits/s/Hz\n', L, jj, iter, w, gBest_fitness);
            end
        end
        
        %% Extract best solution
        phase_transmit_angles = reshape(gBest(1:M*L), M, L);
        phase_receive_angles = reshape(gBest(M*L+1:end), N, K);
        
        phase_transmit = exp(1i*phase_transmit_angles);
        phase_receive = exp(1i*phase_receive_angles);
        
        %% Calculate final TX-SIM response
        P = diag(phase_transmit(:,1))*W_T_1;
        for l=1:L-1
            P = diag(phase_transmit(:,l+1))*W_T*P;
        end
        
        %% Calculate final RX-SIM response
        Q = U_R_1*diag(phase_receive(:,1));
        for k = 1:K-1
            Q = Q*U_R*diag(phase_receive(:,k+1));
        end
        
        H_SIM = Q*G*P;
        H_SIM_vec = H_SIM(:);  % Vectorize matrix (no toolbox needed)
        Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec;
        
        %% Calculate NMSE
        Error_new = norm(Factor*H_SIM-H_true)^2/Norm_H;
        NMSE(jj) = Error_new;
        
        %% Calculate Capacity
        for pp = 1:S
            C_single_stream(pp) = log2(1+PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2/ ...
                (Sigma2+(abs(Factor*H_SIM(pp,:)).^2*PA_WF-PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2)));
        end
        Capacity(jj) = sum(C_single_stream);
        
        fprintf('L=%d, MC=%d, Final: Capacity=%.4f, NMSE=%.6f\n', L, jj, Capacity(jj), NMSE(jj));
    end
    
    NMSE_average(ii) = mean(NMSE);
    Capacity_average(ii) = mean(Capacity);
    toc
end

%% Save convergence history for plotting
if exist('convergence_history', 'var')
    save('PSO_TVIW_convergence.mat', 'convergence_history', 'convergence_time');
    fprintf('\nConvergence history saved to PSO_TVIW_convergence.mat\n');
end

%% Save execution time
algorithm_execution_time = toc(algorithm_start_time);
save('PSO_TVIW_execution_time.mat', 'algorithm_execution_time');
fprintf('Execution time: %.2f seconds (%.2f minutes)\n', algorithm_execution_time, algorithm_execution_time/60);
fprintf('Execution time saved to PSO_TVIW_execution_time.mat\n');
%% Plot Results
figure;
plot(1:Max_L, NMSE_average, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of TX-SIM Layers (L)');
ylabel('NMSE');
title('NMSE vs Number of TX-SIM Layers (PSO-TVIW)');
grid on;
figure;
plot(1:Max_L, Capacity_average, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of TX-SIM Layers (L)');
ylabel('Sum Rate (bits/s/Hz)');
title('Spectral Efficiency vs Number of TX-SIM Layers (PSO-TVIW)');
grid on;
%% Save results
NMSE_K_10_PSO_TVIW = NMSE_average;
save('NMSE_K_10_PSO_TVIW.mat', 'NMSE_K_10_PSO_TVIW');
Capacity_K_10_PSO_TVIW = Capacity_average;
save('Capacity_K_10_PSO_TVIW.mat', 'Capacity_K_10_PSO_TVIW');
fprintf('\n========================================\n');
fprintf('Results saved to:\n');
fprintf('  - NMSE_K_10_PSO_TVIW.mat\n');
fprintf('  - Capacity_K_10_PSO_TVIW.mat\n');
fprintf('========================================\n');
%% Fitness Function for PSO
function fitness = evaluate_fitness_PSO(particle, M, L, N, K, W_T, W_T_1, U_R, U_R_1, G, PA_WF, Sigma2, H_true, H_true_vec, Norm_H)
    % Extract phase angles from particle
    phase_transmit_angles = reshape(particle(1:M*L), M, L);
    phase_receive_angles = reshape(particle(M*L+1:end), N, K);
    
    % Convert to complex phase shifts
    phase_transmit = exp(1i*phase_transmit_angles);
    phase_receive = exp(1i*phase_receive_angles);
    
    % Calculate TX-SIM response
    P = diag(phase_transmit(:,1))*W_T_1;
    for l=1:L-1
        P = diag(phase_transmit(:,l+1))*W_T*P;
    end
    
    % Calculate RX-SIM response
    Q = U_R_1*diag(phase_receive(:,1));
    for k = 1:K-1
        Q = Q*U_R*diag(phase_receive(:,k+1));
    end
    
    % End-to-end channel
    H_SIM = Q*G*P;
    H_SIM_vec = H_SIM(:);
    
    % Get number of streams
    S = size(H_SIM, 1);
    
    % Calculate compensation factor (CORRECT METHOD - match to H_true)
    Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec;
    
    % Check for numerical issues
    if isnan(Factor) || isinf(Factor)
        fitness = -1e10;
        return;
    end
    
    % Calculate sum rate (Spectral Efficiency) as fitness
    sum_rate = 0;
    for pp = 1:S
        SINR = PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2 / ...
               (Sigma2 + (abs(Factor*H_SIM(pp,:)).^2*PA_WF - PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2));
        
        % Check for numerical issues
        if SINR > 0 && ~isnan(SINR) && ~isinf(SINR)
            sum_rate = sum_rate + log2(1 + SINR);
        end
    end
    
    % Return sum rate as fitness (we want to maximize this)
    fitness = sum_rate;
    
    % Penalty for invalid solutions
    if isnan(fitness) || isinf(fitness)
        fitness = -1e10;
    end
end
