%% Simulation parameters
clear;
close all

rng(0); % Set a fixed seed for reproducibility

sim_times = 1e5; % number of simulation trials

swarm_size = 2; % number of IRS UAVs in the swarm (set to 1 for original scenario)
N = 5; % number of reflecting elements per IRS

snr_dB = 0:2:20; % average transmit SNR in dB
R_th = 1; % SE threshold b/s/Hz
snr_th = 2^R_th - 1;

PLE = 2.7; % path-loss exponent

%% Generate distances and channel parameters for each IRS
for irs = 1:swarm_size
    d_Sr(irs) = 1+rand; % random distance S->RIS
    d_rD(irs) = 1+rand; % random distance RIS->D
    
    % Nakagami-m parameters
    m_Sr(irs) = 2.5 + rand; % random shape, S->RIS
    m_rD(irs) = 2.5 + rand; % random shape, RIS->D
    Omega_Sr(irs) = d_Sr(irs)^(-PLE); % random spread, S->RIS
    Omega_rD(irs) = d_rD(irs)^(-PLE); % random spread, RIS->D

    % Inverse Gamma (IG) parameters
    alpha_Sr(irs) = 3.0+rand; % random shape, S->RIS
    alpha_rD(irs) = 3.0+rand; % random shape, RIS->D
    beta_Sr(irs) = 1; % random spread, S->RIS
    beta_rD(irs) = 1; % random spread, RIS->D
end

Z_sim = 0;
kappa = 1; % for RIS

%% Channel modeling
waitbar_h = waitbar(0, 'Generating channel models...');
for irs = 1:swarm_size
    % Nakagami-m fading channel
    G_Sr{irs} = random('Naka', m_Sr(irs), Omega_Sr(irs), [N, sim_times]);
    G_rD{irs} = random('Naka', m_rD(irs), Omega_rD(irs), [N, sim_times]);

    % Inverse Gamma shadowing
    L_Sr{irs} = 1./random('Gamma', alpha_Sr(irs), 1/beta_Sr(irs), [N, sim_times]);
    L_rD{irs} = 1./random('Gamma', alpha_rD(irs), 1/beta_rD(irs), [N, sim_times]);

    Gr{irs} = G_Sr{irs}.*G_rD{irs}; % e2e fading w.r.t. one element
    Lr{irs} = L_Sr{irs}.*L_rD{irs}; % e2e shadowing w.r.t. one element

    W_sim{irs} = Gr{irs}.*Lr{irs}; % e2e channel w.r.t. one element
    Z_sim = Z_sim + sum(W_sim{irs},1); % e2e channel w.r.t. the whole swarm
    
    waitbar(irs/swarm_size, waitbar_h);
end
close(waitbar_h);

Z2_sim = Z_sim.^2; % e2e squared magnitude

%% Optimal phase shift
waitbar_h = waitbar(0, 'Calculating optimal phase shift...');
for irs = 1:swarm_size
    % phase of channels
    phase_Sr{irs} = 2*pi*rand(N, sim_times); % domain [0,2pi)
    phase_rD{irs} = 2*pi*rand(N, sim_times); % domain [0,2pi)

    % Channel modeling
    G_Sr_complex_fading{irs} = G_Sr{irs} .* exp(1i*phase_Sr{irs});
    G_rD_complex_fading{irs} = G_rD{irs} .* exp(1i*phase_rD{irs});
    
    waitbar(irs/swarm_size, waitbar_h);
end
close(waitbar_h);

% Phase-shift configuration
Z_sim_optimal_phase_shift = zeros(1,sim_times);

waitbar_h = waitbar(0, 'Configuring phase shift...');
for ss = 1:sim_times % loop over simulation trials
    for irs = 1:swarm_size
        for ll = 1:N % loop over each element of the IRS
            % unknown domain phase-shift
            phase_shift_element_temp = - phase_Sr{irs}(ll,ss) - phase_rD{irs}(ll,ss);
            
            % convert to domain of [0, 2pi)
            phase_shift_element = wrapTo2Pi(phase_shift_element_temp);
            
            Gr_optimal_phase_shift = abs(G_Sr_complex_fading{irs}(ll,ss) * ...
                exp(1i*phase_shift_element) * G_rD_complex_fading{irs}(ll,ss));
            W_optimal_phase_shift{irs}(ll,ss) = Gr_optimal_phase_shift*Lr{irs}(ll,ss);
        end
    end
    Z_sim_optimal_phase_shift(ss) = sum(cellfun(@(x) sum(x(:,ss)), W_optimal_phase_shift));
    
    if mod(ss, sim_times/100) == 0
        waitbar(ss/sim_times, waitbar_h);
    end
end
close(waitbar_h);

Z2_sim_optimal_phase_shift = Z_sim_optimal_phase_shift.^2;

%% Analysis (adapted for multiple IRS)
L_W_fit = @(s) 1;
for irs = 1:swarm_size
    % STEP-1: MOMENT MATCHING Gr -> Gamma
    Upsilon_G = m_Sr(irs)*m_rD(irs)/Omega_Sr(irs)/Omega_rD(irs);
    Omega_G = gamma(m_Sr(irs) + 1/2) * gamma(m_rD(irs)+1/2) / gamma(m_Sr(irs)) / gamma(m_rD(irs)) * Upsilon_G^(-1/2);
    m_G = Omega_G^2/ (Omega_Sr(irs)*Omega_rD(irs) - Omega_G^2);

    % STEP-2: MOMENT MATCHING 1./sqrt(Lr) -> Gamma
    Omega_L = gamma(alpha_Sr(irs)+1/2)*gamma(alpha_rD(irs)+1/2) / sqrt(beta_Sr(irs)*beta_rD(irs)) / gamma(alpha_Sr(irs)) / gamma(alpha_rD(irs));
    m_L = Omega_L^2 / (alpha_Sr(irs)*alpha_rD(irs)/beta_Sr(irs)/beta_rD(irs) - Omega_L^2);

    % STEP-3: APPROXIMATE THE PDF of \tilde{R}_{r} USING MG
    K = 5; % number of terms in G-L quadrature
    [ z_W, w_W, ~ ] = gengausslegquadrule(K); % G-L abscissas and weights
    f_W_n = @(r) 0;

    for kk = 1:K
        zeta_W = @(x) m_G/Omega_G*(z_W(x)*Omega_L/m_L).^2;
        theta_W = @(x) w_W(x)/gamma(m_L)*zeta_W(x)^m_G*z_W(x)^(m_L-1);%-> psi_k
        alpha_W = @(x) 0; %-> xi_k
        
        for ii = 1:K
            alpha_W = @(x) alpha_W(x) + theta_W(ii)*zeta_W(ii)^(-m_G);
        end
        
        alpha_W = @(x) theta_W(x)/alpha_W(x);
        
        f_W_n = @(r) f_W_n(r) + alpha_W(kk)/gamma(m_G).*r.^(m_G-1).*exp(-zeta_W(kk)*r);
    end

    % Laplace transform for individual IRS
    L_W_fit_individual = @(s) ( integral(@(r) exp(-s*r).*f_W_n(r), 0, Inf) )^N;
    
    % Multiply individual Laplace transforms
    L_W_fit = @(s) L_W_fit(s) .* L_W_fit_individual(s);
end

L_W_exact = @(s) mean( exp(-s*Z_sim) );
 
zz = linspace(0, max(max(Z_sim)), 100);
waitbar_h = waitbar(0, 'Calculating Laplace Transform...');
for tt = 1:length(zz)
    point_L_W_exact(tt) = L_W_exact(zz(tt));
    point_L_W_fit(tt) = L_W_fit(zz(tt));
    waitbar(tt/length(zz), waitbar_h);
end
close(waitbar_h);

figure(1);
plot(zz, point_L_W_exact, 'DisplayName', 'Exact'); hold on;
plot(zz, point_L_W_fit, '+', 'DisplayName', 'Approximated'); hold on;
ylabel('Laplace Transform');
xlabel('z');
legend('location', 'best');
title(['Laplace Transform for Swarm of ' num2str(swarm_size) ' IRS UAVs']);

% STEP-5: DERIVE CDF and PDF pf \sum_{r=1}^N( \tilde{W}_{r} )
[~, setInd] = nsumk(K, N*swarm_size);

F_Z_fit = @(x) 0;

waitbar_h = waitbar(0, 'Deriving CDF and PDF...');
for caseIndex = 1:size(setInd, 1)
    indices = setInd(caseIndex, :);
    
    prodAlp = 1;
    for kk = 1:K
        vk = indices(kk);
        prodAlp = prodAlp * alpha_W(kk)^vk;
    end
    
    F_Z_fit = @(x) F_Z_fit(x) + factorial(N*swarm_size)/prod(factorial(indices)) * prodAlp * Phi2(indices.*m_G, 1+m_G*N*swarm_size, -zeta_W(1:K)*x, 50);
    % Phi2 is Humbert function
    
    waitbar(caseIndex/size(setInd, 1), waitbar_h);
end
close(waitbar_h);

F_Z_fit = @(x) F_Z_fit(x).*x.^(m_G*N*swarm_size)/gamma(1+m_G*N*swarm_size);
F_Z2_fit= @(y) F_Z_fit(sqrt(y)); % F_Y (y) = F_X (sqrt(y))

%% Outage probability
waitbar_h = waitbar(0, 'Calculating Outage Probability...');
for ii = 1:length(snr_dB)
    snr = 10^(snr_dB(ii)/10);
    OP_sim(ii) = mean( Z2_sim < snr_th/snr );
    OP_sim_new(ii) = mean( Z2_sim_optimal_phase_shift < snr_th/snr );
    OP_ana(ii) = F_Z2_fit( snr_th/snr );
    waitbar(ii/length(snr_dB), waitbar_h);
end
close(waitbar_h);

%% Plotting results
figure(2);
semilogy(snr_dB, OP_sim, 'b+:', 'DisplayName', 'Simulation (Random Phase)'); hold on;
semilogy(snr_dB, OP_sim_new, 'r-', 'DisplayName', 'Simulation (Optimal Phase)'); hold on;
semilogy(snr_dB, OP_ana, 'g:o', 'DisplayName', 'Analysis (Approx.)'); hold on;

xlabel('Transmit SNR [dB]');
ylabel('Outage Probability');
legend('location', 'southwest');
title(['Outage Probability for Swarm of ' num2str(swarm_size) ' IRS UAVs']);
axis([-Inf Inf 10^(-4) 1]);