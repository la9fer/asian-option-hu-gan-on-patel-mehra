clear; clc;

% --- Option Parameters (Table 4) ---
Smin = 0;
Smax = 45;
sigma = 0.3;
X = 15;
r = 0.02;
T = 0.5;

grids = [10, 20, 40];

fprintf('=====================================================\n');
fprintf('                     Table 4                        \n');
fprintf('=====================================================\n');

fprintf('Scheme: Crank-Nicolson\n');
fprintf('Grid\t\tError\t\tOrder\n');
fprintf('-----------------------------------------\n');

cn_errors = zeros(length(grids), 1);

for idx = 1:length(grids)
    N = grids(idx);
    M = grids(idx);
    
    h = (Smax - Smin) / N;
    dt = T / M;
    
    S = linspace(Smin, Smax, N+1)';
    S_int = S(2:end-1);
    
    % --- Matrices A and B (Eq 3.9) ---
    A = zeros(N-1, N-1);
    B = zeros(N-1, N-1);
    
    for i = 2:N-2
        A(i, i-1:i+1) = [-8, 0, 8];
        B(i, i-1:i+1) = [16, -30, 16];
        if i > 2, A(i, i-2) = 1; B(i, i-2) = -1; end
        if i < N-2, A(i, i+2) = -1; B(i, i+2) = -1; end
    end
    
    A(1, 1:4)       = [-10, 18, -6, 1];
    A(N-1, N-4:N-1) = [-1, 6, -18, 10]; 
    
    B(1, 1:5)       = [-15, -4, 14, -6, 1];
    B(N-1, N-5:N-1) = [1, -6, 14, -4, -15];
    
    A = A / (12 * h);
    B = B / (12 * h^2);
    
    S_diag = diag(S_int);
    S2_diag = diag(S_int.^2);
    I = eye(N-1);
    
    P = 0.5 * sigma^2 * S2_diag * B + r * S_diag * A - r * I;
    
    % Initial condition for European Put at tau = 0
    V = max(X - S_int, 0);
    
    % Time Loop: Crank-Nicolson
    for j = 1:M
        tau_prev = (j - 1) * dt;
        tau_curr = j * dt;
        
        V0_prev = X * exp(-r * tau_prev);
        V0_curr = X * exp(-r * tau_curr);
        VN_prev = 0;
        VN_curr = 0;
        
        b1_prev = zeros(N-1, 1); b1_curr = zeros(N-1, 1);
        b2_prev = zeros(N-1, 1); b2_curr = zeros(N-1, 1);
        
        b1_prev(1) = -3 * V0_prev / (12*h); b1_prev(2) = V0_prev / (12*h);
        b1_curr(1) = -3 * V0_curr / (12*h); b1_curr(2) = V0_curr / (12*h);
        
        b2_prev(1) = 10 * V0_prev / (12*h^2); b2_prev(2) = -V0_prev / (12*h^2);
        b2_curr(1) = 10 * V0_curr / (12*h^2); b2_curr(2) = -V0_curr / (12*h^2);
        
        Q_prev = 0.5 * sigma^2 * S2_diag * b2_prev + r * S_diag * b1_prev;
        Q_curr = 0.5 * sigma^2 * S2_diag * b2_curr + r * S_diag * b1_curr;
        
        LHS = I - 0.5 * dt * P;
        RHS = (I + 0.5 * dt * P) * V + 0.5 * dt * (Q_curr + Q_prev);
        
        V = LHS \ RHS;
    end
    
    % Analytical Solution
    d1 = (log(S_int / X) + (r + 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    d2 = (log(S_int / X) + (r - 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    V_exact = X * exp(-r * T) * normcdf(-d2) - S_int .* normcdf(-d1); 
    
    cn_errors(idx) = max(abs(V - V_exact));
    
    if idx == 1
        fprintf('%dx%d\t\t%.4f\t\t-\n', N, M, cn_errors(idx));
    else
        order = log2(cn_errors(idx-1) / cn_errors(idx));
        fprintf('%dx%d\t\t%.4f\t\t%.4f\n', N, M, cn_errors(idx), order);
    end
end

fprintf('\n');

fprintf('Scheme: BDF4\n');
fprintf('Grid\t\tError\t\tOrder\n');
fprintf('-----------------------------------------\n');

bdf_errors = zeros(length(grids), 1);

for idx = 1:length(grids)
    N = grids(idx);
    M = grids(idx);
    
    h = (Smax - Smin) / N;
    dt = T / M;
    
    S = linspace(Smin, Smax, N+1)';
    S_int = S(2:end-1);
    
    % --- Matrices A and B (Eq 3.9) ---
    A = zeros(N-1, N-1);
    B = zeros(N-1, N-1);
    
    for i = 2:N-2
        A(i, i-1:i+1) = [-8, 0, 8];
        B(i, i-1:i+1) = [16, -30, 16];
        if i > 2, A(i, i-2) = 1; B(i, i-2) = -1; end
        if i < N-2, A(i, i+2) = -1; B(i, i+2) = -1; end
    end
    
    A(1, 1:4)       = [-10, 18, -6, 1];
    A(N-1, N-4:N-1) = [-1, 6, -18, 10]; 
    
    B(1, 1:5)       = [-15, -4, 14, -6, 1];
    B(N-1, N-5:N-1) = [1, -6, 14, -4, -15];
    
    A = A / (12 * h);
    B = B / (12 * h^2);
    
    S_diag = diag(S_int);
    S2_diag = diag(S_int.^2);
    I = eye(N-1);
    
    P = 0.5 * sigma^2 * S2_diag * B + r * S_diag * A - r * I;
    
    % Helper function to calculate Q(tau)
    get_Q = @(tau) 0.5 * sigma^2 * S2_diag * [10*X*exp(-r*tau)/(12*h^2); -X*exp(-r*tau)/(12*h^2); zeros(N-3,1)] ...
                 + r * S_diag * [-3*X*exp(-r*tau)/(12*h); X*exp(-r*tau)/(12*h); zeros(N-3,1)];
             
    % Storage for history: [u^(j-3), u^(j-2), u^(j-1), u^(j)]
    V_hist = zeros(N-1, 4);
    V_hist(:, 1) = max(X - S_int, 0); % u^(0)
    
    % Step 1 & 2: Crank-Nicolson for u^(1) and u^(2)
    for step = 1:2
        tau_prev = (step - 1) * dt;
        tau_curr = step * dt;
        LHS_cn = I - 0.5 * dt * P;
        RHS_cn = (I + 0.5 * dt * P) * V_hist(:, step) + 0.5 * dt * (get_Q(tau_curr) + get_Q(tau_prev));
        V_hist(:, step+1) = LHS_cn \ RHS_cn;
    end
    
    % Step 3: BDF3 for u^(3) (Eq 4.10)
    tau_3 = 3 * dt;
    LHS_bdf3 = (11/6) * I - dt * P;
    RHS_bdf3 = 3 * V_hist(:, 3) - 1.5 * V_hist(:, 2) + (1/3) * V_hist(:, 1) + dt * get_Q(tau_3);
    V_hist(:, 4) = LHS_bdf3 \ RHS_bdf3;
    
    % Step 4 to M: BDF4 for u^(4) to u^(M) (Eq 4.9)
    for j = 4:M
        tau_next = j * dt;
        LHS_bdf4 = (25/12) * I - dt * P;
        RHS_bdf4 = 4 * V_hist(:, 4) - 3 * V_hist(:, 3) + (4/3) * V_hist(:, 2) - 0.25 * V_hist(:, 1) + dt * get_Q(tau_next);
        
        V_next = LHS_bdf4 \ RHS_bdf4;
        
        % Shift window
        V_hist(:, 1:3) = V_hist(:, 2:4);
        V_hist(:, 4) = V_next;
    end
    
    V_final = V_hist(:, 4);
    
    % Analytical Solution
    d1 = (log(S_int / X) + (r + 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    d2 = (log(S_int / X) + (r - 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    V_exact = X * exp(-r * T) * normcdf(-d2) - S_int .* normcdf(-d1); 
    
    bdf_errors(idx) = max(abs(V_final - V_exact));
    
    if idx == 1
        fprintf('%dx%d\t\t%.4f\t\t-\n', N, M, bdf_errors(idx));
    else
        order = log2(bdf_errors(idx-1) / bdf_errors(idx));
        fprintf('%dx%d\t\t%.4f\t\t%.4f\n', N, M, bdf_errors(idx), order);
    end
end