clear; clc;

% --- Option Parameters (Table 5) ---
Smin = 0;
Smax = 45;
sigma = 0.3;
X = 15;
r = 0.02;
T = 0.5;

fprintf('=====================================================\n');
fprintf('                     Table 5                        \n');
fprintf('=====================================================\n');

grids_cn = [20, 40, 80];
cn_errors = zeros(length(grids_cn), 1);

fprintf('Scheme: Crank-Nicolson\n');
fprintf('Grid\t\tError\t\tOrder\n');
fprintf('-----------------------------------------\n');

for idx = 1:length(grids_cn)
    N = grids_cn(idx);
    M = grids_cn(idx);
    
    h = (Smax - Smin) / N;
    dt = T / M;
    
    S = linspace(Smin, Smax, N+1)';
    S_int = S(2:end-1);
    
    % --- Matrices F, G, U, W for Compact Scheme (Eq 3.21 & 3.22) ---
    F = zeros(N-1, N-1);
    G = zeros(N-1, N-1);
    U = zeros(N-1, N-1);
    W = zeros(N-1, N-1);
    
    for i = 2:N-2
        F(i, i-1:i+1) = [1/4, 1, 1/4];
        G(i, i-1:i+1) = [-3/(4*h), 0, 3/(4*h)];
        U(i, i-1:i+1) = [1/10, 1, 1/10];
        W(i, i-1:i+1) = [6/(5*h^2), -12/(5*h^2), 6/(5*h^2)];
    end
    
    F(1, 1:2) = [1, 3];
    G(1, 1:4) = [-17/(6*h), 3/(2*h), 3/(2*h), -1/(6*h)];
    U(1, 1:2) = [1, 10];
    W(1, 1:5) = [145/(12*h^2), -76/(3*h^2), 29/(2*h^2), -4/(3*h^2), 1/(12*h^2)];
    
    F(N-1, N-2:N-1) = [3, 1];
    G(N-1, N-4:N-1) = [1/(6*h), -3/(2*h), -3/(2*h), 17/(6*h)];
    U(N-1, N-2:N-1) = [10, 1];
    W(N-1, N-5:N-1) = [1/(12*h^2), -4/(3*h^2), 29/(2*h^2), -76/(3*h^2), 145/(12*h^2)];
    
    D = F \ G;
    D2 = U \ W;
    
    S_diag = diag(S_int);
    S2_diag = diag(S_int.^2);
    I = eye(N-1);
    
    P = 0.5 * sigma^2 * S2_diag * D2 + r * S_diag * D - r * I;
    
    % Initial condition for European Put at tau = 0
    V = max(X - S_int, 0);
    
    % Time Loop: Crank-Nicolson
    for j = 1:M
        % Compact scheme encapsulates boundaries internally for this setup
        LHS = I - 0.5 * dt * P;
        RHS = (I + 0.5 * dt * P) * V;
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

%% -------------------------------------------------------------------
%  PART 2: BDF4 METHOD (Grids: 10, 20, 40)
%  -------------------------------------------------------------------
grids_bdf = [20, 40, 80];
bdf_errors = zeros(length(grids_bdf), 1);

fprintf('Scheme: BDF4\n');
fprintf('Grid\t\tError\t\tOrder\n');
fprintf('-----------------------------------------\n');

for idx = 1:length(grids_bdf)
    N = grids_bdf(idx);
    M = grids_bdf(idx);
    
    h = (Smax - Smin) / N;
    dt = T / M;
    
    S = linspace(Smin, Smax, N+1)';
    S_int = S(2:end-1);
    
    % --- Matrices F, G, U, W for Compact Scheme (Eq 3.21 & 3.22) ---
    F = zeros(N-1, N-1);
    G = zeros(N-1, N-1);
    U = zeros(N-1, N-1);
    W = zeros(N-1, N-1);
    
    for i = 2:N-2
        F(i, i-1:i+1) = [1/4, 1, 1/4];
        G(i, i-1:i+1) = [-3/(4*h), 0, 3/(4*h)];
        U(i, i-1:i+1) = [1/10, 1, 1/10];
        W(i, i-1:i+1) = [6/(5*h^2), -12/(5*h^2), 6/(5*h^2)];
    end
    
    F(1, 1:2) = [1, 3];
    G(1, 1:4) = [-17/(6*h), 3/(2*h), 3/(2*h), -1/(6*h)];
    U(1, 1:2) = [1, 10];
    W(1, 1:5) = [145/(12*h^2), -76/(3*h^2), 29/(2*h^2), -4/(3*h^2), 1/(12*h^2)];
    
    F(N-1, N-2:N-1) = [3, 1];
    G(N-1, N-4:N-1) = [1/(6*h), -3/(2*h), -3/(2*h), 17/(6*h)];
    U(N-1, N-2:N-1) = [10, 1];
    W(N-1, N-5:N-1) = [1/(12*h^2), -4/(3*h^2), 29/(2*h^2), -76/(3*h^2), 145/(12*h^2)];
    
    D = F \ G;
    D2 = U \ W;
    
    S_diag = diag(S_int);
    S2_diag = diag(S_int.^2);
    I = eye(N-1);
    
    P = 0.5 * sigma^2 * S2_diag * D2 + r * S_diag * D - r * I;
    
    % Storage for history: [u^(j-3), u^(j-2), u^(j-1), u^(j)]
    V_hist = zeros(N-1, 4);
    V_hist(:, 1) = max(X - S_int, 0); % u^(0)
    
    % Step 1 & 2: Crank-Nicolson for u^(1) and u^(2)
    for step = 1:2
        LHS_cn = I - 0.5 * dt * P;
        RHS_cn = (I + 0.5 * dt * P) * V_hist(:, step);
        V_hist(:, step+1) = LHS_cn \ RHS_cn;
    end
    
    % Step 3: BDF3 for u^(3) (Eq 4.10)
    LHS_bdf3 = (11/6) * I - dt * P;
    RHS_bdf3 = 3 * V_hist(:, 3) - 1.5 * V_hist(:, 2) + (1/3) * V_hist(:, 1);
    V_hist(:, 4) = LHS_bdf3 \ RHS_bdf3;
    
    % Step 4 to M: BDF4 for u^(4) to u^(M) (Eq 4.9)
    for j = 4:M
        LHS_bdf4 = (25/12) * I - dt * P;
        RHS_bdf4 = 4 * V_hist(:, 4) - 3 * V_hist(:, 3) + (4/3) * V_hist(:, 2) - 0.25 * V_hist(:, 1);
        
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