clear; clc;

% --- Option Parameters (Table 6) ---
Smin = 0;
Smax = 45;
sigma = 0.3;
X = 15;
r = 0.02;
T = 0.5;

% --- Coordinate Transformation Parameter (Example 1 uses xi = 12) ---
xi = 12; 
c1 = asinh(xi * (Smin - X));
c2 = asinh(xi * (Smax - X));

grids = [20, 40, 80];

% --- Print Header with Clean Column Widths ---
fprintf('=========================================================================================\n');
fprintf('                                         Table 6                                         \n');
fprintf('=========================================================================================\n');
fprintf('%-48s %-12s %-16s %-10s %-12s\n', 'Scheme', 'Grid', 'Error', 'Order', 'Computing time');
fprintf('%s\n', repmat('-', 1, 98));

%% -------------------------------------------------------------------
%  PART 1: FOURTH-ORDER CENTRAL FINITE DIFFERENCE METHOD
%  -------------------------------------------------------------------
fprintf('%-48s\n', 'Fourth-order central finite difference method');

for idx = 1:length(grids)
    N = grids(idx);
    M = grids(idx);
    
    h_y = 1 / N;
    dt = T / M;
    
    y = linspace(0, 1, N+1)';
    y_int = y(2:end-1);
    
    S_int = X + (1/xi) * sinh(y_int * (c2 - c1) + c1);
    J_int = ((c2 - c1) / xi) * cosh(y_int * (c2 - c1) + c1);
    H_int = ((c2 - c1)^2 / xi) * sinh(y_int * (c2 - c1) + c1);
    
    tic;
    
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
    
    A = A / (12 * h_y);
    B = B / (12 * h_y^2);
    
    Alpha = 0.5 * sigma^2 * (S_int.^2) ./ (J_int.^2);
    Beta = r * S_int ./ J_int - 0.5 * sigma^2 * (S_int.^2) .* H_int ./ (J_int.^3);
    
    Alpha_diag = diag(Alpha);
    Beta_diag = diag(Beta);
    I = eye(N-1);
    
    P = Alpha_diag * B + Beta_diag * A - r * I;
    
    get_Q = @(tau) Alpha_diag * [10*X*exp(-r*tau)/(12*h_y^2); -X*exp(-r*tau)/(12*h_y^2); zeros(N-3,1)] ...
                 + Beta_diag * [-3*X*exp(-r*tau)/(12*h_y); X*exp(-r*tau)/(12*h_y); zeros(N-3,1)];
             
    V_hist = zeros(N-1, 4);
    V_hist(:, 1) = max(X - S_int, 0);
    
    for step = 1:2
        tau_prev = (step - 1) * dt;
        tau_curr = step * dt;
        LHS = I - 0.5 * dt * P;
        RHS = (I + 0.5 * dt * P) * V_hist(:, step) + 0.5 * dt * (get_Q(tau_curr) + get_Q(tau_prev));
        V_hist(:, step+1) = LHS \ RHS;
    end
    
    tau_3 = 3 * dt;
    LHS = (11/6) * I - dt * P;
    RHS = 3 * V_hist(:, 3) - 1.5 * V_hist(:, 2) + (1/3) * V_hist(:, 1) + dt * get_Q(tau_3);
    V_hist(:, 4) = LHS \ RHS;
    
    for j = 4:M
        tau_next = j * dt;
        LHS = (25/12) * I - dt * P;
        RHS = 4 * V_hist(:, 4) - 3 * V_hist(:, 3) + (4/3) * V_hist(:, 2) - 0.25 * V_hist(:, 1) + dt * get_Q(tau_next);
        
        V_next = LHS \ RHS;
        V_hist(:, 1:3) = V_hist(:, 2:4);
        V_hist(:, 4) = V_next;
    end
    
    V_final = V_hist(:, 4);
    comp_time = toc;
    
    d1 = (log(S_int / X) + (r + 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    d2 = (log(S_int / X) + (r - 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    V_exact = X * exp(-r * T) * normcdf(-d2) - S_int .* normcdf(-d1); 
    
    err_val = max(abs(V_final - V_exact));
    
    grid_str = sprintf('%dx%d', N, M);
    time_str = sprintf('%.4f s', comp_time);
    
    if err_val < 1e-3
        err_str = sprintf('%.4e', err_val);
    else
        err_str = sprintf('%.4f', err_val);
    end
    
    if idx == 1
        order_str = '-';
        central_prev = err_val;
    else
        order = log2(central_prev / err_val);
        order_str = sprintf('%.4f', order);
        central_prev = err_val;
    end
    
    fprintf('%-48s %-12s %-16s %-10s %-12s\n', '', grid_str, err_str, order_str, time_str);
end

fprintf('\n');

%% -------------------------------------------------------------------
%  PART 2: FOURTH-ORDER COMPACT FINITE DIFFERENCE METHOD
%  -------------------------------------------------------------------
fprintf('%-48s\n', 'Fourth-order compact finite difference method');

for idx = 1:length(grids)
    N = grids(idx);
    M = grids(idx);
    
    h_y = 1 / N;
    dt = T / M;
    
    y = linspace(0, 1, N+1)';
    y_int = y(2:end-1);
    
    S_int = X + (1/xi) * sinh(y_int * (c2 - c1) + c1);
    J_int = ((c2 - c1) / xi) * cosh(y_int * (c2 - c1) + c1);
    H_int = ((c2 - c1)^2 / xi) * sinh(y_int * (c2 - c1) + c1);
    
    tic;
    
    F = zeros(N-1, N-1); G = zeros(N-1, N-1);
    U = zeros(N-1, N-1); W = zeros(N-1, N-1);
    
    for i = 2:N-2
        F(i, i-1:i+1) = [1/4, 1, 1/4];
        G(i, i-1:i+1) = [-3/(4*h_y), 0, 3/(4*h_y)];
        U(i, i-1:i+1) = [1/10, 1, 1/10];
        W(i, i-1:i+1) = [6/(5*h_y^2), -12/(5*h_y^2), 6/(5*h_y^2)];
    end
    F(1, 1:2) = [1, 3]; G(1, 1:4) = [-17/(6*h_y), 3/(2*h_y), 3/(2*h_y), -1/(6*h_y)];
    U(1, 1:2) = [1, 10]; W(1, 1:5) = [145/(12*h_y^2), -76/(3*h_y^2), 29/(2*h_y^2), -4/(3*h_y^2), 1/(12*h_y^2)];
    F(N-1, N-2:N-1) = [3, 1]; G(N-1, N-4:N-1) = [1/(6*h_y), -3/(2*h_y), -3/(2*h_y), 17/(6*h_y)];
    U(N-1, N-2:N-1) = [10, 1]; W(N-1, N-5:N-1) = [1/(12*h_y^2), -4/(3*h_y^2), 29/(2*h_y^2), -76/(3*h_y^2), 145/(12*h_y^2)];
    
    D = F \ G;
    D2 = U \ W;
    
    Alpha = 0.5 * sigma^2 * (S_int.^2) ./ (J_int.^2);
    Beta = r * S_int ./ J_int - 0.5 * sigma^2 * (S_int.^2) .* H_int ./ (J_int.^3);
    
    Alpha_diag = diag(Alpha);
    Beta_diag = diag(Beta);
    I = eye(N-1);
    
    P = Alpha_diag * D2 + Beta_diag * D - r * I;
    
    V_hist = zeros(N-1, 4);
    V_hist(:, 1) = max(X - S_int, 0);
    
    for step = 1:2
        LHS = I - 0.5 * dt * P;
        RHS = (I + 0.5 * dt * P) * V_hist(:, step);
        V_hist(:, step+1) = LHS \ RHS;
    end
    
    LHS = (11/6) * I - dt * P;
    RHS = 3 * V_hist(:, 3) - 1.5 * V_hist(:, 2) + (1/3) * V_hist(:, 1);
    V_hist(:, 4) = LHS \ RHS;
    
    for j = 4:M
        LHS = (25/12) * I - dt * P;
        RHS = 4 * V_hist(:, 4) - 3 * V_hist(:, 3) + (4/3) * V_hist(:, 2) - 0.25 * V_hist(:, 1);
        
        V_next = LHS \ RHS;
        V_hist(:, 1:3) = V_hist(:, 2:4);
        V_hist(:, 4) = V_next;
    end
    
    V_final = V_hist(:, 4);
    comp_time = toc;
    
    d1 = (log(S_int / X) + (r + 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    d2 = (log(S_int / X) + (r - 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    V_exact = X * exp(-r * T) * normcdf(-d2) - S_int .* normcdf(-d1); 
    
    err_val = max(abs(V_final - V_exact));
    
    grid_str = sprintf('%dx%d', N, M);
    time_str = sprintf('%.4f s', comp_time);
    
    if err_val < 1e-3
        err_str = sprintf('%.4e', err_val);
    else
        err_str = sprintf('%.4f', err_val);
    end
    
    if idx == 1
        order_str = '-';
        compact_prev = err_val;
    else
        order = log2(compact_prev / err_val);
        order_str = sprintf('%.4f', order);
        compact_prev = err_val;
    end
    
    fprintf('%-48s %-12s %-16s %-10s %-12s\n', '', grid_str, err_str, order_str, time_str);
end