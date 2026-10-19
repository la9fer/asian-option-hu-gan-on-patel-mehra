% Table 1: Black-Scholes 2nd Order Crank-Nicolson on Equidistant Grid
clear; clc;

% 1. Option Parameters (from Section 5.1)
K = 15;        % Exercise price
sigma = 0.3;   % Volatility
r = 0.05;      % Interest rate
rd = 0.03;     % Continuous dividend yield
T = 0.5;       % Time to maturity
Smax = 3 * K;  % Outer boundary placed at 3 times strike

% Grid setups for Table 1 (10x10, 20x20, 40x40)
N_values = [10, 20, 40]; 

% Storage for errors to calculate reduction factor (c_inf)
errors_U = zeros(length(N_values), 1);
errors_D = zeros(length(N_values), 1);
errors_G = zeros(length(N_values), 1);

fprintf('Grid \t ||u_h - u_ex|| \t c_inf \t ||D_h - D_ex|| \t c_inf \t ||G_h - G_ex|| \t c_inf\n');

for idx = 1:length(N_values)
    N = N_values(idx);
    M = N; % Space steps = Time steps
    
    % 2. Discretization Steps
    h = Smax / N; % Spatial step size
    k = T / M;    % Time step size
    
    S = linspace(0, Smax, N+1)'; % Stock price grid
    
    % 3. Exact Analytical Solution for comparison (at t=0)
    d1 = (log(S(2:end-1)/K) + (r - rd + 0.5*sigma^2)*T) / (sigma*sqrt(T));
    d2 = d1 - sigma*sqrt(T);
    
    U_exact = S(2:end-1) * exp(-rd*T) .* normcdf(d1) - K * exp(-r*T) * normcdf(d2);
    Delta_exact = exp(-rd*T) * normcdf(d1);
    Gamma_exact = exp(-rd*T) * normpdf(d1) ./ (S(2:end-1) * sigma * sqrt(T));
    
    % 4. Setting up the Tridiagonal Matrix Operator (L)
    % Central difference coefficients
    alpha = 0.5 * sigma^2 * S.^2 / h^2;
    beta  = (r - rd) * S / (2*h);
    
    % Diagonals for the operator matrix
    lower_diag = alpha(3:N) - beta(3:N);
    main_diag  = -2*alpha(2:N) - r;
    upper_diag = alpha(2:N-1) + beta(2:N-1);
    
    L = diag(main_diag) + diag(lower_diag, -1) + diag(upper_diag, 1);
    
    % Crank-Nicolson Matrices: (I - 0.5*k*L) * U_new = (I + 0.5*k*L) * U_old
    I = eye(N-1);
    A = I - 0.5 * k * L;
    B = I + 0.5 * k * L;
    
    % 5. Initial Condition (Payoff at Maturity)
    % We solve backwards, so t=T is our starting point (tau = 0)
    U = max(S - K, 0);
    U_interior = U(2:N); % We only solve for interior points
    
    % 6. Time Stepping Loop
    for j = 1:M
        tau = j * k;
        
        % Boundary Conditions at current and next time step
        % Left boundary (S=0): U = 0 for a Call
        % Right boundary (S=Smax): U = Smax*exp(-rd*tau) - K*exp(-r*tau)
        bc_old = zeros(N-1, 1);
        bc_new = zeros(N-1, 1);
        
        bc_old(end) = (alpha(N) + beta(N)) * (Smax*exp(-rd*(tau-k)) - K*exp(-r*(tau-k)));
        bc_new(end) = (alpha(N) + beta(N)) * (Smax*exp(-rd*tau) - K*exp(-r*tau));
        
        % Right-hand side of the linear system
        rhs = B * U_interior + 0.5 * k * (bc_old + bc_new);
        
        % Solve the system
        U_interior = A \ rhs;
    end
    
    % 7. Calculate Greeks (Numerical Derivatives at t=0)
    % Pad the interior solution with the final boundary conditions
    U_final = [0; U_interior; Smax*exp(-rd*T) - K*exp(-r*T)];
    
    Delta_num = zeros(N-1, 1);
    Gamma_num = zeros(N-1, 1);
    
    for i = 2:N
        Delta_num(i-1) = (U_final(i+1) - U_final(i-1)) / (2*h);
        Gamma_num(i-1) = (U_final(i+1) - 2*U_final(i) + U_final(i-1)) / h^2;
    end
    
    % 8. Calculate L-infinity Errors
    errors_U(idx) = max(abs(U_interior - U_exact));
    errors_D(idx) = max(abs(Delta_num - Delta_exact));
    errors_G(idx) = max(abs(Gamma_num - Gamma_exact));
    
    % 9. Calculate Error Reduction Factor (c_inf) and Print
    if idx == 1
        fprintf('%dx%d \t %.1e \t \t - \t %.1e \t \t - \t %.1e \t \t -\n', ...
            N, M, errors_U(idx), errors_D(idx), errors_G(idx));
    else
        c_U = errors_U(idx-1) / errors_U(idx);
        c_D = errors_D(idx-1) / errors_D(idx);
        c_G = errors_G(idx-1) / errors_G(idx);
        fprintf('%dx%d \t %.1e \t \t %.1f \t %.1e \t \t %.1f \t %.1e \t \t %.1f\n', ...
            N, M, errors_U(idx), c_U, errors_D(idx), c_D, errors_G(idx), c_G);
    end
end