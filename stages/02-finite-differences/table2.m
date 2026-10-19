% Error and Convergence Rates for Test 2

clc; clear;

N_exact = 2048;
M_exact = 1024;
[x_exact, U_exact] = run_black_scholes_test2(N_exact, M_exact);

% Define the N values to test
N_list = [128, 256, 512, 1024];
num_tests = length(N_list);

errors = zeros(num_tests, 1);
rates = zeros(num_tests, 1);

for k = 1:num_tests
    N = N_list(k);
    [x_coarse, U_coarse] = run_black_scholes_test2(N, M_exact);
    
    U_exact_interp = zeros(N+1, M_exact+1);
    
    for j = 1:(M_exact+1)
        U_exact_interp(:, j) = interp1(x_exact, U_exact(:, j), x_coarse, 'linear');
    end
    
    errors(k) = max(max(abs(U_coarse - U_exact_interp)));
end

for k = 1:(num_tests - 1)
    rates(k) = log2(errors(k) / errors(k+1));
end

% Print the final table
fprintf('\n--- Table 2 ---\n');
fprintf('%-6s %-6s %-15s %-10s\n', 'M', 'N', 'Error e^{N,M}', 'Rate r^{N,M}');
fprintf('----------------------------------------\n');
for k = 1:num_tests
    if k < num_tests
        fprintf('%-6d %-6d %-15.4e %-10.3f\n', M_exact, N_list(k), errors(k), rates(k));
    else
        fprintf('%-6d %-6d %-15.4e %-10s\n', M_exact, N_list(k), errors(k), '-');
    end
end
fprintf('----------------------------------------\n');

% LOCAL HELPER FUNCTION: The Test 2 Solver
function [x, U] = run_black_scholes_test2(N, M)
    % 1. Define Test 2 Parameters
    T = 1.0;            
    r = 0.06;           
    K = 25;             
    Smax = 100;         
    epsillon = 1e-4;     
    
    % UPDATED: Alpha bound for Test 2
    alpha = 0.16;       % min(sigma)^2 = (0.4 * (2 - 1))^2 = 0.16
    beta_star = 0.06;   
    
    tau = T / M;                      
    t_grid = linspace(0, T, M+1);     

    % 2. Initialize Grids
    x = space_mesh(N, K, Smax, epsillon, alpha, beta_star);
    U = zeros(N+1, M+1); 

    % 3. Apply Terminal Condition
    U(:, M+1) = pi_epsillon(x - K, epsillon);

    % 4. Backward Time-Stepping Loop
    for j = M:-1:1
        t_current = t_grid(j);
        
        A_diag = zeros(N+1, 1); 
        B_diag = ones(N+1, 1);  
        C_diag = zeros(N+1, 1); 
        b = U(:, j+1); 
        
        for i = 2:N
            hi = x(i) - x(i-1);
            hip1 = x(i+1) - x(i);
            denom = hi + hip1;
            
            % UPDATED: Test 2 Volatility formula
            sigma = 0.4 * (2 + sin(x(i)));
            sigma2_x2 = (sigma^2) * (x(i)^2);
            
            A_diag(i) = -tau * ( sigma2_x2 / (hi * denom) - (r * x(i)) / denom );
            C_diag(i) = -tau * ( sigma2_x2 / (hip1 * denom) + (r * x(i)) / denom );
            B_diag(i) = 1 + tau * ( sigma2_x2 / (hi * hip1) + r );
        end
        
        A_matrix = sparse(diag(B_diag) + diag(A_diag(2:end), -1) + diag(C_diag(1:end-1), 1));
        
        % 5. Apply Boundary Conditions
        A_matrix(1, :) = 0; A_matrix(1, 1) = 1; b(1) = 0;
        A_matrix(N+1, :) = 0; A_matrix(N+1, N+1) = 1; 
        b(N+1) = Smax - K * exp(-r * (T - t_current));
        
        % 6. Solve
        U(:, j) = A_matrix \ b;
    end
end