% Error and Convergence Rates for Test 1

clc; clear;

N_exact = 2048;
M_exact = 1024;
[x_exact, U_exact] = run_black_scholes(N_exact, M_exact);

% Define the N values we want to test (M is fixed at 1024)
N_list = [128, 256, 512, 1024];
num_tests = length(N_list);

errors = zeros(num_tests, 1);
rates = zeros(num_tests, 1);

for k = 1:num_tests
    N = N_list(k);
    [x_coarse, U_coarse] = run_black_scholes(N, M_exact);
    
    % Preallocate array for the interpolated exact solution
    U_exact_interp = zeros(N+1, M_exact+1);
    
    % Interpolate the exact solution onto the coarse spatial grid for each time step
    for j = 1:(M_exact+1)
        U_exact_interp(:, j) = interp1(x_exact, U_exact(:, j), x_coarse, 'linear');
    end
    
    % Calculate the maximum absolute error across all space and time
    errors(k) = max(max(abs(U_coarse - U_exact_interp)));
end

% Calculate convergence rates: log2( Error(N) / Error(2N) )
for k = 1:(num_tests - 1)
    rates(k) = log2(errors(k) / errors(k+1));
end

% Print the final table to the Command Window
fprintf('\n--- Table 1 ---\n');
fprintf('%-6s %-6s %-15s %-10s\n', 'M', 'N', 'Error e^{N,M}', 'Rate r^{N,M}');
fprintf('----------------------------------------\n');
for k = 1:num_tests
    if k < num_tests
        fprintf('%-6d %-6d %-15.4e %-10.3f\n', M_exact, N_list(k), errors(k), rates(k));
    else
        % The last row doesn't have a rate calculation
        fprintf('%-6d %-6d %-15.4e %-10s\n', M_exact, N_list(k), errors(k), '-');
    end
end
fprintf('----------------------------------------\n');


% LOCAL HELPER FUNCTION: The Core Solver
function [x, U] = run_black_scholes(N, M)
    % 1. Define Test 1 Parameters
    T = 1.0;            
    r = 0.06;           
    K = 25;             
    Smax = 100;         
    epsillon = 1e-4;     
    alpha = 0.04;       
    beta_star = 0.06;   
    tau = T / M;                      
    t_grid = linspace(0, T, M+1);     

    % 2. Initialize Grids
    % Make sure the filename matches what you saved the mesh function as!
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
            
            % Test 1 Volatility
            sigma = 0.2 * (1 + t_current * exp(-x(i)));
            sigma2_x2 = (sigma^2) * (x(i)^2);
            
            A_diag(i) = -tau * ( sigma2_x2 / (hi * denom) - (r * x(i)) / denom );
            C_diag(i) = -tau * ( sigma2_x2 / (hip1 * denom) + (r * x(i)) / denom );
            B_diag(i) = 1 + tau * ( sigma2_x2 / (hi * hip1) + r );
        end
        
        % The corrected M-Matrix assembly
        A_matrix = sparse(diag(B_diag) + diag(A_diag(2:end), -1) + diag(C_diag(1:end-1), 1));
        
        % 5. Apply Boundary Conditions
        A_matrix(1, :) = 0; A_matrix(1, 1) = 1; b(1) = 0;
        A_matrix(N+1, :) = 0; A_matrix(N+1, N+1) = 1; 
        b(N+1) = Smax - K * exp(-r * (T - t_current));
        
        % 6. Solve
        U(:, j) = A_matrix \ b;
    end
end