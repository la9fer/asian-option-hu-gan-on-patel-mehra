% Parameters
T = 1.0;                % Time to maturity
r = 0.06;               % Risk-free interest rate
K = 25;                 % Strike price
Smax = 100;        % Maximum asset price
epsillon = 1e-4;    % smoothing polynomial parameter
N = 64;                 % Spatial steps
M = 64;                 % Time steps

% bounds for the piecewise mesh
alpha = 0.04;           % sigma = 0.2(1+t*e^-x), min sigma is 0.2
beta_star = 0.06;     % Maximum r (which is constant at 0.06)

% Initialize Grids
% Call the mesh generation function
x = space_mesh(N, K, Smax, epsillon, alpha, beta_star);

tau = T / M;                               % Time step size
t_grid = linspace(0, T, M+1);     % Time vector from 0 to T

% Solution matrix U (Rows = Spatial nodes (asset prices), Columns = Time steps)
U = zeros(N+1, M+1); 

% Apply the Terminal Condition
% At t = T (the last column), the value is the smoothed payoff function
U(:, M+1) = pi_epsillon(x - K, epsillon);

% 4. Backward Time Stepping Loop
% Start from the end (T) and step backward to today (t=0)
for j = M:-1:1
    t_current = t_grid(j); % The time step we are solving for
    
    % Preallocate vectors for the tridiagonal matrix
    A_diag = zeros(N+1, 1); % Sub-diagonal (U_{i-1} coefficients)
    B_diag = ones(N+1, 1);  % Main-diagonal (U_i coefficients)
    C_diag = zeros(N+1, 1); % Super-diagonal (U_{i+1} coefficients)
    
    % The right-hand side is the known option values from the "future" step
    b = U(:, j+1); 
    
    % Assemble matrix coefficients for interior nodes
    for i = 2:N
        hi = x(i) - x(i-1);
        hip1 = x(i+1) - x(i);
        denom = hi + hip1;
        
        % Volatility for Test 1 (space and time dependent)
        sigma = 0.2 * (1 + t_current * exp(-x(i)));
        sigma2_x2 = (sigma^2) * (x(i)^2);
        
        % Central Difference Scheme
        A_diag(i) = -tau * ( sigma2_x2 / (hi * denom) - (r * x(i)) / denom );
        C_diag(i) = -tau * ( sigma2_x2 / (hip1 * denom) + (r * x(i)) / denom );
        B_diag(i) = 1 + tau * ( sigma2_x2 / (hi * hip1) + r );
    end
    
    % Create the Sparse Tridiagonal Matrix (M-Matrix)
    % We format the diagonals for MATLAB's spdiags function
    A_matrix = sparse(diag(B_diag) + diag(A_diag(2:end), -1) + diag(C_diag(1:end-1), 1));
    
    % 5. Apply Boundary Conditions (Equation 2.11 and 2.12)
    
    % Lower Boundary (x = 0): Option value is 0
    A_matrix(1, :) = 0;
    A_matrix(1, 1) = 1;
    b(1) = 0;
    
    % Upper Boundary (x = Smax): Value approaches Smax - K * discount_factor
    A_matrix(N+1, :) = 0;
    A_matrix(N+1, N+1) = 1;
    b(N+1) = Smax - K * exp(-r * (T - t_current));
    
    % 6. Solve the Linear System (A * U^j = U^{j+1})
    U(:, j) = A_matrix \ b;
end

% 7. Plot the Results 
[T_mesh, X_mesh] = meshgrid(t_grid, x);
figure;
surf(X_mesh, T_mesh, U, 'EdgeColor', 'none');
colormap jet;
xlabel('Share Price x');
ylabel('Time t');
zlabel('Option Value U');
title('Computed Option Value U for Test 1');
view(-45, 30);