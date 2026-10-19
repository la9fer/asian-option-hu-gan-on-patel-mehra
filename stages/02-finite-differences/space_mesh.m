function x = space_mesh(N, K, Smax, epsilon, alpha, beta_star)
    % SPACE_MESH Generates the piecewise uniform grid.
    x = zeros(N+1, 1);
    x(1) = 0; % Initial boundary at x = 0

    h = (K - epsilon) / (1 + (alpha / beta_star) * (N/4 - 1));

    for i = 1:N
        idx = i + 1; % Shift by 1 for MATLAB's 1 based indexing
        
        if i == 1
            x(idx) = h;
            
        elseif i >= 2 && i <= N/4
            step_size = h * (alpha / beta_star);
            x(idx) = x(idx-1) + step_size;
            
        elseif i == (N/4) + 1
            x(idx) = K + epsilon;
            
        elseif i >= (N/4) + 2 && i <= N
            step_size = (Smax - K - epsilon) / (3*N/4 - 1);
            x(idx) = x(idx-1) + step_size;
        end
    end
end