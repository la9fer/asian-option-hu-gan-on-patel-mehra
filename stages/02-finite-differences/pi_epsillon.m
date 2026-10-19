function result = pi_epsillon(y, epsillon)
    % PI_EPSILLON smooths the max(y, 0) function around y = 0.
    % y: Can be a scalar or a vector (e.g., the array of x_i - K)
    % epsillon: The small transition parameter (e.g., 1e-4)
    
    % Initialize output array of the same size as input y
    result = zeros(size(y));
    
    % The 10 polynomial coefficients derived in the paper
    c0 = (35/256) * epsillon;
    c1 = 1/2;
    c2 = 35 / (64 * epsillon);
    % c3 = 0
    c4 = -35 / (128 * epsillon^3);
    % c5 = 0
    c6 = 7 / (64 * epsillon^5);
    % c7 = 0
    c8 = -5 / (256 * epsillon^7);
    % c9 = 0
    
    % Vectorized Logical Indexing
    % 1. Where y >= epsillon, the function is just y
    idx_high = (y >= epsillon);
    result(idx_high) = y(idx_high);
    
    % 2. Where y <= -epsillon, the function is 0 (already handled by zeros() init)
    
    % 3. -epsillon < y < epsillon
    idx_mid = (y > -epsillon) & (y < epsillon);
    ym = y(idx_mid); % Extract just the points in the transition zone
    
    % Apply the 9th degree polynomial
    result(idx_mid) = c0 + c1.*ym + c2.*(ym.^2) + c4.*(ym.^4) + ...
                      c6.*(ym.^6) + c8.*(ym.^8);
end