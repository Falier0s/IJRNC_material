function r_output = pw_ref(bnd, T, f, m, method, options)
    % pw_ref generates a piecewise constant reference signal using different distributions.
    %
    % INPUTS:
    %   bnd    - n×2 matrix with lower and upper bounds for each signal
    %   T      - Total duration of the signal
    %   f      - Sampling frequency
    %   m      - Number of time intervals (segments)
    %   method - Distribution type: 'uniform', 'exp', 'beta', 'custom' (default: 'uniform')
    %   options - Struct with additional parameters (only needed for 'custom' method)
    %       options.c - Closeness factor (default: 0.05, used only in 'custom' method)
    %
    % OUTPUT:
    %   r_output - n×(T*f) matrix containing the generated reference signals

    if nargin < 5
        method = "uniform"; % Default distribution
    end

    if nargin < 6
        options.c = 0.05; % Default closeness for 'custom' method
    elseif ~isfield(options, 'c')
        options.c = 0.05;
    end

    % Number of signals (rows)
    n = size(bnd, 1);

    % Extract lower and upper bounds
    l_bnd = bnd(:, 1);
    u_bnd = bnd(:, 2);

    % Generate values based on the selected method
    switch lower(method)
        case "uniform"
            r = (u_bnd - l_bnd) .* rand(n, m) + l_bnd;

        case "exp"
            lambda = 2;  
            x = exprnd(1/lambda, n, m);
            signs = randi([0 1], n, m) * 2 - 1;
            x = signs .* x;
            x = tanh(x); 
            r = l_bnd + (u_bnd - l_bnd) .* (0.5 * (sign(x) .* (1 - exp(-abs(x) * 3)) + 1));

        case "beta"
            a = 0.3; b = 0.3;  
            x = betarnd(a, b, n, m);
            r = x .* (u_bnd - l_bnd) + l_bnd;

        case "custom"
            % Generate uniform random values
            r = (u_bnd - l_bnd) .* rand(n, m) + l_bnd;
            
            % Select 10% of the elements randomly
            num_elements = numel(r);
            num_to_adjust = round(0.2 * num_elements);
            idx = randperm(num_elements, num_to_adjust);
            
            % Compute boundary ranges
            range = options.c * (u_bnd - l_bnd);  % Closeness range
            
            % Generate random boundary values
            boundary_values = zeros(size(idx));
            for k = 1:length(idx)
                [row, col] = ind2sub(size(r), idx(k));
                if rand < 0.5
                    boundary_values(k) = l_bnd(row) + range(row) * rand;
                else
                    boundary_values(k) = u_bnd(row) - range(row) * rand;
                end
            end

            % Assign boundary values to selected indices
            r(idx) = boundary_values;

        otherwise
            error("Invalid method. Choose from 'uniform', 'exp', 'beta', or 'custom'.");
    end

    % Total number of samples
    total_samples = round(T * f);

    % Samples per interval
    samples_per_interval = round(total_samples / m);

    % Initialize output matrix
    r_output = zeros(n, total_samples);

    % Populate the reference signal matrix
    for i = 1:m
        start_idx = (i - 1) * samples_per_interval + 1;
        end_idx = min(i * samples_per_interval, total_samples);
        r_output(:, start_idx:end_idx) = repmat(r(:, i), 1, end_idx - start_idx + 1);
    end
end




