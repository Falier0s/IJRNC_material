function r_output = var_ref(bnd, T, f, p)
    % var_ref generates a reference signal using harmonic functions.
    %
    % INPUTS:
    %   bnd          - n×2 matrix with lower and upper bounds for each signal
    %   T            - Total duration of the signal
    %   f            - Sampling frequency
    %   p - Number of harmonic functions to superimpose
    %
    % OUTPUT:
    %   r_output - n×(T*f) matrix containing the generated reference signals

    % Number of signals (rows)
    n = size(bnd, 1);

    % Extract lower and upper bounds
    l_bnd = bnd(:, 1);
    u_bnd = bnd(:, 2);

    % Total number of samples
    total_samples = round(T * f);

    % Time vector
    t = linspace(0, T, total_samples);

    % Initialize output matrix
    r_output = zeros(n, total_samples);

    for i = 1:n
        % Generate random amplitudes, frequencies, and phases
        a = rand(1, p);
        f_scale = 0.5 + (0.8 * p - 0.5) * rand(1, p);
        phi = (pi / 3) * rand(1, p);

        % Generate the signal for each harmonic and sum them up
        signal = zeros(1, total_samples);
        for j = 1:p
            signal = signal + a(j) * sin(2 * pi * (f_scale(j) / T) * t + phi(j));
        end

        % Scale the signal to fit within the desired bounds
        signal = signal / max(abs(signal)) * (1* abs(u_bnd(i) - l_bnd(i)) / 2);

        % Shift the signal to be within the bounds
        signal = signal + (u_bnd(i) + l_bnd(i)) / 2;

        % Store the generated signal in the output matrix
        r_output(i, :) = signal;
    end
end
