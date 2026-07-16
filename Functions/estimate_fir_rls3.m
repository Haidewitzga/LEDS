function theta_hat = estimate_fir_rls3(u, y, n,  alpha, theta0, step_callback)

% Estimates the parameters of a discrete-time FIR model using the recursive
% least-squares algorithm implemented in rls3_step.
%
% Model:
%   y(t) = b1*u(t-1) + ... + bn*u(t-n) + e(t)
%
% INPUT
% u: measured input signal samples
% y: measured output signal samples
% order: FIR model order and number of parameters to estimate
% alpha: optional initial scaling factor of the inverse information matrix
% theta0: optional initial estimate of the FIR parameter vector
% step_callback: optional function handle called after each recursive update
%
% OUTPUT
% theta_hat: final estimate of the FIR parameter vector

% Although the complete input and output datasets are provided to the 
% function at initialization, the measurements are processed sequentially. 
% At each iteration, only the current output sample and the previously 
% available input samples are used. The implementation therefore simulates 
% an online identification scenario in which measurements become available 
% progressively, as would sensor data in a real-time application.

    arguments 
        u 
        y 
        n 
        alpha (1,1) {mustBePositive} = 1
        theta0 = zeros(n,1)
        step_callback function_handle = @(~) []
    end

    %% Prepare data
    % Force arguments to be column vectors
    u = u(:);
    y = y(:);
    theta0 = theta0(:);

    N = length(y);

    if length(u) ~= N
        error('Input and output must have the same length.');
    end

    if length(theta0) ~= n
        error('theta0 must contain exactly "order" elements.');
    end

    if alpha <= 0
        error('alpha must be strictly positive.');
    end

    %% Initial conditions
    theta_hat = theta0;
    Sinv = alpha * eye(n);

    %% Process measurements sequentially
    for t = n+1:N

        % FIR regressor:
        % phi(t) = [u(t-1), ..., u(t-order)]'
        phi = u(t-1:-1:t-n);

        % Perform one RLS III step (For further description go to rls3_step.mlx)
        [theta_hat, Sinv, prediction_error, ~] = rls3_step(phi, y(t), theta_hat, Sinv);

        
        % update animation based on correction
        state.t = t; 
        state.theta_hat = theta_hat; 
        state.prediction_error = prediction_error;
        step_callback(state);
        
    end
end

