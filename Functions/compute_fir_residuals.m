function residuals = compute_fir_residuals(u, y, theta_hat)
% Computes the one-step prediction residuals of an estimated FIR model.
%
% INPUT
% u: measured input signal samples
% y: measured output signal samples
% theta_hat: estimated FIR parameter vector
%
% OUTPUT
% residuals: one-step prediction residuals of the estimated FIR model

    u = u(:);
    y = y(:);
    theta_hat = theta_hat(:);

    n = numel(theta_hat);
    N = numel(y);

    if numel(u) ~= N
        error('Input and output signals must have the same length.');
    end

    residuals = zeros(N - n, 1);

    for t = n+1:N
        phi = u(t-1:-1:t-n);

        y_hat = phi' * theta_hat;

        residuals(t-n) = y(t) - y_hat;
    end
end