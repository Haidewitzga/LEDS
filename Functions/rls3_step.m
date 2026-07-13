function [theta_new, Sinv_new, prediction_error, K] = ...
    rls3_step(phi, y_current, theta_old, Sinv_old)

    phi = phi(:);

    % One-step prediction error
    prediction_error = ...
        y_current - phi' * theta_old;

    % Denominator from the matrix inversion lemma
    denominator = 1 + phi' * Sinv_old * phi;

    % Recursive gain
    K = (Sinv_old * phi) / denominator;

    % Parameter estimate update
    theta_new = theta_old + K * prediction_error;

    % Recursive update of S^(-1)(t)
    Sinv_new = Sinv_old  - (Sinv_old * phi * phi' * Sinv_old) / denominator;
end