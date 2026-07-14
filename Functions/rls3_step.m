function [theta_new, Sinv_new, prediction_error, K] = rls3_step(phi, y_current, theta_old, Sinv_old)
% Performs one Recursive Least Squares III update using the matrix inversion lemma.
% The function processes a single measurement sample and updates the
% parameter estimate and inverse information matrix

    % Force phi to be a column vector
    phi = phi(:);

    % One-step-ahead prediction error
    prediction_error = y_current - phi' * theta_old;

    % Denominator obtained from the matrix inversion lemma
    denominator = 1 + phi' * Sinv_old * phi;

    % Recursive gain
    K = (Sinv_old * phi) / denominator;

    % Parameter estimate update
    theta_new = theta_old + K * prediction_error;

    % Recursive update of the inverse information matrix S^(-1)(t)
    Sinv_new = Sinv_old  - (Sinv_old * phi * phi' * Sinv_old) / denominator;
end

