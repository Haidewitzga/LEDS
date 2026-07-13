function [theta_hat] = ...
    estimate_fir_rls3(u, y, order, theta0, alpha, draw )
    arguments
        u
        y
        order
        theta0
        alpha
        draw = false
    end
% Estimates an FIR model by repeatedly calling rls3_step.
%
% Model:
%   y(t) = b1*u(t-1) + ... + bn*u(t-n) + e(t)

    %% Prepare data
    u = u(:);
    y = y(:);
    theta0 = theta0(:);

    N = length(y);

    if length(u) ~= N
        error('Input and output must have the same length.');
    end

    if length(theta0) ~= order
        error('theta0 must contain exactly "order" elements.');
    end

    if alpha <= 0
        error('alpha must be strictly positive.');
    end

    %% Initial conditions
    theta_hat = theta0;
    Sinv = alpha * eye(order);

    theta_history = NaN(order, N);
    error_history = NaN(1, N);

    theta_history(:,1:order) = ...
        repmat(theta0, 1, order);

 %% Create live plot
    if draw
       
        figure;
    
        layout = tiledlayout(2,1, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');
    
        % Parameter estimates
        ax1 = nexttile(layout);
        hold(ax1, 'on');
        grid(ax1, 'on');
    
        theta_lines = gobjects(order,1);
        colors = lines(order);
        
        for k = 1:order
            theta_lines(k) = animatedline( ...
                ax1, ...
                'Color', colors(k,:), ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('\\theta_%d', k));
        end
        xlabel(ax1, 'Sample t');
        ylabel(ax1, 'Parameter estimate');
        title(ax1, 'Online RLS parameter estimates');
    
        % Predetermined x-axis
        xlim(ax1, [1 N]);
    
        % Fixed legend
        legend(ax1, ...
            'Location', 'northeast', ...
            'AutoUpdate', 'off');
    
        % Prediction error
        ax2 = nexttile(layout);
        hold(ax2, 'on');
        grid(ax2, 'on');
    
       error_line = animatedline( ...
        ax2, ...
        'Color', [0.7 0.12 0.15], ...
        'LineWidth', 1.2);
    
        xlabel(ax2, 'Sample t');
        ylabel(ax2, 'Prediction error e(t)');
        title(ax2, 'Online prediction error');
    
        % Predetermined x-axis
        xlim(ax2, [1 N]);
    end


    %% Process measurements sequentially
    for t = order+1:N

        % FIR regressor:
        % phi(t) = [u(t-1), ..., u(t-order)]'
        phi = u(t-1:-1:t-order);

        % Perform one RLS III step
        [theta_hat, Sinv, prediction_error, ~] = rls3_step(phi,y(t),theta_hat, Sinv);

  
        % Update parameter plot
        if draw
             % Store results
            theta_history(:,t) = theta_hat;
            error_history(t) = prediction_error;

            % Draw
            for k = 1:order
                addpoints( ...
                    theta_lines(k), ...
                    t, ...
                    theta_hat(k));
            end
    
            % Update error plot
            addpoints( ...
                error_line, ...
                t, ...
                prediction_error);
    
            drawnow limitrate;
        end
        
    end
end