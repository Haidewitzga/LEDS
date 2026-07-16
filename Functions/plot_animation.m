function update_plot = plot_animation(order, N)
% Creates an RLS animation and returns its update function.
% Needs the order and the number of samples to initialize the plot axis and
% legends

    figure;

    layout = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    %% Parameter estimates
    ax1 = nexttile(layout);
    hold(ax1, 'on');
    grid(ax1, 'on');

    theta_lines = gobjects(order, 1);
    colors = lines(order);

    for k = 1:order
        theta_lines(k) = animatedline(  ax1, 'Color', colors(k, :),  'LineWidth', 1.5,  'DisplayName', sprintf('\\theta_%d', k));
    end

    xlabel(ax1, 'Sample t');
    ylabel(ax1, 'Parameter estimate');
    title(ax1, 'Online RLS parameter estimates');

    xlim(ax1, [1, N]);

    legend(ax1, ...
        'Location', 'northeast', ...
        'AutoUpdate', 'off');

    %% Prediction error
    ax2 = nexttile(layout);
    hold(ax2, 'on');
    grid(ax2, 'on');

    error_line = animatedline( ...
        ax2, ...
        'Color', [0.7, 0.12, 0.15], ...
        'LineWidth', 1.2);

    xlabel(ax2, 'Sample t');
    ylabel(ax2, 'Prediction error e(t)');
    title(ax2, 'Online prediction error');

    xlim(ax2, [1, N]);

    %% Return the nested update function
    update_plot = @update_animation;
    % For every new sample this function is going to be called to draw a
    % new point
    function update_animation(state)

        for k = 1:order
            addpoints(theta_lines(k), state.t, state.theta_hat(k));
        end

        addpoints(error_line, state.t, state.prediction_error);

        drawnow limitrate;
    end
end

