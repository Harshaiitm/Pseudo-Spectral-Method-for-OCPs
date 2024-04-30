function [nodes, weights, D_LG] = LG_computations(N)
    %% Calculation of the LG-Nodes 
    % These nodes are defined by the roots of (P_N(x))

    syms z;
    expr = legendreP(N, z);         % Legendre Polynomial of order N

    %  finds the roots of the equation P_N = 0, 
    % The vpasolve function is used to solve the equation symbolically.
    nodes = double(vpasolve(expr));

    % Calculation of weights of LG-Nodes
    weights = zeros(N, 1);
    for i = 1:N
        % Substituting symbolic value with numeric value
        node_value = nodes(i);
        num = subs(diff(legendreP(N + 1, z), z), node_value);
        den = (N + 1) * subs(legendreP(N, z), node_value);
        weights(i, 1) = 2 / (1 - node_value ^ 2) * (num / den) ^ 2;
    end

    %% Calculation of D_matrix of LG-Nodes
    D_LG = zeros(N, N);

    % Diagonal elements
    for i = 1:N
        D_LG(i, i) = nodes(i) / (1 - nodes(i) ^ 2);
    end

    % Non-diagonal elements
    for i = 1:N
        for j = 1:N
            if i ~= j
                % Substituting symbolic values with numeric values
                node_i = nodes(i);
                node_j = nodes(j);
                D_LG(i, j) = subs(diff(legendreP(N, z), z), node_i) * node_i / ...
                             ((node_i - node_j) * subs(diff(legendreP(N, z), z), node_j) * node_j);
            end
        end
    end

    %% Plotting the nodes
    close all;
    figure('Position', [200, 200, 400, 100]);
    x = nodes;
    y = zeros(size(x)); % y-coordinate will be zero for all nodes
    plot(x, y, '*', 'MarkerSize', 7, 'MarkerEdgeColor', 'b','LineWidth', 1.5); % Plotting nodes with grey circles
    xlabel('Nodes(\tau)','FontWeight', 'bold');
    % ylabel('LGL','FontWeight', 'bold');
    title(['Legendre-Gauss Nodes for N = ', num2str(N)]);
    grid on;
    set(gca, 'YTickLabel', []);
    ylim([-0.1, 0.1]);
end
