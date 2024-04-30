function [nodes, weights, D_CGL] = CGL_computations(N)
%% Calculation of the LG-Nodes 
%--------------------------------------------------------------------------
% CGL_nodes.m
% determines Chebyshev-Gauss-Lobatto (CGL) nodes
%--------------------------------------------------------------------------
% tau = CGL_nodes(M)  
%   M: number of nodes(N) minus 1, should be an integer greater than 0
% tau: CGL nodes
%--------------------------------------------------------------------------
% Examples: M is the order of the polynomial.M =N-1
% tau = CGL_nodes(1)
% -1     1
% tau = CGL_nodes(2)
% -1     0     1
% tau = CGL_nodes(3)
% -1  -0.5   0.5   
syms n
    % calculate node locations
    for k = 0:N
        nodes(k+1,1) = -cos(pi*k/n); % symbolically to maintain precision
    end
    nodes = double(subs(nodes,'n',N));

    %% Calculation of weights of LG-Nodes
    N = length(nodes)-1;
    theta = pi*(0:N)'/N; 
    x = cos(theta);
    w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
    if mod(N,2)==0
    w(1) = 1/(N^2-1); w(N+1) = w(1);
    for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    v = v - cos(N*theta(ii))/(N^2-1);
    else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    end
    w(ii) = 2*v/N;
    
    weights = flipud(w)'; % since tau from [-1,1] not [1 -1]

    % %% Calculation of D_matrix of LG-Nodes
    % D_LG = zeros(N, N);
    % 
    % % Diagonal elements
    % for i = 1:N
    %     D_LG(i, i) = nodes(i) / (1 - nodes(i) ^ 2);
    % end
    % 
    % % Non-diagonal elements
    % for i = 1:N
    %     for j = 1:N
    %         if i ~= j
    %             % Substituting symbolic values with numeric values
    %             node_i = nodes(i);
    %             node_j = nodes(j);
    %             D_LG(i, j) = subs(diff(legendreP(N, z), z), node_i) * node_i / ...
    %                          ((node_i - node_j) * subs(diff(legendreP(N, z), z), node_j) * node_j);
    %         end
    %     end
    % end

    %% Plotting the nodes
    close all;
    figure('Position', [200, 200, 400, 100]);
    x = nodes;
    y = zeros(size(x)); % y-coordinate will be zero for all nodes
    plot(x, y, '*', 'MarkerSize', 7, 'MarkerEdgeColor', 'r','LineWidth', 1.5); % Plotting nodes with grey circles
    xlabel('Nodes(\tau)','FontWeight', 'bold');
    % ylabel('LGL','FontWeight', 'bold');
    title(['Chebyshev-Gauss-Lobatto Nodes for N = ', num2str(N)]);
    grid on;
    set(gca, 'YTickLabel', []);
    ylim([-0.1, 0.1]);
end
