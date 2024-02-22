function [nodes,weights,D_LGR] = LGR_computations(N)
%% Calculation of the LGR-Nodes 
% These nodes are defined by the roots of (P_(N-1)(x)+P_(N)(x))
syms z;
expr_N = legendreP(N,z);         % Legendre Polynomial of order N
expr_n = legendreP(N-1,z);         % Legendre Polynomial of order N-1

%  finds the roots of the equation (P_(N-1)(x)+P_(N)(x))= 0, 
% The vpasolve function is used to solve the equation symbolically.
nodes = double(vpasolve(expr_n+expr_N));
% nodes(N+1) = 1;
% Calculation of weights of LGR-Nodes
weights = zeros(N,1);
for i = 1:N
    weights(i,1) = ((1-nodes(i))/((N+1)^2))*((legendreP(N,nodes(i)))^-2);
end

 %% Calculation of D_matrix of LGL-Nodes
    D_LGR = zeros(N, N);

    % Diagonal elements
    D_LGR(1, 1) = -N * (N + 2) / 4;
    for i = 2:N
        D_LGR(i, i) = 1 / (2 * (1 - nodes(i)));
    end

    % Non-diagonal elements
    for i = 1:N
        for j = 1:N
            if i ~= j
                D_LGR(i, j) = ((legendreP(N,nodes(i))*((1 -nodes(j))) / (legendreP(N,nodes(j)))*(1 - nodes(i)))*(nodes(i) - nodes(j)));
            end
        end
    end

%% Plotting the nodes
close all;
figure('Position', [200, 200, 400, 100]);
x = nodes;
y = zeros(size(x)); % y-coordinate will be zero for all nodes
plot(x, y, '*', 'MarkerSize', 7, 'MarkerEdgeColor','b','LineWidth',1.5); % Plotting nodes with grey circles
xlabel('Nodes(\tau)','FontWeight', 'bold');
% ylabel('LGL','FontWeight', 'bold');
title(['Legendre-Gauss-Radau Nodes for N = ', num2str(N)]);
grid on;
set(gca, 'YTickLabel', []);
ylim([-0.1, 0.1]);
end
