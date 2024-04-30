function [nodes,weights] = LGR_node_t(N)
% Calculation of the LGL-Nodes 
% These nodes are defined by the roots of (1-x^2)L_N'(x)
% where L_N'(x) is the derivative of Nth order Legendre Function

syms z;
expr_N = legendreP(N,z);         % Legendre Polynomial of order N
expr_n = legendreP(N-1,z);         % Legendre Polynomial of order N

%  finds the roots of the equation (1-x^2)P_N'(x) = 0, 
% The vpasolve function is used to solve the equation symbolically.
nodes = double(vpasolve(expr_n+expr_N));
% nodes(N+1) = 1;
% Calculation of weights of LGL-Nodes
weights = zeros(N,1);
for i = 1:N
    weights(i,1) = ((1-nodes(i))/((N+1)^2))*((legendreP(N,nodes(i)))^-2);
end

% Plotting the nodes
close all;
figure('Position', [200, 200, 400, 100]);
x = nodes;
y = zeros(size(x)); % y-coordinate will be zero for all nodes
plot(x, y, '*', 'MarkerSize', 7, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5, 0.5, 0.5]); % Plotting nodes with grey circles
xlabel('Nodes(\tau)','FontWeight', 'bold');
% ylabel('LGL','FontWeight', 'bold');
title(['Legendre-Gauss-Radau Nodes for N = ', num2str(N)]);
grid on;
set(gca, 'YTickLabel', []);
ylim([-0.1, 0.1]);
end
