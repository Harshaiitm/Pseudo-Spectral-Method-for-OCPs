function [nodes,weights,D_LGL] = LGL_computations(N)
%% Calculation of the LGL-Nodes 
% These nodes are defined by the roots of (1-x^2)L_N'(x)
% where L_N'(x) is the derivative of Nth order Legendre Function

syms z;
expr = legendreP(N,z);         % Legendre Polynomial of order N

%  finds the roots of the equation (1-x^2)P_N'(x) = 0, 
% The vpasolve function is used to solve the equation symbolically.
nodes = double(vpasolve((1-z^2)*(diff(expr,z,1))));

%% Calculation of weights of LGL-Nodes
weights = zeros(N+1,1);
for i = 1:N+1
    weights(i,1) = ((2/(N*(N+1)))*((legendreP(N,nodes(i)))^-2));
end
%% Calculation of D_matrix of LGL-Nodes
    D_LGL = zeros(N+1, N+1);

    % Diagonal elements
    D_LGL(N+1, N+1) = N*(N+1)/4;
    D_LGL(1, 1) = -N*(N+1)/4;

    % Non-diagonal elements
    for i = 1:N+1
        for j = 1:N+1
            if i ~= j
                D_LGL(i, j) = (legendreP(N,nodes(i)))/((legendreP(N,nodes(j)))*(nodes(i)-nodes(j)));
            end
        end
    end

  % elements corresponding to i = j and not equal to 1 or N
    for k = 2:N-1
        D_LGL(k, k) = 0;
    end

%% Plotting the nodes
close all;
figure('Position', [200, 200, 400, 100]);
x = nodes;
y = zeros(size(x)); % y-coordinate will be zero for all nodes
plot(x, y, '*', 'MarkerSize', 7, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5, 0.5, 0.5]); % Plotting nodes with grey circles
xlabel('Nodes(\tau)','FontWeight', 'bold');
% ylabel('LGL','FontWeight', 'bold');
title(['Legendre-Gauss-Lobatto Nodes for N = ', num2str(N)]);
grid on;
set(gca, 'YTickLabel', []);
ylim([-0.1, 0.1]);
end
