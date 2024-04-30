clear all; close all; clc;

N = 10;
nodes(1) = -1;
[nodes(2:N+1), weights] = LG_nodes(N, -1, 1);
nodes(N+2) = 1;
x = nodes;
syms x;
w = zeros(size(nodes));

for i = 1:length(nodes)
    xi = nodes(i);
    Li = 1;

    for j = 1:length(nodes)
        if j ~= i
            Li = Li * (x - nodes(j)) / (xi - nodes(j));
        end
    end
    
    % Compute the weight using numerical integration (quad function)
    w(i) = double(int(Li, -1, 1));
end

disp(w);
