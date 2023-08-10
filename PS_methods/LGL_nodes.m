function [nodes,weights] = LGL_nodes(N)
%% Calculation of the LGL-Nodes 
% These nodes are defined by the roots of (1-x^2)L_N'(x)
% where L_N'(x) is the derivative of Nth order Legendre Function

syms z;
expr = legendreP(N,z);         % Legendre Polynomial of order N

%  finds the roots of the equation (1-x^2)P_N'(x) = 0, 
% where P_N'(x) is the derivative of the Nth order Legendre polynomial. 
% The vpasolve function is used to solve the equation symbolically.
nodes = double(vpasolve((1-z^2)*(diff(expr,z,1))));

%% Calculation of weights of LGL-Nodes
weights = zeros(N+1,1);
for i = 1:N+1
 weights(i,1) = ((2/(N*(N+1)))*((legendreP(N,nodes(i)))^-2));
end

end