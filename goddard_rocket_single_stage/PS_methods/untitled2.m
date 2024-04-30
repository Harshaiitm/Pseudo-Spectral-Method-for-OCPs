N =10;
[nodes,weights]=LG_nodes(N,-1,1);

% Range for integration
a = 1;
b = length(nodes); % Assuming you want to integrate from 1 to the number of nodes

% Define the Lagrange polynomial
syms x;
lagrange_poly = 1;
for j = 1:length(nodes)
    if  j  ~= i
         lagrange_poly = lagrange_poly * (x - nodes(i)) / (nodes(1) - nodes(1i));
         w(i) = x * lagrange_poly;
         disp(w(i));
    end
end

% Define the function to be integrate



