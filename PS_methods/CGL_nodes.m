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
% -1  -0.5   0.5   1

function tau = CGL_nodes(M)
    syms n
    % calculate node locations
    for k = 0:M
        tau(k+1,1) = -cos(pi*k/n); % symbolically to maintain precision
    end
    tau = double(subs(tau,'n',M));
end