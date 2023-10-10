function Objective = two_dim_rocket_single_stage_objective_func(x,N,m0)
% mass(end) = x(5*N+5);
mp = m0-x(5*N+5);
F = mp/m0;        
Objective = -F;
end
