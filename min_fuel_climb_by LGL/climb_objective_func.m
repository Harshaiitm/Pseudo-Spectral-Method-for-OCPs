function Objective = climb_objective_func(x,N)
mass = x(4*N+5);         
Objective = -mass;
end