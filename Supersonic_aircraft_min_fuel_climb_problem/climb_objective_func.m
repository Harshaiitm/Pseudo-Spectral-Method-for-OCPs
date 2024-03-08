function Objective = climb_objective_func(x,M)
mass = x(4*M);         
Objective = -mass;
end