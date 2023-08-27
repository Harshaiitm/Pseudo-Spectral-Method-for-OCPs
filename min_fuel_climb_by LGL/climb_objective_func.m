function Objective = climb_objective_func(x,N)
mass = x(3*N+4:4*N+4);         
Objective = -mass(end);
end