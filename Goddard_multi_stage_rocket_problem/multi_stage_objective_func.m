function Objective = multi_stage_objective_func(x,N)
h_f = x(5*N+5);        
Objective = -h_f;
end