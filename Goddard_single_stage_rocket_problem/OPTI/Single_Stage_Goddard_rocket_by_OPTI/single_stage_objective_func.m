function Objective = single_stage_objective_func(x,N)
h_f = x(N+1);        
Objective = -h_f;
end