function Objective = single_stage_objective_func(x,M)
h_f = x(M);        
Objective = -h_f;
end