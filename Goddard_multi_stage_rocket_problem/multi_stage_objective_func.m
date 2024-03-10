function Objective = multi_stage_objective_func(x,M)
h_f = x(5*M);        
Objective = -h_f;
end