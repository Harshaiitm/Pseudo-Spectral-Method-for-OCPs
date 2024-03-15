function Objective = multi_stage_objective_func_OPTI(x,N)
h_f = x(5*N+5);        
Objective = -h_f;
end