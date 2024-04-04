function Objective = Three_dimensional_objective_func(x,M,m0)
mp = m0 - (x(7*M)+x(17*M));
F = mp/m0;        
Objective = F;
end