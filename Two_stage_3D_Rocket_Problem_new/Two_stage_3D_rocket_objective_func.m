function Objective = Two_stage_3D_rocket_objective_func(x,M,problem)

m0 = problem.m0;
m0_2 = problem.m0_2;

mp2 = m0 - x(7*M);
F = mp2/m0;
Objective = F;

end