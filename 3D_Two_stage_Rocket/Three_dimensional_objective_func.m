function Objective = Three_dimensional_objective_func(x,M,problem)

m0 = problem.m0;
m0_2 = problem.m0_2;

% mp1 = m0 - (x(7*M));
mp2 = m0 - (x(22*M));
% F1 = mp1/m0;
F = mp2/m0;
Objective = F;
end