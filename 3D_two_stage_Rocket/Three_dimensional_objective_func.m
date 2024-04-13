function Objective = Three_dimensional_objective_func(x,M,m0,m0_2)
% mp1 = m0 - (x(7*M));
mp2 = m0 - (x(17*M));
% F1 = mp1/m0;
F = mp2/m0;
Objective = F;
end