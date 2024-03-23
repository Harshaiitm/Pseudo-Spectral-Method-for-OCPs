function Objective = Normal_Tangential_Polar_Rocket_objective_func(x,M,m0)
% mass(end) = x(5*M);
mp = m0 - x(5*M);
F = mp/m0;        
Objective = F;
end