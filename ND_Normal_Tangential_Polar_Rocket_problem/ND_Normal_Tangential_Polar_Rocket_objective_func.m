function Objective =ND_Normal_Tangential_Polar_Rocket_objective_func(x,M,m0,n_mass)
% mass(end) = x(5*M);
mp = m0*n_mass - x(5*M);
F = mp/(m0*n_mass);        
Objective = -F;
end
