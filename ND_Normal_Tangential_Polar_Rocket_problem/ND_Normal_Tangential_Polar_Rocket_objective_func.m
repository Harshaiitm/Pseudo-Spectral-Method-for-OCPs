function Objective =ND_Normal_Tangential_Polar_Rocket_objective_func(x,N,m0,n_mass)
% mass(end) = x(5*N+5);
mp = m0*n_mass-x(5*N+5);
F = mp/(m0*n_mass);        
Objective = -F;
end