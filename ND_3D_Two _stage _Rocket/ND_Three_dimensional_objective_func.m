function Objective = ND_Three_dimensional_objective_func(x,M,problem)


m0 = problem.m0;
m0_2 = problem.m0_2;
Re = problem.Re;
mu = problem.mu;
g0 = problem.g0;

n_length = 1/Re;
n_velocity = 1/sqrt(mu/Re);
n_time = n_length/n_velocity;
n_mass = 1/m0;
n_thrust = 1/(m0*g0);

% mp1 = m0 - (x(7*M));
mp2 = m0*n_mass - x(21*M);
% F1 = mp1/m0;
F = mp2/(m0*n_mass);
Objective = F;
end