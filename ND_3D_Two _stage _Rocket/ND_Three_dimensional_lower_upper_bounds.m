function [lb, ub] = ND_Three_dimensional_lower_upper_bounds(M,problem)

m0 = problem.m0;
g0 = problem.g0;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
Thrust_max_3 = problem.Thrust_max_3;
m0_2 = problem.m0_2;
Re = problem.Re;
hi = problem.hi;
Vi = problem.Vi;
mu = problem.mu;

n_length = 1/Re;
n_velocity = 1/sqrt(mu/Re);
n_time = n_length/n_velocity;
n_mass = 1/m0;
n_thrust = 1/(m0*g0);

lb(0*M+1:1*M) = -inf;
lb(1*M+1:2*M) = -inf;
lb(2*M+1:3*M) = -inf;
lb(3*M+1:4*M) = -inf;
lb(4*M+1:5*M) = -inf;
lb(5*M+1:6*M) = -inf;
lb(6*M+1:7*M) = m0_2*n_mass;
lb(7*M+1:8*M) = 0;
lb(8*M+1:9*M) = 0;
lb(9*M+1:10*M) = 0;
lb(10*M+1:11*M) = -1;
lb(11*M+1:12*M) = -1;
lb(12*M+1:13*M) = -1;
lb(13*M+1:14*M) = -1;
lb(14*M+1:15*M) = -inf;
lb(15*M+1:16*M) = -inf;
lb(16*M+1:17*M) = -inf;
lb(17*M+1:18*M) = -inf;
lb(18*M+1:19*M) = -inf;
lb(19*M+1:20*M) = -inf;
lb(20*M+1:21*M) = 0*n_mass;
lb(21*M+1:22*M) = 0;
lb(22*M+1:13*M) = 0;
lb(23*M+1:24*M) = 0;
lb(24*M+1:25*M) = -1;
lb(25*M+1:26*M) = -1;
lb(26*M+1:27*M) = -1;
lb(27*M+1:28*M) = -1;
lb(28*M+1) = 0*n_time;
lb(28*M+2) = 0*n_time;

ub(0*M+1:1*M) = inf;
ub(1*M+1:2*M) = inf;
ub(2*M+1:3*M) = inf;
ub(3*M+1:4*M) = inf;
ub(4*M+1:5*M) = inf;
ub(5*M+1:6*M) = inf;
ub(6*M+1:7*M) = m0*n_mass;
ub(7*M+1:8*M) = Thrust_max*n_thrust;
ub(8*M+1:9*M) = Thrust_max*n_thrust;
ub(9*M+1:10*M) = Thrust_max*n_thrust;
ub(10*M+1:11*M) = 1;
ub(11*M+1:12*M) = 1;
ub(12*M+1:13*M) = 1;
ub(13*M+1:14*M) = 1;
ub(14*M+1:15*M) = inf;
ub(15*M+1:16*M) = inf;
ub(16*M+1:17*M) = inf;
ub(17*M+1:18*M) = inf;
ub(18*M+1:19*M) = inf;
ub(19*M+1:20*M) = inf;
ub(20*M+1:21*M) = m0_2*n_mass;
ub(21*M+1:22*M) = Thrust_max_2*n_thrust;
ub(22*M+1:23*M) = Thrust_max_2*n_thrust;
ub(23*M+1:24*M) = Thrust_max_2*n_thrust;
ub(24*M+1:25*M) = 1;
ub(25*M+1:26*M) = 1;
ub(26*M+1:27*M) = 1;
ub(27*M+1:28*M) = 1;
ub(28*M+1) = 150*n_time;
ub(28*M+2) = 1750*n_time;
end