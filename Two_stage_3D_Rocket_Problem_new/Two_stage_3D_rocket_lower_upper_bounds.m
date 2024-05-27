function [lb, ub] = Two_stage_3D_rocket_lower_upper_bounds(M,problem)

m0 = problem.m0;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
Thrust_max_3 = problem.Thrust_max_3;
m0_2 = problem.m0_2;
Re = problem.Re;
hi = problem.hi;
Vi = problem.Vi;
mu = problem.mu;
mass1_f = problem.mass1_f;
mass2_f = problem.mass2_f;
g0 = problem.g0;
M1 = problem.M1;
M2 = problem.M2;

% Lower Bounds 
lb(0*M+1:1*M) = -inf;                      % Rx_1
lb(1*M+1:2*M) = -inf;                  % Ry_1
lb(2*M+1:3*M) = -inf;              % Rz_1
lb(3*M+1:4*M) = -inf;              % Vx_1
lb(4*M+1:5*M) = -inf;              % Vy_1
lb(5*M+1:6*M) = -inf;              % Vz_1
lb(6*M+1:6*M+M1) = mass1_f;           % mass_1
lb(6*M+M1+1:7*M) = mass2_f;           % mass_2
lb(7*M+1:7*M+M1) = mass1_f*g0*1.2;    % Thrust_1
lb(7*M+M1+1:8*M) = mass2_f*g0*1.2;    % Thrust_2
lb(8*M+1:9*M) = -1;                   % uTx1
lb(9*M+1:10*M) = -1;                % uTy1
lb(10*M+1:11*M) = -1;              % uTz1
lb(11*M+1:12*M) = -1;              % q11
lb(12*M+1:13*M) = -1;              % q12
lb(13*M+1:14*M) = -1;              % q13
lb(14*M+1:15*M) = -1;              % q14
lb(15*M+1) = 0;                       % stage_time
lb(15*M+2) = 0;                       % final_time

% Lower Bounds 
ub(0*M+1:1*M) = inf;                       % Rx_1
ub(1*M+1:2*M) = inf;                   % Ry_1
ub(2*M+1:3*M) = inf;               % Rz_1
ub(3*M+1:4*M) = inf;               % Vx_1
ub(4*M+1:5*M) = inf;               % Vy_1
ub(5*M+1:6*M) = inf;               % Vz_1
ub(6*M+1:6*M+M1) = m0;                % mass_1
ub(6*M+M1+1:7*M) = m0_2;              % mass_2
ub(7*M+1:7*M+M1) = m0*g0*1.2;         % Thrust_1
ub(7*M+M1+1:8*M) = m0_2*g0*1.2;       % Thrust_2
ub(8*M+1:9*M) = 1;                 % uTx1
ub(9*M+1:10*M) = 1;                 % uTy1
ub(10*M+1:11*M) = 1;               % uTz1
ub(11*M+1:12*M) = 1;               % q11
ub(12*M+1:13*M) = 1;               % q12
ub(13*M+1:14*M) = 1;               % q13
ub(14*M+1:15*M) = 1;               % q14
ub(15*M+1) = inf;                     % stage_time
ub(15*M+2) = inf;                     % final_time

end