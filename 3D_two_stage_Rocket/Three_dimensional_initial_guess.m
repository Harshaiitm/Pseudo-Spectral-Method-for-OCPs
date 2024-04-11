function [x0] = Three_dimensional_initial_guess(M,problem)

m0 = problem.m0;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
m0_2 = problem.m0_2;
Re = problem.Re;
hi = problem.hi;
Vi = problem.Vi;

x0(1:M) = (Re + hi) * cos(deg2rad(28));                          % Rx_1
x0(M+1:2*M) = 0;                                                 % Ry_1
x0(2*M+1:3*M) = (Re + hi) * sin(deg2rad(28));                    % Rz_1
x0(3*M+1:4*M) = Vi * cos(deg2rad(28));                           % Vx_1
x0(4*M+1:5*M) = 0;                                               % Vy_1
x0(5*M+1:6*M) = Vi * sin(deg2rad(28));                           % Vz_1
x0(6*M+1:7*M)  = m0;                   % mass_1
x0(7*M+1:8*M) = Thrust_max;            % Thrust_x1
x0(8*M+1:9*M) = Thrust_max;            % Thrust_y1
x0(9*M+1:10*M) = Thrust_max;           % Thrust_z1
x0(10*M+1:11*M) = 50000;               % Rx_2
x0(11*M+1:12*M) = 50000;               % Ry_2 
x0(12*M+1:13*M) = 50000;               % Rz_2 
x0(13*M+1:14*M) = 400;                 % Vx_2 
x0(14*M+1:15*M) = 400;                 % Vy_2
x0(15*M+1:16*M) = 400;                 % Vz_2
x0(16*M+1:17*M) = m0_2;                % mass_2
x0(17*M+1:18*M) = Thrust_max_2;        % Thrust_x2
x0(18*M+1:19*M) = Thrust_max_2;        % Thrust_y2
x0(19*M+1:20*M) = Thrust_max_2;        % Thrust_z2
x0(20*M+1:21*M) = 0.5;                 % q1
x0(21*M+1:22*M) = 0.5;                 % q2
x0(22*M+1:23*M) = 0.5;                 % q3
x0(23*M+1:24*M) = 0.5;                 % q4
x0(24*M+1) = 100;                      % stage_time
x0(24*M+2) = 1000;                     % final_time

end