function [x0] = Three_dimensional_initial_guess(M,problem)

m0 = problem.m0;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
m0_2 = problem.m0_2;
Re = problem.Re;
hi = problem.hi;
Vi = problem.Vi;
mu = problem.mu;
mass2_f = problem.mass2_f;

x0(0*M+1:1*M) = linspace((Re + hi) * cos(deg2rad(28)),6421000/sqrt(3),M);                          % Rx_1
x0(1*M+1:2*M) = linspace(0,6421000/sqrt(3),M);                                                % Ry_1

x0(2*M+1:3*M) = linspace((Re + hi) * sin(deg2rad(28)),6421000/sqrt(3),M);                    % Rz_1
x0(3*M+1:4*M) = linspace(Vi * cos(deg2rad(28)) ,sqrt(mu/(50000/sqrt(3))),M) ;                           % Vx_1
                                             % Vy_1
x0(4*M+1:5*M) = linspace(0,sqrt(mu/(50000/sqrt(3))),M) ;
x0(5*M+1:6*M) = linspace(Vi * sin(deg2rad(28)) ,sqrt(mu/(50000/sqrt(3))),M);                            % Vz_1
x0(6*M+1:7*M)  = linspace(m0,m0_2,M);                                                % mass_1
x0(7*M+1:8*M) = linspace(Thrust_max * cos(deg2rad(28)),Thrust_max_2,M);                      % Thrust_x1
                                                                                        % Thrust_y1
x0(8*M+1:9*M) = linspace(0,Thrust_max_2,M);
x0(9*M+1:10*M) = linspace(Thrust_max * sin(deg2rad(28)),Thrust_max_2,M);                     % Thrust_z1
x0(10*M+1:11*M) = 0.5;                 % q11
x0(11*M+1:12*M) = 0.5;                 % q12
x0(12*M+1:13*M) = 0.5;                 % q13
x0(13*M+1:14*M) = 0.5;                 % q14
x0(14*M+1:15*M) = linspace(6421000/sqrt(3),6771000/sqrt(3),M);                      % Rx_2
x0(15*M+1:16*M) = linspace(6421000/sqrt(3),6771000/sqrt(3),M);               % Ry_2 
x0(16*M+1:17*M) = linspace(6421000/sqrt(3),6771000/sqrt(3),M);               % Rz_2 
x0(17*M+1:18*M) = linspace(sqrt(mu/(50000*3)),sqrt(mu/(400000*3)),M);                      % Vx_2 
x0(18*M+1:19*M) = linspace(sqrt(mu/(50000*3)),sqrt(mu/(400000*3)),M);                      % Vy_2
x0(19*M+1:20*M) = linspace(sqrt(mu/(50000*3)),sqrt(mu/(400000*3)),M);                      % Vz_2
x0(20*M+1:21*M) = linspace(m0_2,mass2_f,M);                % mass_2
x0(21*M+1:22*M) = linspace(Thrust_max_2/sqrt(3),10,M);        % Thrust_x2
x0(22*M+1:23*M) = linspace(Thrust_max_2/sqrt(3),10,M);                   % Thrust_y2
x0(23*M+1:24*M) = linspace(Thrust_max_2/sqrt(3),10,M);        % Thrust_z2
x0(24*M+1:25*M) = 0.5;                 % q21
x0(25*M+1:26*M) = 0.5;                 % q22
x0(26*M+1:27*M) = 0.5;                 % q23
x0(27*M+1:28*M) = 0.5;                 % q24
x0(28*M+1) = 150;                      % stage_time
x0(28*M+2) = 1750;                    % final_time

end