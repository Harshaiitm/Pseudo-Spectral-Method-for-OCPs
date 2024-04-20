function [x0] = Three_dimensional_initial_guess(M,problem)

m0 = problem.m0;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
m0_2 = problem.m0_2;
Re = problem.Re;
mu = problem.mu;
mass2_f = problem.mass2_f;

lat_i = problem.lat_i;
long_i = problem.long_i;
Rx_i = problem.Rx_i;
Ry_i = problem.Ry_i;
Rz_i = problem.Rz_i;
Vx_i = problem.Vx_i;
Vy_i = problem.Vy_i;
Vz_i = problem.Vz_i;
hi = problem.hi;
Vi = problem.Vi;
Thrust_xi = problem.Thrust_xi; 
Thrust_yi = problem.Thrust_yi;
Thrust_zi = problem.Thrust_zi;

lat_s = problem.lat_s;
long_s = problem.long_s;
Rx_s = problem.Rx_s;
Ry_s = problem.Ry_s;
Rz_s = problem.Rz_s;
Vx_s = problem.Vx_s;
Vy_s = problem.Vy_s;
Vz_s = problem.Vz_s;
hf_s = problem.hf_s;
Vf_s = problem.Vf_s;
Thrust_xs = problem.Thrust_xs; 
Thrust_ys = problem.Thrust_ys;
Thrust_zs = problem.Thrust_zs;

lat_f = problem.lat_f;
long_f = problem.long_f;
Rx_f = problem.Rx_f;
Ry_f = problem.Ry_f;
Rz_f = problem.Rz_f;
Vx_f = problem.Vx_f;
Vy_f = problem.Vy_f;
Vz_f = problem.Vz_f;
hf = problem.hf;
Vf = problem.Vf;
Thrust_xf = problem.Thrust_xf; 
Thrust_yf = problem.Thrust_yf;
Thrust_zf = problem.Thrust_zf;


% 
x0(0*M+1:1*M) = linspace(Rx_i,Rx_s,M);                          % Rx_1
x0(1*M+1:2*M) = linspace(Ry_i,Ry_s,M);                          % Ry_1
x0(2*M+1:3*M) = linspace(Rz_i,Rz_s,M);                          % Rz_1
x0(3*M+1:4*M) = linspace(Vx_i,Vx_s,M);                          % Vx_1                                           
x0(4*M+1:5*M) = linspace(Vy_i,Vy_s,M);                          % Vy_1
x0(5*M+1:6*M) = linspace(Vz_i,Vz_s,M);                          % Vz_1    
x0(6*M+1:7*M)  = linspace(m0,m0_2,M);                           % mass_1
x0(7*M+1:8*M) = linspace(Thrust_xi,Thrust_xs,M);                % Thrust_x1                                
x0(8*M+1:9*M) = linspace(Thrust_yi,Thrust_ys,M);                % Thrust_y1
x0(9*M+1:10*M) = linspace(Thrust_zi,Thrust_zs,M);               % Thrust_z1
x0(10*M+1:11*M) = 0.5;                                          % q11
x0(11*M+1:12*M) = 0.5;                                          % q12
x0(12*M+1:13*M) = 0.5;                                          % q13
x0(13*M+1:14*M) = 0.5;                                          % q14
x0(14*M+1:15*M) = linspace(Rx_s,Rx_f,M);                        % Rx_2
x0(15*M+1:16*M) = linspace(Ry_s,Ry_f,M);                        % Ry_2 
x0(16*M+1:17*M) = linspace(Rz_s,Rz_f,M);                        % Rz_2 
x0(17*M+1:18*M) = linspace(Vx_s,Vx_f,M);                        % Vx_2 
x0(18*M+1:19*M) = linspace(Vy_s,Vx_f,M);                        % Vy_2
x0(19*M+1:20*M) = linspace(Vz_s,Vx_f,M);                        % Vz_2
x0(20*M+1:21*M) = linspace(m0_2,mass2_f,M);                     % mass_2
x0(21*M+1:22*M) = linspace(Thrust_xs,Thrust_xf,M);              % Thrust_x2
x0(22*M+1:23*M) = linspace(Thrust_ys,Thrust_yf,M);              % Thrust_y2
x0(23*M+1:24*M) = linspace(Thrust_zs,Thrust_zf,M);              % Thrust_z2
x0(24*M+1:25*M) = 0.5;                                          % q21
x0(25*M+1:26*M) = 0.5;                                          % q22
x0(26*M+1:27*M) = 0.5;                                          % q23
x0(27*M+1:28*M) = 0.5;                                          % q24
x0(28*M+1) = 137;                                               % stage_time
x0(28*M+2) = 1750;                                              % final_time

end