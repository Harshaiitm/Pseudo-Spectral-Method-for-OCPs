function [x0] = Three_dimensional_initial_guess_new_vec(M,problem)

m0 = problem.m0;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
m0_1 = problem.m0_2;
mp1 = problem.mp1;
mass1_f = problem.mass1_f;
m0_2 = problem.m0_2;
Re = problem.Re;
mu = problem.mu;
mass2_f = problem.mass2_f;

lat_i = problem.lat_i;
long_i = problem.long_i;
Elev_i = problem.Elev_i;
Azim_i = problem.Azim_i;
hi = problem.hi;
Vi = problem.Vi;

lat_s = problem.lat_s;
long_s = problem.long_s;
Elev_s = problem.Elev_s;
Azim_s = problem.Azim_s;
hf_s = problem.hf_s;
Vf_s = problem.Vf_s;

lat_f = problem.lat_f;
long_f = problem.long_f;
Elev_f = problem.Elev_f;
Azim_f = problem.Azim_f;
hf_f = problem.hf_f;
Vf_f = problem.Vf_f;

[R1_I,V1_I,T1_I,qn1] = lat_long_elev_azi_vec1(M,problem);
Rx0_1 = R1_I(1,1:M);
Ry0_1 = R1_I(2,1:M);
Rz0_1 = R1_I(3,1:M);
 
Vx0_1 = V1_I(1,1:M);
Vy0_1 = V1_I(2,1:M);
Vz0_1 = V1_I(3,1:M);

Thrustx0_1 = T1_I(1,1:M);
Thrusty0_1 = T1_I(2,1:M);
Thrustz0_1 = T1_I(3,1:M);

uTx01 = 1;
uTy01 = 0;
uTz01 = 0;

q11 = qn1(1,1:M);
q12 = qn1(2,1:M);
q13 = qn1(3,1:M);
q14 = qn1(4,1:M);

[R2_I,V2_I,T2_I,qn2] = lat_long_elev_azi_vec2(M,problem); 

Rx0_2 = R2_I(1,1:M);
Ry0_2 = R2_I(2,1:M);
Rz0_2 = R2_I(3,1:M);
 
Vx0_2 = V2_I(1,1:M);
Vy0_2 = V2_I(2,1:M);
Vz0_2 = V2_I(3,1:M);

Thrustx0_2 = T2_I(1,1:M);
Thrusty0_2 = T2_I(2,1:M);
Thrustz0_2 = T2_I(3,1:M);

uTx02 = 1;
uTy02 = 0;
uTz02 = 0;

q21 = qn2(1,1:M);
q22 = qn2(2,1:M);
q23 = qn2(3,1:M);
q24 = qn2(4,1:M);

Omega_z = 2*pi/(24*60*60);
t0 = 0;
stage_time = 137;
final_time = 1650;
addpath('../PS_methods')                    % add the PS_method file directory
N = M-1;                            % Order of the polynomial
[nodes,~] = LGL_nodes(N);     % calculate scaled node locations and weights
  
t_1 = ((stage_time-t0)/2).*nodes+(stage_time+t0)/2;
t_2 = ((final_time-stage_time)/2).*nodes+(final_time+stage_time)/2;
Vx1_I = 10.*cos(deg2rad(28));
Vy1_I = Omega_z*(Re+hi)*cos(deg2rad(28)); 
Vz1_I = 10.*sin(deg2rad(28));

Vx2_I = sqrt(mu/(Re+hf_f));
Vy2_I = 0;
Vz2_I = 0;

Vx_I = interp1([0;1650], [Vx1_I;Vx2_I], [t_1;t_2], 'spline');
Vy_I = interp1([0;1650], [Vy1_I;Vy2_I], [t_1;t_2], 'spline');
Vz_I = interp1([0;1650], [Vz1_I;Vz2_I], [t_1;t_2], 'spline');


x0(0*M+1:1*M) = Rx0_1;                          % Rx_1
x0(1*M+1:2*M) = Ry0_1;                          % Ry_1
x0(2*M+1:3*M) = Rz0_1;                          % Rz_1
x0(3*M+1:4*M) = Vx_I(1:M);                          % Vx_1                                           
x0(4*M+1:5*M) = Vy_I(1:M);                          % Vy_1
x0(5*M+1:6*M) = Vz_I(1:M);                          % Vz_1    
x0(6*M+1:7*M) = linspace(m0,m0-mass1_f,M);                       % mass_1
x0(7*M+1:8*M) = linspace(Thrust_max,Thrust_max_2,M);                % Thrust_x1                                
x0(8*M+1:9*M) = uTx01;
x0(9*M+1:10*M) = uTy01;
x0(10*M+1:11*M) = uTz01;              % Thrust_z1
x0(11*M+1:12*M) = q11;                                          % q11
x0(12*M+1:13*M) = q12;                                          % q12
x0(13*M+1:14*M) = q13;                                          % q13
x0(14*M+1:15*M) = q14;                                          % q14
x0(15*M+1:16*M) = Rx0_2;                        % Rx_2
x0(16*M+1:17*M) = Ry0_2;                        % Ry_2 
x0(17*M+1:18*M) = Rz0_2;                        % Rz_2 
x0(18*M+1:19*M) = Vx_I(M+1:2*M);                        % Vx_2 
x0(19*M+1:20*M) = Vy_I(M+1:2*M);                        % Vy_2
x0(20*M+1:21*M) = Vz_I(M+1:2*M);                        % Vz_2
x0(21*M+1:22*M) = linspace(m0_2,mass2_f,M);                     % mass_2
x0(22*M+1:23*M) = linspace(Thrust_max_2,10,M);              % Thrust_x2
x0(23*M+1:24*M) = uTx02;
x0(24*M+1:25*M) = uTy02;
x0(25*M+1:26*M) = uTz02;            % Thrust_z2
x0(26*M+1:27*M) = q21;                                          % q21
x0(27*M+1:28*M) = q22;                                          % q22
x0(28*M+1:29*M) = q23;                                          % q23
x0(29*M+1:30*M) = q24;                                          % q24
x0(30*M+1) = 137;                                               % stage_time
x0(30*M+2) = 1750;                                              % final_time


end