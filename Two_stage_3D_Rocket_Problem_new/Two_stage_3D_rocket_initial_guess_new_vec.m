function [x0] = Two_stage_3D_rocket_initial_guess_new_vec(M,problem)

m0 = problem.m0;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
m0_1 = problem.m0_2;
mp1 = problem.mp1;
mass1_f = problem.mass1_f;
m0_2 = problem.m0_2;
Re = problem.Re;
mu = problem.mu;
mass2_i = problem.mass2_i;
mass2_f = problem.mass2_f;
g0 = problem.g0;
M1 = problem.M1;
M2 = problem.M2;


lat_i = problem.lat_i;
long_i = problem.long_i;
Elev_i = problem.Elev_i;
Azim_i = problem.Azim_i;
hi = problem.hi;
Vi = problem.Vi;


hf_f = problem.hf_f;
Vf_f = problem.Vf_f;

latitude_i = lat_i;
longitude_i = long_i;


Omega_z = 2*pi/(24*60*60);
t0 = 0;
stage_time = 137;
final_time = 1650;
addpath('../PS_methods')              % add the PS_method file directory
N = M-1;                              % Order of the polynomial
[nodes,~] = LGL_nodes(N);             % calculate scaled node locations and weights
  
t_1 = ((stage_time-t0)/2).*nodes(1:M1)+(stage_time+t0)/2;
t_2 = ((final_time-stage_time)/2).*nodes(M1+1:M)+(final_time+stage_time)/2;

% Rx1_I = (Re+10)*cos(deg2rad(28))*cos(deg2rad(0));
% Ry1_I = (Re+10)*cos(deg2rad(28))*sin(deg2rad(0));
% Rz1_I = (Re+10)*sin(deg2rad(28));
% 
% Rx2_I = (Re+hf_f)*cos(deg2rad(1))*cos(deg2rad(87));
% Ry2_I = (Re+hf_f)*cos(deg2rad(1))*sin(deg2rad(87));
% Rz2_I = (Re+hf_f)*sin(deg2rad(1));
% 
% Rx_I = interp1([0;1650], [Rx1_I;Rx2_I], [t_1;t_2], 'spline');
% Ry_I = interp1([0;1650], [Ry1_I;Ry2_I], [t_1;t_2], 'spline');
% Rz_I = interp1([0;1650], [Rz1_I;Rz2_I], [t_1;t_2], 'spline');

[R_I,V_I,T_I,qn] = lat_long_elev_azi_vec(M,problem);
Rx_I = R_I(1,1:M);
Ry_I = R_I(2,1:M);
Rz_I = R_I(3,1:M);

% Vx_I = V_I(1,1:M);
% Vy_I = V_I(2,1:M);
% Vz_I = V_I(3,1:M);

Thrustx_I = T_I(1,1:M);
Thrusty_I = T_I(2,1:M);
Thrustz_I = T_I(3,1:M);

q1 = qn(1,1:M);
q2 = qn(2,1:M);
q3 = qn(3,1:M);
q4 = qn(4,1:M);
% 
Vx1_I = 10.*cos(deg2rad(28));
Vy1_I = Omega_z*(Re+hi)*cos(deg2rad(28)); 
Vz1_I = 10.*sin(deg2rad(28));

Vx2_I = sqrt(mu/(Re+hf_f))*cos(deg2rad(1))*cos(deg2rad(87));
Vy2_I = sqrt(mu/(Re+hf_f))*cos(deg2rad(1))*sin(deg2rad(87));
Vz2_I = sqrt(mu/(Re+hf_f))*sin(deg2rad(1));

Vx_I = interp1([0;1650], [Vx1_I;Vx2_I], [t_1;t_2], 'spline');
Vy_I = interp1([0;1650], [Vy1_I;Vy2_I], [t_1;t_2], 'spline');
Vz_I = interp1([0;1650], [Vz1_I;Vz2_I], [t_1;t_2], 'spline');

x0(0*M+1:1*M) = Rx_I;
x0(1*M+1:2*M) = Ry_I;
x0(2*M+1:3*M) = Rz_I;
x0(3*M+1:4*M) = Vx_I;
x0(4*M+1:5*M) = Vy_I;
x0(5*M+1:6*M) = Vz_I;
x0(6*M+1:6*M+M1) = linspace(m0,mass1_f,M1); 
x0(6*M+M1+1:7*M) = linspace(mass2_i,mass2_f,M2);
x0(7*M+1:7*M+M1) = linspace(Thrust_max,mass1_f*g0*1.2,M1);
x0(7*M+M1+1:8*M) = linspace(Thrust_max_2,mass2_f*g0*1.2,M2);
x0(8*M+1:9*M) = 1;
x0(9*M+1:10*M) = 0;
x0(10*M+1:11*M) = 0;
x0(11*M+1:12*M) = q1;
x0(12*M+1:13*M) = q2;
x0(13*M+1:14*M) = q3;
x0(14*M+1:15*M) = q4;
x0(15*M+1) = 135;
x0(15*M+2) = 1650;

end