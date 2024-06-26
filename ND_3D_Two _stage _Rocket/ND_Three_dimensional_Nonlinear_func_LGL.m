function [c,ceq,dc,dceq]= ND_Three_dimensional_Nonlinear_func_LGL(x,M,D,problem)

dc = [];
dceq = [];

n_length = problem.n_length;
n_velocity = problem.n_velocity;
n_mass = problem.n_mass;
n_thrust = problem.n_thrust;
n_time = problem.n_time;

Re = problem.Re;
h_scale = problem.h_scale;
mu = problem.mu;
Omega_x = problem.Omega_x;
Omega_y = problem.Omega_y; 
Omega_z = problem.Omega_z;
rho0 = problem.rho0;
g0 = problem.g0;


m0_1 = problem.m0_1;
m0_2 = problem.m0_2;
m0 = problem.m0;
mass1_i = problem.mass1_i;
mass1_i = mass1_i*n_mass;
mass1_f = problem.mass1_f;
mass1_f = mass1_f*n_mass;
mass2_i = problem.mass2_i;
mass2_i = mass2_i*n_mass;
mass2_f = problem.mass2_f;

A_ref = problem.A_ref;
Cbx1 = problem.Cbx1;
Cbx2 = problem.Cbx2;
dCby1_by_dbeta = problem.dCby1_by_dbeta;
dCby2_by_dbeta = problem.dCby2_by_dbeta;
dCbz1_by_dalpha = problem.dCbz1_by_dalpha;
dCbz2_by_dalpha = problem.dCbz2_by_dalpha;


Isp = problem.Isp;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
Thrust_max_3 = problem.Thrust_max_3;


q_max = problem.q_max;
a_sen_max = problem.a_sen_max;

hi = problem.hi;
Vi = problem.Vi;
t0 = problem.t0;
t0 = t0*n_length;
hf_f = problem.hf_f;
Vf_f = problem.Vf_f;
gamma_f = problem.gamma_f;
inclin_f = problem.inclin_f;

% Decision veriables
Rx_1 = x(0*M+1:1*M);
Ry_1 = x(1*M+1:2*M);
Rz_1 = x(2*M+1:3*M);
Vx_1 = x(3*M+1:4*M);
Vy_1 = x(4*M+1:5*M);
Vz_1 = x(5*M+1:6*M);
mass_1 = x(6*M+1:7*M);
Thrust_1 = x(7*M+1:8*M);
uTx1 = x(8*M+1:9*M);
uTy1 = x(9*M+1:10*M);
uTz1 = x(10*M+1:11*M);
q11 = x(11*M+1:12*M);
q12 = x(12*M+1:13*M);
q13 = x(13*M+1:14*M);
q14 = x(14*M+1:15*M);
Rx_2 = x(15*M+1:16*M);
Ry_2 = x(16*M+1:17*M);
Rz_2 = x(17*M+1:18*M);
Vx_2 = x(18*M+1:19*M);
Vy_2 = x(19*M+1:20*M);
Vz_2 = x(20*M+1:21*M);
mass_2 = x(21*M+1:22*M);
Thrust_2 = x(22*M+1:23*M);
uTx2 = x(23*M+1:24*M);
uTy2 = x(24*M+1:25*M);
uTz2 = x(25*M+1:26*M);
q21 = x(26*M+1:27*M);
q22 = x(27*M+1:28*M);
q23 = x(28*M+1:29*M);
q24 = x(29*M+1:30*M);
stage_time = x(30*M+1);
final_time = x(30*M+2);

% Attitude matrix for stage_1
Q111 = q11.^2 - q12.^2 - q13.^2 + q14.^2;
Q112 = 2*(q11.*q12 + q13.*q14);
Q113 = 2*(q11.*q13 - q12.*q14);
Q121 = 2*(q11.*q12 - q13.*q14);
Q122 = -q11.^2 + q12.^2 - q13.^2 + q14.^2;
Q123 = 2*(q12.*q13 + q11.*q14);
Q131 = 2*(q11.*q13 + q12.*q14);
Q132 = 2*(q12.*q13 - q11.*q14);
Q133 = -q11.^2 - q12.^2 + q13.^2 + q14.^2;


Q1 = [Q111 Q112 Q113; Q121 Q122 Q123; Q131 Q132 Q133];

% Attitude matrix for stage_2
Q211 = q21.^2 - q22.^2 - q23.^2 + q24.^2;
Q212 = 2*(q21.*q22 + q23.*q24);
Q213 = 2*(q21.*q23 - q22.*q24);
Q221 = 2*(q21.*q22 - q23.*q24);
Q222 = -q21.^2 + q22.^2 - q23.^2 + q24.^2;
Q223 = 2*(q22.*q23 + q21.*q24);
Q231 = 2*(q21.*q23 + q22.*q24);
Q232 = 2*(q22.*q23 - q21.*q24);
Q233 = -q21.^2 - q22.^2 + q23.^2 + q24.^2;


Q2 = [Q211 Q212 Q213; Q221 Q222 Q223; Q231 Q232 Q233];

% Dimensionlization
Rx_1 = Rx_1/n_length;
Ry_1 = Ry_1/n_length;
Rz_1 = Rz_1/n_length;
Rx_2 = Rx_2/n_length;
Ry_2 = Ry_2/n_length;
Rz_2 = Rz_2/n_length;
Vx_1 = Vx_1/n_velocity;
Vy_1 = Vy_1/n_velocity;
Vz_1 = Vz_1/n_velocity;
Vx_2 = Vx_2/n_velocity;
Vy_2 = Vy_2/n_velocity;
Vz_2 = Vz_2/n_velocity;
mass_1 = mass_1/n_mass;
mass_2 = mass_2/n_mass;
Thrust_1 = Thrust_1/n_thrust;
Thrust_2 = Thrust_2/n_thrust;
stage_time = stage_time/n_time;
final_time = final_time/n_time;

Thrust_x1 = linspace(Thrust_max,Thrust_max_2,M).*(Q111.*uTx1 + Q112.*uTy1 + Q113.*uTz1);
Thrust_y1 = linspace(Thrust_max,Thrust_max_2,M).*(Q121.*uTx1 + Q122.*uTy1 + Q123.*uTz1);
Thrust_z1 = linspace(Thrust_max,Thrust_max_2,M).*(Q131.*uTx1 + Q132.*uTy1 + Q133.*uTz1);

Thrust_x2 = linspace(Thrust_max_2,Thrust_max_3,M).*(Q211.*uTx2 + Q212.*uTy2 + Q213.*uTz2);
Thrust_y2 = linspace(Thrust_max_2,Thrust_max_3,M).*(Q221.*uTx2 + Q222.*uTy2 + Q223.*uTz2);
Thrust_z2 = linspace(Thrust_max_2,Thrust_max_3,M).*(Q231.*uTx2 + Q232.*uTy2 + Q233.*uTz2);

% Gravity
g_x1 = (-mu*Rx_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_y1 = (-mu*Ry_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_z1 = (-mu*Rz_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_x2 = (-mu*Rx_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);
g_y2 = (-mu*Ry_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);
g_z2 = (-mu*Rz_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);

Vrel_x1 = Vx_1 - Rz_1.*Omega_y + Ry_1.*Omega_z; 
Vrel_y1 = Vy_1 - Rx_1.*Omega_z + Rz_1.*Omega_x;
Vrel_z1 = Vz_1 - Ry_1.*Omega_x + Rx_1.*Omega_y;
Vrel_x2 = Vx_2 - Rz_2.*Omega_y + Ry_2.*Omega_z;
Vrel_y2 = Vy_2 - Rx_2.*Omega_z + Rz_2.*Omega_x;
Vrel_z2 = Vz_2 - Ry_2.*Omega_x + Rx_2.*Omega_y;

Vrel_1 = sqrt(Vrel_x1.^2 + Vrel_y1.^2 + Vrel_z1.^2);
Vrel_2 = sqrt(Vrel_x2.^2 + Vrel_y2.^2 + Vrel_z2.^2);

Vbx1 = Q111.*Vrel_x1 + Q121.*Vrel_y1 + Q131.*Vrel_z1;
Vby1 = Q112.*Vrel_x1 + Q122.*Vrel_y1 + Q132.*Vrel_z1;
Vbz1 = Q113.*Vrel_x1 + Q123.*Vrel_y1 + Q133.*Vrel_z1;
Vbx2 = Q211.*Vrel_x2 + Q221.*Vrel_y2 + Q231.*Vrel_z2;
Vby2 = Q212.*Vrel_x2 + Q222.*Vrel_y2 + Q232.*Vrel_z2;
Vbz2 = Q213.*Vrel_x2 + Q223.*Vrel_y2 + Q233.*Vrel_z2;

Vb_1 = sqrt(Vbx1.^2 + Vby1.^2 + Vbz1.^2);
Vb_2 = sqrt(Vbx2.^2 + Vby2.^2 + Vbz2.^2);

R_1 = sqrt(Rx_1.^2 + Ry_1.^2 + Rz_1.^2);
h_1 = R_1-Re;
R_2 = sqrt(Rx_2.^2 + Ry_2.^2 + Rz_2.^2);
h_2 = R_2-Re;

rho_1 = rho0*exp(-(R_1-Re)./h_scale);
rho_2 = rho0*exp(-(R_2-Re)./h_scale);

q_mag1 = 0.5* rho_1.* Vb_1.^2;
q_mag2 = 0.5* rho_2.* Vb_2.^2;
q_mag = [q_mag1 q_mag2];


% Inertial Aerodynamic Coefficients
alpha_1 = deg2rad(atan(Vbz1./Vbx1));
alpha_2 = deg2rad(atan(Vbz2./Vbx2));
beta_1 = deg2rad(atan(Vby1./sqrt(Vbx1.^2 + Vbz1.^2)));
beta_2 = deg2rad(atan(Vby2./sqrt(Vbx2.^2 + Vbz2.^2)));
phi_1 = deg2rad(atan(Vby1./Vbx1));
phi_2 = deg2rad(atan(Vby2./Vbx2));

% % Cbx1 = -Ca;
% Cby1 = -Cn*sin(phi_1);
% Cbz1 = -Cn*cos(phi_1);
% % Cbx2 = -Ca;
% Cby2 = -Cn*sin(phi_2);
% Cbz2 = -Cn*cos(phi_2);

Cby1 = dCby1_by_dbeta*beta_1;
Cbz1 = dCbz1_by_dalpha*alpha_1;
Cby2 = dCby2_by_dbeta*beta_2;
Cbz2 = dCbz2_by_dalpha*alpha_2;

Cbx1 = problem.Cbx1;
Cbx2 = problem.Cbx2;
Cby1 = problem.Cby1;
Cby2 = problem.Cby2;
Cbz1 = problem.Cbz1;
Cbz2 = problem.Cbz2;

Cx1 = Q111.*Cbx1 + Q112.*Cby1 + Q113.*Cbz1;
Cy1 = Q121.*Cbx1 + Q122.*Cby1 + Q123.*Cbz1;
Cz1 = Q131.*Cbx1 + Q132.*Cby1 + Q133.*Cbz1;
Cx2 = Q211.*Cbx2 + Q212.*Cby2 + Q213.*Cbz2;
Cy2 = Q221.*Cbx2 + Q222.*Cby2 + Q223.*Cbz2;
Cz2 = Q231.*Cbx2 + Q232.*Cby2 + Q233.*Cbz2;


A_x1 = q_mag1.* Cx1 * A_ref;
A_y1 = q_mag1.* Cy1 * A_ref;
A_z1 = q_mag1.* Cz1 * A_ref;
A_x2 = q_mag2.* Cx2 * A_ref;
A_y2 = q_mag2.* Cy2 * A_ref;
A_z2 = q_mag2.* Cz2 * A_ref;

Rx_1 = Rx_1*n_length;
Ry_1 = Ry_1*n_length;
Rz_1 = Rz_1*n_length;
Rx_2 = Rx_2*n_length;
Ry_2 = Ry_2*n_length;
Rz_2 = Rz_2*n_length;
Vx_1 = Vx_1*n_velocity;
Vy_1 = Vy_1*n_velocity;
Vz_1 = Vz_1*n_velocity;
Vx_2 = Vx_2*n_velocity;
Vy_2 = Vy_2*n_velocity;
Vz_2 = Vz_2*n_velocity;
mass_1 = mass_1*n_mass;
mass_2 = mass_2*n_mass;
Thrust_1 = Thrust_1*n_thrust;
Thrust_2 = Thrust_2*n_thrust;
Thrust_x1 = Thrust_x1*n_thrust;
Thrust_y1 = Thrust_y1*n_thrust;
Thrust_z1 = Thrust_z1*n_thrust;
Thrust_x2 = Thrust_x2*n_thrust;
Thrust_y2 = Thrust_y2*n_thrust;
Thrust_z2 = Thrust_z2*n_thrust;
A_x1 = A_x1*n_thrust;
A_x2 = A_x2*n_thrust;
A_y1 = A_y1*n_thrust;
A_y2 = A_y2*n_thrust;
A_z1 = A_z1*n_thrust;
A_z2 = A_z2*n_thrust;
stage_time = stage_time*n_time;
final_time = final_time*n_time;
t0 = t0*n_time;
Isp = Isp*n_time;

ceq = zeros(19*M,1);
% system dynamics
ceq(1:M,1) = D*Rx_1' - ((stage_time-t0)/2)*(Vx_1)';
ceq(M+1:2*M,1) = D*Ry_1' - ((stage_time-t0)/2)*(Vy_1)';
ceq(2*M+1:3*M,1) = D*Rz_1' - ((stage_time-t0)/2)*(Vz_1)';
ceq(3*M+1:4*M,1) = D*Vx_1' - ((stage_time-t0)/2)*((Thrust_x1 + A_x1)./mass_1 + 1./Rx_1.^2)';
ceq(4*M+1:5*M,1) = D*Vy_1' - ((stage_time-t0)/2)*((Thrust_y1 + A_y1)./mass_1 + 1./Ry_1.^2)';
ceq(5*M+1:6*M,1) = D*Vz_1' - ((stage_time-t0)/2)*((Thrust_z1 + A_z1)./mass_1 + 1./Rz_1.^2)';
ceq(6*M+1:7*M,1) = D*mass_1' + ((stage_time-t0)/2)*((Thrust_1)./(g0.*Isp))';
ceq(7*M+1:8*M,1) = D*Rx_2' - ((final_time-stage_time)/2)*(Vx_2)';
ceq(8*M+1:9*M,1) = D*Ry_2' - ((final_time-stage_time)/2)*(Vy_2)';
ceq(9*M+1:10*M,1) = D*Rz_2' - ((final_time-stage_time)/2)*(Vz_2)';
ceq(10*M+1:11*M,1) = D*Vx_2' - ((final_time-stage_time)/2)*((Thrust_x2 + A_x2)./mass_2 + 1./Rx_2.^2 )';
ceq(11*M+1:12*M,1) = D*Vy_2' - ((final_time-stage_time)/2)*((Thrust_y2 + A_y2)./mass_2 + 1./Ry_2.^2)';
ceq(12*M+1:13*M,1) = D*Vz_2' - ((final_time-stage_time)/2)*((Thrust_z2 + A_z2)./mass_2 + 1./Rz_2.^2)';
ceq(13*M+1:14*M,1) = D*mass_2' + ((final_time-stage_time)/2)*((Thrust_2)./(g0.*Isp))';
% 
[R1_I,V1_I,~,~] = lat_long_elev_azi_vec1(M,problem);
Rx0_1 = R1_I(1,1:M)*n_length;
Ry0_1 = R1_I(2,1:M)*n_length;
Rz0_1 = R1_I(3,1:M)*n_length;

Vx0_1 = V1_I(1,1:M)*n_velocity;
Vy0_1 = V1_I(2,1:M)*n_velocity;
Vz0_1 = V1_I(3,1:M)*n_velocity;

% % Initial and final constraints
ceq(14*M+1) = Rx_1(1) - Rx0_1(1);
ceq(14*M+2) = Vx_1(1) - Vx0_1(1);
ceq(14*M+3) = Ry_1(1) - Ry0_1(1);
ceq(14*M+4) = Vy_1(1) - Vy0_1(1);
ceq(14*M+5) = Rz_1(1) - Rz0_1(1);
ceq(14*M+6) = Vz_1(1) - Vz0_1(1);
ceq(14*M+7) = (Rx_1(1)^2 + Ry_1(1)^2 + Rz_1(1)^2) - ((Re*n_length + hi*n_length)^2);
ceq(14*M+8) = (Vx_1(1)^2 + Vy_1(1)^2 + Vz_1(1)^2) - (Vi*n_velocity)^2;


% Knotting constraints
ceq(14*M+11) = Rx_2(1) - Rx_1(end);
ceq(14*M+12) = Ry_2(1) - Ry_1(end);
ceq(14*M+13) = Rz_2(1) - Rz_1(end);
ceq(14*M+14) = Vx_2(1) - Vx_1(end);
ceq(14*M+15) = Vy_2(1) - Vy_1(end);
ceq(14*M+16) = Vz_2(1) - Vz_1(end);
ceq(14*M+17) = mass_1(1) - mass1_i;
ceq(14*M+18) = mass_1(end) - mass1_f;
ceq(14*M+19) = mass_2(1) - mass2_i;

ceq(14*M+20) = (Rx_2(end)^2 + Ry_2(end)^2 + Rz_2(end)^2) - ((Re*n_length + hf_f*n_length)^2);
ceq(14*M+21) = (Vx_2(end)^2 + Vy_2(end)^2 + Vz_2(end)^2) - (Vf_f*n_velocity^2);
ceq(14*M+22) = (Rx_2(end)*Vx_2(end) + Ry_2(end)*Vy_2(end) + Rz_2(end)*Vz_2(end)) - ((Re*n_length + hf_f*n_length) * Vf_f*n_velocity * sin(gamma_f));
ceq(14*M+23) = (Rx_2(end)*Vy_2(end) - Ry_2(end)*Vx_2(end)) - ((Re*n_length + hf_f*n_length) * Vf_f*n_velocity * cos(gamma_f) * sin(inclin_f));


% Normalisation constraint for Quaternion elements 
ceq(17*M+1:18*M) = (q11.^2 + q12.^2 + q13.^2 + q14.^2) - 1;
ceq(18*M+1:19*M) = (q21.^2 + q22.^2 + q23.^2 + q24.^2) - 1; 
% find(isnan(ceq))
% find(ceq~=0)



% Senced acceleration calculation
a_sen_x1 = D*(Vx_1/n_velocity)' - (g_x1/g0)';
a_sen_y1 = D*(Vy_1/n_velocity)' - (g_y1/g0)';
a_sen_z1 = D*(Vz_1/n_velocity)' - (g_z1/g0)';
a_sen_mag1 = sqrt((a_sen_x1).^2 + (a_sen_y1).^2 + (a_sen_z1).^2);

a_sen_x2 = D*(Vx_2/n_velocity)' - (g_x2/g0)';
a_sen_y2 = D*(Vy_2/n_velocity)' - (g_y2/g0)';
a_sen_z2 = D*(Vz_2/n_velocity)' - (g_z2/g0)';
a_sen_mag2 = sqrt((a_sen_x2).^2 + (a_sen_y2).^2 + (a_sen_z2).^2);


% Inequality_constraints
% c = [];
c = zeros(7*M,1);

q_mag1 = 0.5* rho_1.* (Vb_1).^2;
q_mag2 = 0.5* rho_2.* (Vb_2).^2;
q_mag = [q_mag1 q_mag2];

c(1:M,1) = q_mag1 - q_max;
c(M+1:2*M,1) = q_mag2 - q_max;

c(2*M+1:3*M,1) = a_sen_mag1.^2 - a_sen_max^2;
c(3*M+1:4*M,1) = a_sen_mag2.^2 - a_sen_max^2;
% 
c(4*M+1:5*M,1) = Thrust_1 - Thrust_max;
c(5*M+1:6*M,1) = Thrust_2 - Thrust_max_2;

mass_1 = mass_1/n_mass;
mass_2 = mass_2/n_mass;
%Mass change constraint
for i = 1:(M-1)
    c(6*M+i) = (mass_1(i+1) - mass_1(i)) - 0;
end
for i = 1:(M-1)  
    c(6*M+(M-1)+i) = (mass_2(i+1) - mass_2(i)) - 0;
end


end

