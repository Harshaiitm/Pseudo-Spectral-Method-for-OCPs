function [c,ceq,dc,dceq]= Two_stage_3D_rocket_Nonlinear_func_LGL(x,M,D,problem)

dc = [];
dceq = [];

Re = problem.Re;
h_scale = problem.h_scale;
mu = problem.mu;
Omega_x = problem.Omega_x;
Omega_y = problem.Omega_y; 
Omega_z = problem.Omega_z;
rho0 = problem.rho0;
g0 = problem.g0;

M1 = problem.M1;
M2 = problem.M2;
nodes = problem.nodes;

m0_1 = problem.m0_1;
m0_2 = problem.m0_2;
m0 = problem.m0;
mass1_i = problem.mass1_i;
mass1_f = problem.mass1_f;
mass2_i = problem.mass2_i;
mass2_f = problem.mass2_f;

A_ref = problem.A_ref;
Cbx = problem.Cbx;
dCby_by_dbeta = problem.dCby_by_dbeta;
dCbz_by_dalpha = problem.dCbz_by_dalpha;

Isp = problem.Isp;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;

q_max = problem.q_max;
a_sen_max = problem.a_sen_max;

hi = problem.hi;
Vi = problem.Vi;
t0 = problem.t0;

hf_f = problem.hf_f;
Vf_f = problem.Vf_f;
gamma_f = problem.gamma_f;
inclin_f = problem.inclin_f;

Gamma = problem.Gamma;
R = problem.R;
Temp0 = problem.Temp0;

% Decision veriables
Rx = x(0*M+1:1*M);
Ry = x(1*M+1:2*M);
Rz = x(2*M+1:3*M);
Vx = x(3*M+1:4*M);
Vy = x(4*M+1:5*M);
Vz = x(5*M+1:6*M);
mass = x(6*M+1:7*M);
Thrust = x(7*M+1:8*M);
uTx = x(8*M+1:9*M);
uTy = x(9*M+1:10*M);
uTz = x(10*M+1:11*M);
q1 = x(11*M+1:12*M);
q2 = x(12*M+1:13*M);
q3 = x(13*M+1:14*M);
q4 = x(14*M+1:15*M);
stage_time = x(15*M+1);
final_time = x(15*M+2);

% Attitude matrix for stage_1
Q11 = q1.^2 - q2.^2 - q3.^2 + q4.^2;
Q12 = 2*(q1.*q2 + q3.*q4);
Q13 = 2*(q1.*q3 - q2.*q4);
Q21 = 2*(q1.*q2 - q3.*q4);
Q22 = -q1.^2 + q2.^2 - q3.^2 + q4.^2;
Q23 = 2*(q2.*q3 + q1.*q4);
Q31 = 2*(q1.*q3 + q2.*q4);
Q32 = 2*(q2.*q3 - q1.*q4);
Q33 = -q1.^2 - q2.^2 + q3.^2 + q4.^2;


Q = [Q11 Q12 Q13; Q21 Q22 Q23; Q31 Q32 Q33];


% Inertial Thrust vector
Thrust_x = linspace(Thrust_max,mass2_f*g0*1.2,M).*(Q11.*uTx + Q12.*uTy + Q13.*uTz);
Thrust_y = linspace(Thrust_max,mass2_f*g0*1.2,M).*(Q21.*uTx + Q22.*uTy + Q23.*uTz);
Thrust_z = linspace(Thrust_max,mass2_f*g0*1.2,M).*(Q31.*uTx + Q32.*uTy + Q33.*uTz);

% Gravity
g_x = (-mu*Rx)./(Rx.^2 + Ry.^2 + Rz.^2).^(3/2);
g_y = (-mu*Ry)./(Rx.^2 + Ry.^2 + Rz.^2).^(3/2);
g_z = (-mu*Rz)./(Rx.^2 + Ry.^2 + Rz.^2).^(3/2);

Vrel_x = Vx - Rz.*Omega_y + Ry.*Omega_z;
Vrel_y = Vy - Rx.*Omega_z + Rz.*Omega_x;
Vrel_z = Vz - Ry.*Omega_x + Rx.*Omega_y;

Vrel = sqrt(Vrel_x.^2 + Vrel_y.^2 + Vrel_z.^2);

R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
h = R-Re;

rho = rho0*exp(-(h)./h_scale);

Vbx = Q11.*Vrel_x + Q21.*Vrel_y + Q31.*Vrel_z;
Vby = Q12.*Vrel_x + Q22.*Vrel_y + Q32.*Vrel_z;
Vbz = Q13.*Vrel_x + Q23.*Vrel_y + Q33.*Vrel_z;

Vb = sqrt(Vbx.^2 + Vby.^2 + Vbz.^2);

q_mag = 0.5* rho.* Vb.^2;

% Inertial Aerodynamic Coefficients
alpha = atan(Vbz./Vbx);
beta = atan(Vby./sqrt(Vbx.^2 + Vbz.^2));
phi = atan(Vby./Vbx);


Cby = dCby_by_dbeta*beta;
Cbz = dCbz_by_dalpha*alpha;

Cx = Q11.*Cbx + Q12.*Cby + Q13.*Cbz;
Cy = Q21.*Cbx + Q22.*Cby + Q23.*Cbz;
Cz = Q31.*Cbx + Q32.*Cby + Q33.*Cbz;

A_x = q_mag.* Cx * A_ref;
A_y = q_mag.* Cy * A_ref;
A_z = q_mag.* Cz * A_ref;

ceq = zeros(8*M+18,1);
% system dynamics
ceq(1:M,1) = D*Rx' - ((final_time-t0)/2)*(Vx)';
ceq(M+1:2*M,1) = D*Ry' - ((final_time-t0)/2)*(Vy)';
ceq(2*M+1:3*M,1) = D*Rz' - ((final_time-t0)/2)*(Vz)';
ceq(3*M+1:4*M,1) = D*Vx' - ((final_time-t0)/2)*((Thrust_x + A_x)./mass + g_x)';
ceq(4*M+1:5*M,1) = D*Vy' - ((final_time-t0)/2)*((Thrust_y + A_y)./mass + g_y)';
ceq(5*M+1:6*M,1) = D*Vz' - ((final_time-t0)/2)*((Thrust_z + A_z)./mass + g_z)';
ceq(6*M+1:7*M,1) = D*mass' + ((final_time-t0)/2)*((Thrust)./(g0.*Isp))';

Rx_I = (Re+10)*cos(deg2rad(28))*cos(deg2rad(0));
Ry_I = (Re+10)*cos(deg2rad(28))*sin(deg2rad(0));
Rz_I = (Re+10)*sin(deg2rad(28));

Vx_I = 10.*cos(deg2rad(28));
Vy_I = Omega_z*(Re+hi)*cos(deg2rad(28)); 
Vz_I = 10.*sin(deg2rad(28));

% Initial and final constraints
ceq(7*M+1) = Rx(1) - Rx_I;
ceq(7*M+2) = Vx(1) - Vx_I;
ceq(7*M+3) = Ry(1) - Ry_I;
ceq(7*M+4) = Vy(1) - Vy_I;
ceq(7*M+5) = Rz(1) - Rz_I;
ceq(7*M+6) = Vz(1) - Vz_I;

% time span
t_1= ((stage_time-t0)/2).*nodes(1:M1)+(stage_time+t0)/2;
t_2= ((final_time-stage_time)/2).*nodes(M1+1:M)+(final_time+stage_time)/2;

% Knotting constraints
ceq(7*M+6) = t_2(1) - t_1(end);
ceq(7*M+7) = Rx(M1+1) - Rx(M1);
ceq(7*M+8) = Ry(M1+1) - Ry(M1);
ceq(7*M+9) = Rz(M1+1) - Rz(M1);
ceq(7*M+10) = Vx(M1+1) - Vx(M1);
ceq(7*M+11) = Vy(M1+1) - Vy(M1);
ceq(7*M+12) = Vz(M1+1) - Vz(M1);
ceq(7*M+13) = mass(1) - mass1_i;
ceq(7*M+14) = mass(M1+1) - mass2_i;

ceq(7*M+15) = (Rx(end)^2 + Ry(end)^2 + Rz(end)^2) - ((Re + hf_f)^2);
ceq(7*M+16) = (Vx(end)^2 + Vy(end)^2 + Vz(end)^2) - (Vf_f^2);
ceq(7*M+17) = (Rx(end)*Vx(end) + Ry(end)*Vy(end) + Rz(end)*Vz(end)) - ((Re + hf_f) * Vf_f * sin(gamma_f));
ceq(7*M+18) = (Rx(end)*Vy(end) - Ry(end)*Vx(end)) - ((Re + hf_f) * Vf_f * cos(gamma_f) * sin(inclin_f));

% Normalisation constraint for Quaternion elements 
ceq(7*M+19:8*M+18) = (q1.^2 + q2.^2 + q3.^2 + q4.^2) - 1;

% Senced acceleration calculation
a_sen_x = (Thrust_x + A_x)./mass;
a_sen_y = (Thrust_y + A_y)./mass;
a_sen_z = (Thrust_z + A_z)./mass;
a_sen_mag = sqrt((a_sen_x).^2 + (a_sen_y).^2 + (a_sen_z).^2);

% Inequality_constraints
% c = [];
c = zeros(4*M,1);
% 
c(1:M,1) = q_mag - q_max;
c(1*M+1:2*M,1) = a_sen_mag.^2 - a_sen_max^2;
c(2*M+1:3*M,1) = Thrust - Thrust_max;

%Mass change constraint
for i = 1:(M-1)
    c(3*M+i) = (mass(i+1) - mass(i)) - 0;
end
if any(isnan(c) | isinf(c))
   % pause;

   inf_indices = find(isnan(c) | isinf(c));
   disp(inf_indices);
end

% ceq_inf_indices = find(isnan(ceq) | isinf(ceq));
end