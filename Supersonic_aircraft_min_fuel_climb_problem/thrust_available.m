function Thrust = thrust_available(h,Mach)
% Lookup Table 3: Propulsion Data (Thrust as function of mach number and altitude)
MachLT = [0; 0.2; 0.4; 0.6; 0.8; 1; 1.2; 1.4; 1.6; 1.8];
AltLT = 304.8 * [0 5 10 15 20 25 30 40 50 70];
ThrustLT = 4448.2 * [
    24.2 24.0 20.3 17.3 14.5 12.2 10.2 5.7 3.4 0.1;
    28.0 24.6 21.1 18.1 15.2 12.8 10.7 6.5 3.9 0.2;
    28.3 25.2 21.9 18.7 15.9 13.4 11.2 7.3 4.4 0.4;
    30.8 27.2 23.8 20.5 17.3 14.7 12.3 8.1 4.9 0.8;
    34.5 30.3 26.6 23.2 19.8 16.8 14.1 9.4 5.6 1.1;
    37.9 34.3 30.4 26.8 23.3 19.8 16.8 11.2 6.8 1.4;
    36.1 38.0 34.9 31.3 27.3 23.6 20.1 13.4 8.3 1.7;
    36.1 36.6 38.5 36.1 31.6 28.1 24.2 16.2 10.0 2.2;
    36.1 35.2 42.1 38.7 35.7 32.0 28.1 19.3 11.9 2.9;
    36.1 33.8 45.7 41.3 39.8 34.6 31.1 21.7 13.3 3.1];

% 2-D Interpolation
Thrust = interp2(AltLT,MachLT,ThrustLT,h,Mach,"spline");

end

