function [Clalpha, CD0, eta] = aero_data(Mach)

% Lookup Table 2: Aerodynamic Data 
MachLTAero = [0 0.4 0.8 0.9 1.0 1.2 1.4 1.6 1.8];
ClalphaLT = [3.44 3.44 3.44 3.58 4.44 3.44 3.01 2.86 2.44];
CD0LT = [0.013 0.013 0.013 0.014 0.031 0.041 0.039 0.036 0.035];
etaLT = [0.54 0.54 0.54 0.75 0.79 0.78 0.89 0.93 0.93];

% Fitting of Aerodynamic data with piecewise splines
Clalpha_pp = pchip(MachLTAero, ClalphaLT);
CD0_pp = pchip(MachLTAero, CD0LT);
eta_pp = pchip(MachLTAero, etaLT);

% piecewise polynomial at the specified Mach values
Clalpha = ppval(Clalpha_pp, Mach);
CD0 = ppval(CD0_pp, Mach);
eta = ppval(eta_pp, Mach);

end
