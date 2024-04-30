% Objective
fun = @(x) norm([x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
                 x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))]);

% Bounds
lb = [-4;-4];
ub = [4;4];

% Initial Solution
x0 = [-4;-4];

% Build optiprob structure (intermediate structure)
prob = optiprob('fun',fun,'bounds',lb,ub,'x0',x0);

% Setup solver options
opts1 = optiset('solver','nomad');
opts2 = optiset('solver','pswarm');
opts3 = optiset('solver','nlopt','solverOpts',nloptset('algorithm','GN_DIRECT'));
opts4 = optiset('solver','ipopt','warnings','off');

% Build OPTI objects
Opt1 = opti(prob,opts1);
Opt2 = opti(prob,opts2);
Opt3 = opti(prob,opts3);
Opt4 = opti(prob,opts4);

% Solve the GNLP problems
[x1,fval1] = solve(Opt1,x0);
[x2,fval2] = solve(Opt2,x0);
[x3,fval3] = solve(Opt3,x0);
[x4,fval4] = solve(Opt4,x0);