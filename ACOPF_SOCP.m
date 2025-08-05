%% Solving the Optimal Power Flow by Second-Order Cone Programming

% Radial Distribution Network

% The bus shunt admittance, line charging susceptance, transformer ratio and phase shift are not considered.

clc
clear

define_constants; % Defines useful constants for indexing data

%% Small Distribution System Test Cases
% mpc = case4_dist; % There is no generator cost
% mpc = case10ba; % Infeasible model
% mpc = case12da; % 9.11427733e+00
% mpc = case15da; % 2.57639063e+01
% mpc = case15nbr; % 2.53601940e+01
% mpc = case16am; % Numerical trouble encountered
% mpc = case16ci; % Model is infeasible or unbounded % non-radial
% mpc = case17me; % Infeasible model
% mpc = case18; % 2.39922174e+02 % There is shunt admittance
% mpc = case18nbr; % 2.93821689e+01
% mpc = case22; % 1.36010776e+01
% mpc = case28da; % Model is infeasible or unbounded
mpc = case33bw; % 7.83535612e+01
% mpc = case33mg; % 7.76904401e+01
% mpc = case34sa; % 6.18102086e+01
% mpc = case38si; % 7.83535688e+01
% mpc = case51ga; % 5.18511293e+01
% mpc = case51he; % 3.91668410e+01
% mpc = case69; % 8.05418740e+01
% mpc = case70da; % Infeasible model % non-radial
% mpc = case74ds; % 1.35242752e+02
% mpc = case85; % Infeasible model
% mpc = case94pi; % Infeasible model
% mpc = case118zh; % Infeasible model
% mpc = case136ma; % Infeasible model
% mpc = case141; % 2.51546441e+02

%% Get Data

mpc0 = mpc;
mpc = ext2int(mpc);
mpc.bus(:, [GS, BS]) = 0; % The shunt admittance is set to 0
mpc.branch(:, BR_B) = 0; % The line charging susceptance is set to 0
mpc.branch(:, [TAP, SHIFT]) = 0; % The transformer off nominal turns ratio and phase shift angle are set to 0

baseMVA = mpc.baseMVA;
bus = mpc.bus;
branch = mpc.branch;
gen = mpc.gen;
gencost = mpc.gencost;

nb = size(bus, 1); %% number of buses
nl = size(branch, 1); %% number of branches
ng = size(gen, 1); %% number of dispatchable injections

Pd = bus(:, PD) / baseMVA; % real power demand (p.u.)
Qd = bus(:, QD) / baseMVA; % reactive power demand (p.u.)

Rl = branch(:, BR_R); % resistance (p.u.)
Xl = branch(:, BR_X); % reactance (p.u.)

Vmu = bus(:, VMAX); % maximum voltage magnitude (p.u.)
Vml = bus(:, VMIN); % minimum voltage magnitude (p.u.)

Pmax = gen(:, PMAX) / baseMVA; % maximum real power output (p.u.)
Pmin = gen(:, PMIN) / baseMVA; % minimum real power output (p.u.)

Qmax = gen(:, QMAX) / baseMVA; % maximum reactive power output (p.u.)
Qmin = gen(:, QMIN) / baseMVA; % minimum reactive power output (p.u.)

% find/prepare polynomial generator costs
Qpg = gencost(:, COST) * baseMVA^2; % quadratic
cpg = gencost(:, COST+1) * baseMVA; % linear
kpg = gencost(:, COST+2); % constant

% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
if ~isempty(il)
    flow_max = branch(il, RATE_A) / baseMVA;
end

% build connection matrices
Cg = sparse(gen(:, GEN_BUS), 1:ng, 1, nb, ng); %% connection matrix for generators & buses

Cf = sparse(branch(:, F_BUS), 1:nl, 1, nb, nl); %% connection matrix for line & from buses
Ct = sparse(branch(:, T_BUS), 1:nl, 1, nb, nl); %% connection matrix for line & to buses

%% Define variables
Pg = sdpvar(ng, 1, 'full'); % real power output (p.u.)
Qg = sdpvar(ng, 1, 'full'); % reactive power output (p.u.)
Pl = sdpvar(nl, 1, 'full'); % real power injected into "from" end of branch (p.u.)
Ql = sdpvar(nl, 1, 'full'); % reactive power injected into "from" end of branch (p.u.)
Vm2 = sdpvar(nb, 1, 'full'); % The square of bus voltage magnitude (p.u.)
Im2 = sdpvar(nl, 1, 'full'); % The square of branch current magnitude (p.u.)

%% Define constraints
con = [];

con = con + [Cg * Pg + Ct * (Pl - Rl .* Im2) - Cf * Pl == Pd]; %#ok<*NBRAK1> % real power balance
con = con + [Cg * Qg + Ct * (Ql - Xl .* Im2) - Cf * Ql == Qd]; % reactive power balance

con = con + [Ct' * Vm2 == Cf' * Vm2 - 2 * (Rl .* Pl + Xl .* Ql) + (Rl.^2 + Xl.^2) .* Im2]; % voltage constraint

con = con + cone([(Im2 + Cf' * Vm2)'; [2 * Pl'; 2 * Ql'; (Im2 - Cf' * Vm2)']]); % second-order cone relaxation
% (2*Pl)^2 + (2*Ql)^2 + (Im2-Cf'*Vm2)^2 <= (Im2+Cf'*Vm2)^2
% ==> ||2*Pl; 2*Ql; Im2-Cf'*Vm2||_2 <= Im2+Cf'*Vm2

con = con + [Vml.^2 <= Vm2 <= Vmu.^2]; %#ok<*CHAIN>
con = con + [Pmin <= Pg <= Pmax];
con = con + [Qmin <= Qg <= Qmax];

if ~isempty(il)
    con = con + cone([flow_max'; [Pl(il)'; Ql(il)']]); % flow limit
    % Pl^2 + Ql^2 <= flow_max^2
end

%% Define the objective function
obj = Pg' * (Qpg .* Pg) + cpg' * Pg + sum(kpg);

%% Solve

solver = 'gurobi';
% solver = 'fmincon'; % Solver of MATLAB
ops = sdpsettings('solver', solver, 'verbose', 2, 'savesolveroutput', 1);
% ops.gurobi.QCPDual = 1;
sol = solvesdp(con, obj, ops);

if sol.problem ~= 0
    disp(sol.info)
    return
end

Pl = value(Pl);
Ql = value(Ql);
Pg = value(Pg);
Qg = value(Qg);
Vm2 = value(Vm2);
Im2 = value(Im2);

%% Check if the second-order cone relaxation is tight
resid = (Im2 + Cf' * Vm2).^2 - (2 * Pl).^2 - (2 * Ql).^2 - (Im2 - Cf' * Vm2).^2;
fprintf(['\n\nThe second-order cone relaxation residual is:', num2str(max(resid)), '\n'])

%% Compare with MATPOWER

opt = mpoption('verbose', 0, 'out.all', 0);
res = runopf(mpc, opt);

if res.success == 0
    disp('MATPOWER runopf failed !!!')
    return
end

res = ext2int(res);

Pl_res = res.branch(:, PF) / baseMVA;
Ql_res = res.branch(:, QF) / baseMVA;

Pg_res = res.gen(:, PG) / baseMVA;
Qg_res = res.gen(:, QG) / baseMVA;

Vm_res = res.bus(:, VM);
Vm2_res = Vm_res.^2;

Im2_res = (Pl_res.^2 + Ql_res.^2) ./ (Cf' * Vm2_res);

fprintf(['\nThe error of the real power injected into "from" end of branch is:', num2str(norm(Pl_res-Pl)), '\n'])
fprintf(['\nThe error of the reactive power injected into "from" end of branch is:', num2str(norm(Ql_res-Ql)), '\n'])
fprintf(['\nThe error of the real power output is:', num2str(norm(Pg_res-Pg)), '\n'])
fprintf(['\nThe error of the reactive power output is:', num2str(norm(Qg_res-Qg)), '\n'])
fprintf(['\nThe error of the square of bus voltage magnitude is:', num2str(norm(Vm2_res-Vm2)), '\n'])
fprintf(['\nThe error of the square of branch current magnitude is:', num2str(norm(Im2_res-Im2)), '\n'])

