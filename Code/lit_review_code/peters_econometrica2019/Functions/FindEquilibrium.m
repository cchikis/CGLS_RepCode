% This function finds the equilibrium innovation outcomes for given
% parameters



function EqOutcomes = FindEquilibrium(param)


% Step 1: Find optimal tau and I 

[I,tau,fval] = SolveEquilibrium_I_tau(param);

% Step 2: Find optimal x

x = (param.phix/param.phiz * param.zeta^(-1))^(1/(param.zeta-1));

% Step 3: Get Entry rate

z = tau - x;

% Save results

EqOutcomes.x = x;
EqOutcomes.tau = tau;
EqOutcomes.z = z;
EqOutcomes.I = I;
EqOutcomes.fval = fval;

end

function [I,tau,fval] = SolveEquilibrium_I_tau(param)


MinProblem = @(y) EquilibriumConditions_I_tau(y,param);

% Set starting value 
y_init(1) = 0.5; %I
y_init(2) = 0.2; %tau

options = optimoptions('fmincon','Display','off');
[y_sol, fval] = fmincon(MinProblem, y_init, [], [], [], [], [0,0], [], [], options);

I = y_sol(1);
tau = y_sol(2);


end


% This function calculates the difference in the equilibrium conditions,
% which we are going to set to zero

function diff = EquilibriumConditions_I_tau(y,param)

I = y(1);
tau  = y(2);

s = tau/I; %The destruction intensity

%% Parameters

l = param.lambda;
zeta = param.zeta;
phiI = param.phiI;
phiz = param.phiz;
rho = param.rho;

%% EQUATION 1
Lambda = s / (l * s + l - 1);
h_shift = (zeta-1)/zeta * (param.phix/(param.phiz^(zeta) * zeta))^(1/(zeta-1));


H = Lambda * ( l/(l-1) * zeta/phiI * I^(zeta-1) * (rho + I * s) + zeta/phiI * I^(zeta) ) + ...
        1/phiz * s * I - h_shift;
    
%% EQUATION 2

G = zeta / phiI * I^(zeta-1) + (rho + I * s)^(-1) * ( (zeta-1)/phiI * I^(zeta) +  h_shift );

%% EQUILIBRIUM SYSTEM

eps1 = param.L - H;
eps2 = 1/phiz - G;

%% Get MSE

diff = eps1^(2) + eps2^(2);



end

