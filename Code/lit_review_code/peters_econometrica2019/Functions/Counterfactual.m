% This function calculates the equilibrium implications for the baseline
% economy and the US calibration


function [EI_Indonesia,EI_US] = Counterfactual(param)

% Calibrate structural parameters for US
NewMoments.EntryRate = 0.08;
NewMoments.LClnEmpl = log(2);
[paramUS, fval] = Calibration_US(NewMoments,param);

% Calculate equilibrium outcomes for both economies
EqOutcomes = FindEquilibrium(param);            % Baseline
EqOutcomesUS = FindEquilibrium(paramUS);        % US

% Equilibrium implications
EI_Indonesia = EquilibriumImplications(EqOutcomes, param);
EI_US = EquilibriumImplications(EqOutcomesUS, paramUS);

end


function EI = EquilibriumImplications(EqOutcomes, param)

tau = EqOutcomes.tau;
x = EqOutcomes.x;
I = EqOutcomes.I;
z = EqOutcomes.z;

% Pareto tail
theta = log(1 + tau/I)/log(param.lambda);

% Number of firms 
EI.NumFirms =  (tau - x)/x * log( tau/(tau - x));
% Output share of small firms 
EI.ShareSmall = 1 - (x/tau);
% Markup lifecycle
EI.LCMU = Avglnmu_Age(7.5, I, tau, x, param.lambda) - Avglnmu_Age(0.5, I, tau, x, param.lambda);
% Average markup
EI.Emu = theta/(theta-1);
% Std dev of log markups across products
EI.sig_mu = theta^(-1);
% Std dev of log markups across firms
NumProducts = 15000;
[Gap, Names, ~] = SimulateCrossSection(I, z, tau, NumProducts);
mu_product = param.lambda.^Gap;                     %product level markups
tmp = grpstats(mu_product.^(-1),Names,'mean');      %firm level markups
mu_firm = tmp.^(-1);
EI.sig_mu_firm = var(log(mu_firm))^(1/2);
% TFP misallocation (M)
EI.M = exp(-1/theta)*(1+theta)/theta;
% Factor price misallocation (Lambda)
EI.L = theta/(1+theta);
% Growth
EI.AggGrowth = log(param.lambda)*(I+tau);

end
