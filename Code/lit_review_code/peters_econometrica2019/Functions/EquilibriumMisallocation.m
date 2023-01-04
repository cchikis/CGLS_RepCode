% This function calculates varies statistics about the equlibrium
% distribution of markups

function m = EquilibriumMisallocation(EqObjects,param)

% Pareto tail
theta = log(1 + EqObjects.tau/EqObjects.I)/log(param.lambda);

% Average markup
Emu = theta/(theta-1);

% std dev of log markups across products
sig_mu = theta^(-1);

% TFP misallocation (M)
M = exp(-1/theta)*(1+theta)/theta;

% Factor price misallocation (Lambda)
L = theta/(1+theta);


% std dev of log markups across firms
NumProducts = 15000;
[Gap, Names, ~] = SimulateCrossSection(EqObjects.I, EqObjects.z, EqObjects.tau, NumProducts);
mu_product = param.lambda.^Gap;                     %product level markups

tmp = grpstats(mu_product.^(-1),Names,'mean');      %firm level markups
mu_firm = tmp.^(-1);

sig_mu_firm = var(log(mu_firm))^(1/2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save results
m = [theta, Emu, sig_mu, sig_mu_firm, M, L];
end