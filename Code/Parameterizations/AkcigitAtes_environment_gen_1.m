%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: AA_QCU_environment_gen_1.m
% Author: Craig A. Chikis
% Date: 12.07/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = AkcigitAtes_environment_gen_1()

	

	two_adjust = 1;
	s_max = 125;
	kappa = 9999;
	elastic_labor = false;
	algorithm_fsolve = 'trust-region-dogleg';
	LMS_bar_s = Inf;
	w_guess = .9;
	v_guess = repelem(1, length(-s_max:s_max));
	rho = power(1.05,1/12)-1;
	xbar = 1;
	eta = repelem(0.0278/12, length(0:s_max));
	eta(1) = 0;
	lab_tol = 1e-4;
	val_tol = 1e-7;
	gamma = 0.35;
	B = (7.179*gamma)^(-gamma);
	B = B/12;
	lambda = 1.044;
	entrance = true;
	subsidy = repelem(.05, length(-s_max:s_max));
	tax_rate = repelem(.3, length(-s_max:s_max));
	subsidyE = zeros(1,s_max+1);
	B_entrant = (.075*gamma)^(-gamma);
	B_entrant = B_entrant/12;
	PE = false;
	entry_type = "Ufuk/Sina";

	prob = zeros(2*s_max+1);
	prob_exog = zeros(s_max, 2*s_max+1);
	probE = zeros(s_max+1, 2*s_max+1);

	phi_wt = .0423;
	phi_wt_exog = 1;
	phi_wt_tilde = .0423;

	for (ii = 1:(s_max-1))
		prob(s_max+1-ii, s_max+1) = .0423;
		prob(s_max+1-ii, s_max+1-ii+1) = prob(s_max+1-ii, s_max+1-ii+1) + 1-.0423;

		prob(s_max+1+ii, s_max+1+ii+1) = 1;

		probE(s_max+1-ii, s_max+1) = .0423;
		probE(s_max+1-ii, s_max+1-ii+1) = probE(s_max+1-ii, s_max+1-ii+1) + 1-.0423;

		prob_exog(s_max+1-ii, s_max+1) = 1;
	end
	probE(s_max+1, s_max+2) = 1;
	probE(1, s_max+1) = .0423;
	probE(1, 2) = 1-.0423;

	prob_exog(1, s_max+1) = 1;

	prob(1, s_max+1) = .0423;
	prob(1, 2) = 1-.0423;

	prob(s_max+1, s_max+2) = 1;
	prob(2*s_max+1, 2*s_max+1) = 1;


	v_max_try = 1e10;


	model_str = "Akcigit_Ates_QCU";


	m.s_max = s_max;
	m.w_guess = w_guess;
	m.v_guess = v_guess;
	m.rho = rho;
	m.xbar = xbar;
	m.eta = eta;
	m.lab_tol = lab_tol;
	m.val_tol = val_tol; 
	m.B = B;
	m.gamma = gamma; 
	m.lambda = lambda;  
	m.entrance = entrance; 
	m.subsidy = subsidy; 
	m.tax_rate = tax_rate;
	m.subsidyE = subsidyE;
	m.B_entrant = B_entrant;
	m.PE = PE;
	m.entry_type = entry_type;
	m.prob = prob; 
	m.probE = probE; 
	m.prob_exog = prob_exog; 
	m.phi_wt = phi_wt;
	m.phi_wt_exog = phi_wt_exog;
	m.phi_wt_tilde = phi_wt_tilde; 
	% m.l = l;
	% m.l_tilde = 0;
	m.kappa = kappa; 
	m.LMS_bar_s = LMS_bar_s; 
	m.elastic_labor = elastic_labor; 
	m.algorithm_fsolve = algorithm_fsolve; 
	% m.alpha_fric = alpha_fric; 
	% m.bridge_method = bridge_method; 
	% m.EIS = EIS; 
 	% m.r_guess = r_guess; 
	% m.v_max_try = v_max_try; 
	% m.two_adjust = two_adjust; 

end