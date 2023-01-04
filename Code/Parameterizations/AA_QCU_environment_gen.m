%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: AA_QCU_environment_gen.m
% Author: Craig A. Chikis
% Date: 12/07/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = AA_QCU_environment_gen()




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
	eta = repelem(0, length(0:s_max));
	eta(1) = 0;
	lab_tol = 1e-4;
	val_tol = 1e-7;
	B = 0.1;
	gamma = 0.35;
	lambda = 1.05;
	entrance = false;
	subsidy = repelem(0, length(-s_max:s_max));
	tax_rate = repelem(0, length(-s_max:s_max));
	subsidyE = subsidy(1:(s_max+1));
	B_entrant = rand;
	PE = false;
	entry_type = "Undirected";
	phi_wt = 1;
	phi_wt_exog = 1;
	phi_wt_tilde = rand;
	l = 0;
	l_tilde = 0;
	l_exog = 0;
	l_opt = 0;

	[prob, probE, prob_exog] = F_mat_main(phi_wt, phi_wt_tilde, phi_wt_exog, l, l_tilde, s_max);


	v_max_try = 1e10;

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
	m.l = l;
	m.l_tilde = 0;
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