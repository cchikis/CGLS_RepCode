%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: gamma33_current.m
% Author: Craig A. Chikis
% Date: 12/30/2020
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = gamma33_current()


	delta = 1;
	two_adjust = 1;
	s_max = 100;
	w_guess = 1;
	v_guess = ones(1,2*s_max+1);
	rho = power(1.02,1/12/delta)-1;
	xbar = 1;
	eta = zeros(1, s_max+1)/delta;
	lab_tol = 1e-4;
	val_tol = 5e-8;
	B = 0.145836674;
	gamma = 0.33;
	lambda = 1.016150266;  
	entrance = false;
	subsidy = zeros(1, 2*s_max+1);
	tax_rate = 0.2*ones(1, 2*s_max+1);  
	subsidyE = subsidy(1:(s_max+1));
	B_entrant = rand/delta;
	PE = false;
	entry_type = "NaN";
	model_str = "Benchmark";
	phi_wt = 0.137574116;
	phi_wt_exog = 1;
	phi_wt_tilde = rand;
	l = 0;
	l_tilde = 0;
	l_opt = 0;
	l_exog = 0;
	kappa = 9999;
	LMS_bar_s = Inf;
	elastic_labor = true;
	algorithm_fsolve = 'trust-region-dogleg';
	alpha_fric = Inf;
	bridge_method = false;
	EIS = 1;
	r_guess = power(1.03,1/12/delta)-1;
	winsor_vec = [1, 99];
	v_max_try = 1e10;
	zeta = 0.478052648/log(lambda); 

	[prob, probE, prob_exog] = F_mat_main(phi_wt, phi_wt_tilde, phi_wt_exog, l, l_tilde, s_max);


	N = 1100;
	T = 60*12; 
	winsor_vec = [1,99];
	fcitcount = 20;
	nbucket = 5;
	nbucket_patent = 5; 
	winsor_vec_drop = [-1, Inf]; 
	qtrtype = "threeqtr"; 
	numCPC = 130; 

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
	m.alpha_fric = alpha_fric; 
	m.bridge_method = bridge_method; 
	m.EIS = EIS; 
 	m.r_guess = r_guess; 
	m.v_max_try = v_max_try; 
	m.two_adjust = two_adjust; 



	
    m.N = N; 
    m.T = T;
    m.zeta = zeta; 
    m.winsor_vec = winsor_vec; 
    m.fcitcount = fcitcount; 
    m.numCPC = numCPC; 
    m.nbucket = nbucket; 
    m.nbucket_patent = nbucket_patent; 
    m.winsor_vec_drop = winsor_vec_drop; 
    m.qtrtype = qtrtype; 
    m.prctile_inp = nan; 




end
