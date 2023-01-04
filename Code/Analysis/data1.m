%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: data1.m
% Author: Craig A. Chikis
% Date: 03/08/2021
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = data1()



	m = benchmark_current();
	m_tmp = m; 
 
	m_tmp = workhorse_nested_robust(m_tmp);

		
	rho_vec = power(1 + sort([[.02,.0025,.05, .01], linspace(1e-3, 0.065, 200)]), 1/12)-1;

	rho_idx = nan(1,3);
	rho_idx(1) = find(abs(rho_vec - power(1.05, 1/12) + 1) == min(abs(rho_vec - power(1.05, 1/12) + 1)), 1);
	rho_idx(2) = find(abs(rho_vec - power(1.02, 1/12) + 1) == min(abs(rho_vec - power(1.02, 1/12) + 1)), 1); 
	rho_idx(3) = find(abs(rho_vec - power(1.0025, 1/12) + 1) == min(abs(rho_vec - power(1.0025, 1/12) + 1)), 1); 

		

	dg_store = zeros(1, length(rho_vec));
	xe_store = zeros(1, length(rho_vec));
	x_store = zeros(length(rho_vec), length(-m.s_max:m.s_max));
	markup_store = zeros(1, length(rho_vec));
	r_store = zeros(1,length(rho_vec));
	mu_store = zeros(length(rho_vec), length(0:m.s_max));

	m_tmp_vec = cell(1, length(rho_vec));
	for (ii = 1:length(m_tmp_vec))
		m_iter = m; 
		m_iter.rho = rho_vec(ii);
		m_tmp_vec{ii} = m_iter; 
	end
	parfor (ii = 1:length(rho_vec)) 
		m_tmp_vec{ii} = workhorse_nested_robust(m_tmp_vec{ii});

		xe_store(ii) = m_tmp_vec{ii}.investment_entrant(2);
		dg_store(ii) = m_tmp_vec{ii}.growth_rate; 
		x_store(ii, :) = m_tmp_vec{ii}.investment;
		mu_store(ii, :) = m_tmp_vec{ii}.firm_distributions; 
		markup_store(ii) = m_tmp_vec{ii}.markup; 
		w_store(ii) = m_tmp_vec{ii}.wage_share; 
		r_store(ii) = m_tmp_vec{ii}.real_interest_rate; 
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	idx_choose = round(linspace(1, length(dg_store), 100));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	gspline = spline(rho_vec(idx_choose), dg_store(idx_choose));
	dgspline = fnder(gspline, 1);

	wspline = spline(rho_vec(idx_choose), w_store(idx_choose));
	dwspline = fnder(wspline,1);

	dgdrho = ppval(dgspline, m.rho);
	dwdrho = ppval(dwspline, m.rho);

	markupspline = spline(100*(power(1+rho_vec(idx_choose), 12)-1), markup_store(idx_choose));
	dmarkupspline = fnder(markupspline, 1);

	rspline = spline(rho_vec(idx_choose), r_store(idx_choose));
	drspline = fnder(rspline, 1);

	drdrho = ppval(drspline, m.rho);
	dmarkupdrho = ppval(dmarkupspline, 100*(power(1+m.rho,12)-1));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	xespline = spline(rho_vec(idx_choose), xe_store(idx_choose));
	dxespline = fnder(xespline,1);

	dxedrho = ppval(dxespline, m.rho);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dxdrho = zeros(1, length(-m.s_max:m.s_max));
	for (ii = 1:size(dxdrho, 2))
		xspline_it = spline(rho_vec(idx_choose), x_store(idx_choose, ii));
		dxspline_it = fnder(xspline_it, 1);
		dxdrho(ii) = ppval(dxspline_it, m.rho);
	end
	m_tmp.dxdrho = dxdrho;

	dmu_drho = zeros(1, length(0:m.s_max));
	mu_model_store = cell(size(dmu_drho));
	for (ii = 1:length(dmu_drho))
		mu_spline_it = spline(rho_vec(idx_choose), mu_store(idx_choose, ii));
		mu_model_store{ii} = mu_spline_it;
		dmuspline_it = fnder(mu_spline_it, 1);
		dmu_drho(ii) = ppval(dmuspline_it, m.rho);
	end
	m_tmp.dmu_drho = dmu_drho;

	rho_near2 = find(abs(rho_vec - m.rho) == min(abs(rho_vec - m.rho)),1);

	Dvs = zeros(1,2*m.s_max+1);
	for (ii = 1:m.s_max)
		Dvs(m.s_max+1+ii) = sum(m_tmp.value_functions.*m.prob(m.s_max+1+ii, :)) - m_tmp.value_functions(m.s_max+1+ii);
		Dvs(m.s_max+1-ii) = sum(m_tmp.value_functions.*m.prob(m.s_max+1-ii, :)) - m_tmp.value_functions(m.s_max+1-ii);
	end
	Dvs(m.s_max+1) = sum(m_tmp.value_functions.*m.prob(m.s_max+1,:)) - m_tmp.value_functions(m.s_max+1);
	Dv_omega = Dvs/m_tmp.wage_share;

	m_tmp.Dv_omega = Dv_omega; 

	m_tmp = PE_sys_eq_w_tp(m_tmp);
	m_tmp = PE_entrant(m_tmp);
	m_tmp = incumbent_entrant_strategic(m_tmp);
	m_tmp = strategic_general_tp(m_tmp);
	m_tmp = Thm1_robust(m_tmp, "all");

	intensive = log(m.lambda)*sum([2, ones(1,m.s_max)].*m_tmp.firm_distributions.*m_tmp.dxdrho((m.s_max+1):end)); 
	extensive = log(m.lambda)*sum([2, ones(1,m.s_max)].*m_tmp.dmu_drho.*m_tmp.investment((m.s_max+1):end)); 


	sout2 = "In our estimated model, if we compute these expressions conditional on $\rho=2\%$ and a decline in the discount rate " + ...
				"$d\rho=1\%$, we obtain that the intensive margin effect is " + num2str(ceil(abs(intensive*100)), '%0.0f') + " basis points " + ...
				"and the extensive margin effect is " + num2str(ceil(abs(extensive*100)), '%0.0f') + " basis point. " + ...
				"Following a decline in the discount rate, both the intensive and extensive margins contribute to higher growth. " + ...
				"The positive intensive-margin effect is perhaps unsurprising, because it reflects, in part, the standard valuation or " + ...
				"cost-of-capital channel of a lower discount rate. " + ...
				"The positive extensive margin effect obtains because laggards increase their innovation rates, and laggard innovations " + ...
				"completely close the competitive gap in their industry with probability $\phi>0$ (see Section \ref{sec:results} for a " + ...
				"detailed discussion).  Thus, one way to understand why we obtain that a lower discount rate leads to higher productivity " + ...
				"growth is to observe that both margins (intensive and extensive) are positive.  (Of course, a positive extensive-margin " + ...
				"effect is not necessary for a lower discount rate to contribute to higher growth.)";  

	esigma = eye(2*m.s_max+1);

	term1 = log(m.lambda)*[2, ones(1,m.s_max)]; 
	term2 = esigma((m.s_max+1):end, :) .*  m_tmp.firm_distributions'; 
	term3 = m_tmp.firm_distributions' .* (m_tmp.M((m.s_max+3):(2*m.s_max+3), :) - esigma((m.s_max+1):end, :)); 
	term4 = m_tmp.M((2*m.s_max+4):end, :) .* m_tmp.investment((m.s_max+1):(end-1))'; 

	direct = sum((term2 .* term1') * m_tmp.parx_parrho');
	strategic = sum((term3 .* term1') * m_tmp.parx_parrho');
	composition = sum((term4 .* term1(1:(end-1))') * m_tmp.parx_parrho');
	composition_alt = sum(term1 .* m_tmp.dmu_drho);

	term5 = 0.05 * m_tmp.M(1, :) * sum(esigma((m.s_max+1-10):m.s_max, :)', 2); 

	sout3 = "To interpret the magnitude of the growth multiplier's elements, note that if laggards in competitive industries (with $s \leq 10$) " + ...
				"increase their annual innovation rate by 5 percentage points in response to a lower discount rate when taking as given their " + ...
				"competitors' strategies, then the effect on annual growth would be " + ...
				"$0.05 \times \mathbb{M}_g \sum_{\sigma \in \{-10,..,-1\}} e_\sigma' = " + num2str(term5, '%0.4f') + "$, or " + ...
				num2str(1e4*term5, '%0.0f') + " basis points."; 


	sout1 = "Next, we quantify the valuation, strategic, and composition channels through which a lower discount rate affects growth. " + ...
			 "From \eqref{eq:Mgcomp}, with only the valuation channel (i.e., absent strategic and composition effects), a 100 basis point " + ...
			 "decline in the discount rate increases aggregate productivity by " + ...
			 "$[\ln \lambda \sum_{\sigma \in S^+} (1 + \mathbbm{1}_{\sigma=0}) \mu_\sigma e_{\sigma}] \bm{\partial x}$, or approximately " + ...
			 	num2str(ceil(abs(100*direct)), '%0.0f') + " basis points. " + ...
			 "However, in general equilibrium, the boost to aggregate growth is only about " + num2str(ceil(abs(dgdrho*100)), '%0.0f') + " basis points. " + ...
			 "Overall, the strategic and composition channels \emph{dampen}, but do not overturn, the boost to growth from the valuation channel. " + ...
			 "On balance, the strategic channel decreases growth by " + ...
			 "$-[\ln \lambda \sum_{\sigma \in S^+} (1 + \mathbbm{1}_{\sigma=0})  (\mu_\sigma (\mathbb{M}_{x_\sigma}-e_{\sigma}))] \bm{\partial x}$, " + ...
			 "or " + num2str(ceil(abs(strategic*100)), '%0.0f') + " basis points. " + ...
			 " The composition channel increases growth by " + ...
			 " $[\ln \lambda \sum_{\sigma \in S^+} (1 + \mathbbm{1}_{\sigma=0})  (\mathbb{M}_{\mu_\sigma} x_\sigma ) ] \bm{\partial x},$ or " + ...
			 num2str(ceil(abs(composition_alt*100)), '%0.0f') + " basis point, because of the pro-competitive shift induced by a lower discount rate.";  

	
	string_collect = [sout1, sout2, sout3]; 

	
	m = m_tmp; 
	m = sim_wrap(m); 
	m = significance_construct_sim(m); 
	m = innovation_output(m); 
	m = CSTAT_sim(m);
	m = transmat_sim(m); 
	m = pcm_sim(m);
	m = FHK_sim(m); 



	
	writematrix(string_collect, "Output/Figures_Paper/IA_string.csv"); 
	save(strcat("Output/Store_Data/data1_", "benchmark", ".mat"), 'm', 'm_tmp_vec', 'rho_idx');




end

