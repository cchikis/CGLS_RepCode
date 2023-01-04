%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: lms_compare.m
% Author: Craig A. Chikis
% Date: 01/15/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = lms_compare()
	

	load('Output/Store_Data/lms_export.mat');

	g2 = LMS_bgp.gvec/100;
	rho_vec2 = LMS_bgp.rvec/100;
	rho2 = LMS_bgp.rho/100;
	x2 = LMS_bgp.investment;
	mu2 = LMS_bgp.firm_distributions;
	rho_idx_LMS = find(abs(rho_vec2 - rho2) == min(abs(rho_vec2 - rho2)), 1);


	s_max = 50;
	m = LMS_param_paramcorrect(s_max);

	rho_vec = sort([power(1 + linspace(1e-3, 0.065, 25), 1/12)-1, power(1 + rho2, 1/12)-1]);


	res_vec = cell(size(rho_vec));
	for (ii = 1:length(rho_vec))
		m_tmp = m;
		m_tmp.rho = rho_vec(ii);
		res_vec{ii} = m_tmp;
	end
	parfor (ii = 1:length(rho_vec))
		res_vec{ii} = workhorse_nested_robust(res_vec{ii});
	end



	rho_idx = find(abs(power(1+rho2,1/12)-1 - rho_vec) == min(abs(power(1+rho2,1/12)-1 - rho_vec)), 1);
	xvec = [res_vec{rho_idx}.investment; x2/1200];
	xEvec = [res_vec{rho_idx}.investment_entrant; zeros(1,s_max+1)];
	muvec = [res_vec{rho_idx}.firm_distributions; mu2];
	gvector = [res_vec{rho_idx}.growth_rate; power(1 + g2(rho_idx_LMS), 1/12)-1];
	vvector = [res_vec{rho_idx}.value_functions; res_vec{rho_idx}.value_functions]; 

	m = LMS_param_paramcorrect(s_max);
	m.N = 50000;
    m.T = 12*12; 
	m.winsor_vec = [1, 99]; 
	sim_res = cell(1, size(xvec, 1)); 
	for (ii = 1:size(xvec, 1))
		m_tmp = m; 

		m_tmp.investment = xvec(ii, :);
		m_tmp.investment_entrant = xEvec(ii, :);
		m_tmp.firm_distributions = muvec(ii, :);
		m_tmp.growth_rate = gvector(ii);
		m_tmp.rho = rho_vec(rho_idx);
		m_tmp.wage_share = res_vec{rho_idx}.wage_share;
		m_tmp.value_functions = vvector(ii, :);
		m_tmp.profit = res_vec{rho_idx}.profit;
		m_tmp.nu_s = res_vec{rho_idx}.nu_s;
		m_tmp.lab_demand = res_vec{rho_idx}.lab_demand;


		[m_tmp.p50_approx, m_tmp.markup_array, ~] = markup_quant(m_tmp.kappa, m_tmp.lambda, m_tmp.s_max, m_tmp.firm_distributions, m_tmp.nu_s, 0.5); 
		m_tmp.markup_wt = sum(m_tmp.markup_array(1, :) .* m_tmp.markup_array(2, :)) - 1; 

		m_tmp = sim_wrap(m_tmp);
		m_tmp = CSTAT_sim(m_tmp); 
		m_tmp = innovation_output(m_tmp); 

		sim_res{ii} = m_tmp; 

		
	end


	 rng(2022+08+30, 'Threefry')
	 draw_Lerner = betarnd(1.36, 8, 1e6, 1);
	 draw_markup = 1./(1-draw_Lerner) - 1; 

	 markup_targ = [mean(draw_markup), prctile(draw_markup, [50, 90])]*100; 
	 target_filter_vec = ["three_dropm1"]; 
	 total_prod = 3.44 + -.41 + .76 + 1.23 + .12; 
	 ENTRANCE_data = 1.23/total_prod;
	 WITHIN_data = 3.44/total_prod;
	 BETWEEN_data = -.41/total_prod;
	 CROSS_data = .76/total_prod;
	 EXIT_data = .12/total_prod;
			 
	 FHK_targ = 100*[WITHIN_data*total_prod/(total_prod - ENTRANCE_data*total_prod - EXIT_data*total_prod), ...
							 BETWEEN_data*total_prod/(total_prod - ENTRANCE_data*total_prod - EXIT_data*total_prod), ...
							 CROSS_data*total_prod/(total_prod - ENTRANCE_data*total_prod - EXIT_data*total_prod), ...
							 1e-8,1e-8];
	 FHK_targ_ent = 100*[WITHIN_data, BETWEEN_data, CROSS_data, ENTRANCE_data, EXIT_data];
	 share10_targ = 26.19;
	 share5_targ = 15.68;
			
	 rdsales_targets = readtable(strcat("Output/Store_Data/rdsales_t_targets_", target_filter_vec(1), ".csv")); 
	 profvol_targets = readtable(strcat("Output/Store_Data/profitvol_t_", target_filter_vec(1), ".csv")); 
	 uncond_IO_targets = readmatrix(strcat("Output/Store_Data/kogan_targets.csv"));
	 
	 
			

	g_vec = cellfun(@(x) x.growth_rate, res_vec);
	g2(~(rho_vec2 >= 0.001 & rho_vec2 <= max(rho_vec*12))) = [];
	rho_vec2(~(rho_vec2 >= 0.001 & rho_vec2 <= max(rho_vec*12))) = [];

	close all;
	p1 = plot((rho_vec+g_vec)*1200, g_vec*1200, '--r', 'LineWidth', 2);
	hold on 
	p2 = plot((rho_vec2+g2)*100, g2*100, ':b', 'LineWidth', 2);
	xlabel('Interest rate, $r$', 'Interpreter', 'latex');
	ylabel('Growth rate, $g$', 'Interpreter', 'latex');
	xlim([0, 1200*(rho_vec(end)+g_vec(end))]);
	l1 = legend("CGLS code", "LMS code", 'Interpreter', 'latex');
	set(l1, 'box', 'off');

	saveas(gcf, "Output/Figures_Paper/LMS_rrhofig.png");
	saveas(gcf, "Output/Figures_Paper/LMS_rrhofig.eps", 'epsc');

	d2 = '%0.2f';
	d3 = '%0.3f';
	s1 = strcat("\begin{tabular}{lccccccc}  \hline\hline & \multicolumn{2}{c}{Liu et al.'s (2022) model} \\ \cline{2-3} Moments & CGLS code ", ...
					" & LMS code & Data \\ ", ...					
					" \hline Mean markup  & ", num2str(100*sim_res{1}.markup_wt, d2), "\%", ...
					" & ", num2str(100*sim_res{2}.markup_wt, d2), "\%", ...
					" & ", num2str(markup_targ(1), d2), "\%", " \\  ", ...
	 				" Profit volatility \\ ", ...
					" \hspace{0.2in} All firms & ", num2str(100*sim_res{1}.std_table.mean_fun1_profitgrowth_w(end), d2), "\%", ...
					" & ", num2str(100*sim_res{2}.std_table.mean_fun1_profitgrowth_w(end), d2), "\%", ...
					" & ", num2str(100*profvol_targets.sd(end), d2), "\%", " \\ ", ...
					" \hspace{0.2in} Top profit quintile & ", num2str(100*sim_res{1}.std_table.mean_fun1_profitgrowth_w(1), d2), "\%", ...
					" & ", num2str(100*sim_res{2}.std_table.mean_fun1_profitgrowth_w(1), d2), "\%", ...
					" & ", num2str(100*profvol_targets.sd(1), d2), "\%", " \\ ", ...
					" R\&D to sales \\ ", ...
					" \hspace{0.2in} All firms & ", num2str(100*sim_res{1}.std_table.mean_fun2_rdsales_w(end), d2), "\%", ...
					" & ", num2str(100*sim_res{2}.std_table.mean_fun2_rdsales_w(end), d2), "\%", ...
					" & ", num2str(100*rdsales_targets.p50(end), d2), "\%", " \\ ", ...
					" \hspace{0.2in} Top profit quintile & ", num2str(100*sim_res{1}.std_table.mean_fun2_rdsales_w(1), d2), "\%", ...
					" & ", num2str(100*sim_res{2}.std_table.mean_fun2_rdsales_w(1), d2), "\%", ...
					" & ", num2str(100*rdsales_targets.p50(1), d2), "\%", " \\ ", ...
					" Innovation output \\ ", ...
					" \hspace{0.2in} Mean & ", num2str(100*sim_res{1}.uncond_IO(1), d2), "\%", ...
					" & ", num2str(100*sim_res{2}.uncond_IO(1), d2), "\%", ...
					" & ", num2str(uncond_IO_targets(1), d2), "\%", " \\ ", ...
					" \hspace{0.2in} 90th percentile & ", num2str(100*sim_res{1}.uncond_IO(9), d2), "\%", ...
					" & ", num2str(100*sim_res{2}.uncond_IO(9), d2), "\%", ...
					" & ", num2str(uncond_IO_targets(9), d2), "\%", " \\ ", ...
					"\hline \end{tabular}"); 





	writematrix(s1, "Output/LaTeX_Output/CGLS_LMS.txt", 'QuoteStrings', false);
end