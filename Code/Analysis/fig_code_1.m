%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: fig_code_1.m
% Author: Craig A. Chikis
% Date: 03/08/2021
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = fig_code_1()

	model_run = "benchmark"; 
	load("Output/Store_Data/data1_benchmark.mat")
	prctile_inp = nan;

	
	rng(2022+08+30, 'Threefry')
	draw_Lerner = betarnd(1.36, 8, 1, 1e6);
	draw_Markup = 1./(1-draw_Lerner);

	if (m.kappa >= 9999)
		lerner_p = 1-m.lambda.^(-min(0:m.s_max,m.LMS_bar_s));
		markup_p = m.lambda.^min(0:m.s_max, m.LMS_bar_s);
	else
		lerner_p = 1-1./m.markup_wt_array;
		markup_p = m.markup_wt_array; 
	end

	mid_point_markup = conv(m.lambda.^(min(0:m.s_max, m.LMS_bar_s)), [.5,.5], 'valid');
	mid_point_lerner = conv(1-m.lambda.^(-min(0:m.s_max, m.LMS_bar_s)), [.5,.5], 'valid');
	frac = zeros(1,m.s_max+1);
	frac_l = zeros(1,m.s_max+1);
	for (ii = 2:length(mid_point_markup)-1)
		frac(ii) = sum(draw_Markup >= mid_point_markup(ii-1) & (draw_Markup < mid_point_markup(ii)));
		frac_l(ii) = sum(draw_Lerner >= mid_point_lerner(ii-1) & (draw_Lerner < mid_point_lerner(ii)));
	end 
	frac_l(1) = sum(draw_Lerner < mid_point_lerner(1));
	frac_l(end) = sum(draw_Lerner >= mid_point_lerner(end));

	frac(1) = sum(draw_Markup < mid_point_markup(1));
	frac(end) = sum(draw_Markup >= mid_point_markup(end));
	frac = frac/sum(frac);

	frac_l = frac_l/sum(frac_l);



	kog_targets = readmatrix("Output/Store_Data/kogan_targets.csv");
	cdf_data = readtable("Output/Store_Data/kpss1.csv");
	cdf_data_nber = readtable("Output/Store_Data/nber_cit.csv");
	profvol_targets = readtable(strcat("Output/Store_Data/profitvol_t_", "three_dropm1", ".csv")); 

	data_transmat_main = readmatrix(strcat("Output/Store_Data/transmat_targets_", "three_dropm1", ".csv")); 
	data_transmat_top10 = readmatrix(strcat("Output/Store_Data/transmat_targets_bottom90_", "three_dropm1.csv"));
	data_transmat_top20 = readmatrix(strcat("Output/Store_Data/transmat_targets_top10top20_", "three_dropm1.csv"));
	data_transmat_halves = readmatrix(strcat("Output/Store_Data/transmat_targets_halves_", "three_dropm1", ".csv")); 

	bottom60top20 = sum(data_transmat_main(3:5, 1)); 
	top10bottom90 = sum(data_transmat_top10(1,2));
	top50bottom50 = sum(data_transmat_halves(1,2)); 
	top20top20 = sum(data_transmat_top20(1,1)); 

	bottom60top20_model = sum(m.transmat_out(3:5, 1)); 
	top10bottom90_model = sum(m.transmat_out_top10(1,2));
	top50bottom50_model = sum(m.transmat_out_halves(1,2));
	top20top20_model = sum(m.transmat_out_top20(1,1));

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

	uncond_IO_targets = readmatrix(strcat("Output/Store_Data/kogan_targets.csv"));




	cdf_data.wku = round(cdf_data.wku); 
	[~, idx_unique_wku] = unique(cdf_data.wku);
	[~, idx_unique_wku_nber] = unique(cdf_data_nber.wku);
	pmf_data = cdf_data(idx_unique_wku, :);
	pmf_data_nber = cdf_data_nber(idx_unique_wku_nber, :);
	pmf_data.fcit020 = round(pmf_data.fcit020); 
	pmf_data_nber.fcit020 = round(pmf_data_nber.fcit020);
	unique_ints = unique(pmf_data.fcit020);
	unique_ints_nber = unique(pmf_data_nber.fcit020); 
	data_pmf = zeros(size(unique_ints));
	for (ii = 1:length(unique_ints))
		data_pmf(ii) = sum(pmf_data.fcit020 == unique_ints(ii));
	end
	data_pmf_nber = zeros(size(unique_ints_nber)); 
	for (ii = 1:length(unique_ints_nber))
		data_pmf_nber(ii) = sum(pmf_data_nber.fcit020 == unique_ints_nber(ii)); 
	end
	data_pmf = data_pmf/sum(data_pmf);
	data_pmf_nber = data_pmf_nber/sum(data_pmf_nber);


	close all 
	figure;

	set(gcf, 'PaperUnits', 'inches');
	x_width=2.4;
	y_width=1.75;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

	
	p1 = plot(unique_ints, 100*data_pmf, '-k', 'LineWidth', 2);
	hold on ;
	p2 = plot(unique_ints_nber, 100*data_pmf_nber, '--r', 'LineWidth', 2);
	xlim([0,50]);
	xlabel('Number of citations received', 'Interpreter', 'latex');
	ylabel('Probability (\%)', 'Interpreter', 'latex');
	l1 = legend('KPST', 'NBER', 'Interpreter', 'latex');
	set(l1, 'box', 'off');
	
	saveas(gcf, "Output/Figures_Paper/pmf_compare.png");
	saveas(gcf, "Output/Figures_Paper/pmf_compare.eps", 'epsc');

	m.citation_track2.wku = round(m.citation_track2.wku); 
	m.citation_track2.fcit020 = round(m.citation_track2.fcit020);
	[~, idx_unique_wku_model] = unique(m.citation_track2.wku); 
	citation_track2_uniquewkus = m.citation_track2(idx_unique_wku_model, :); 
	unique_ints_model = unique(m.citation_track2.fcit020);
	model_pmf = zeros(size(unique_ints_model)); 
	for (ii = 1:length(unique_ints_model))
		model_pmf(ii) = sum(m.citation_track2.fcit020 == unique_ints_model(ii)); 
	end
	model_pmf = model_pmf/sum(model_pmf); 


	close all
	figure; 
	
	set(gcf, 'PaperUnits', 'inches');
	x_width=3.25;
	y_width=2.5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


	p1 = plot(unique_ints, 100*data_pmf, '-b', 'LineWidth', 2);
	hold on ;
	p2 = plot(unique_ints_model, 100*model_pmf, '--r', 'LineWidth', 2);
	xlim([0,50]);
	xticks([0, 25, 50]);
	xlabel('Number of citations received', 'Interpreter', 'latex');
	ylabel('Probability (\%)', 'Interpreter', 'latex');
	l1 = legend([p1,p2], 'Data', 'Model', 'Interpreter', 'latex', 'Location', 'northeast');
	set(l1,'box','off');

	saveas(gcf, "Output/Figures_Paper/citdist.png");
	saveas(gcf, "Output/Figures_Paper/citdist.eps", 'epsc');


	
	close all
	figure;

	set(gcf, 'PaperUnits', 'inches');
	x_width=1.1*6.5;
	y_width=1.1*2.5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

	xprc = [1,5,10,25,50,75,90,95,99]; 
	subplot(1,2,2);
	p1 = plot(xprc, uncond_IO_targets(3:end), '-ob', 'LineWidth', 2);
	hold on ;
	p2 = plot(xprc, m.uncond_IO(3:end)*100, '--or', 'LineWidth', 2);
	title({'Innovation output distribution'; ''}, 'Interpreter', 'latex');
	xlabel('Percentile', 'Interpreter', 'latex');
	ylabel('Innovation output (\%)', 'Interpreter', 'latex');
	

	subplot(1,2,1);
	p2 = plot(100*(markup_p-1), 100*m.firm_distributions, '--r', 'LineWidth', 2);
	hold on ;
	p1 =plot(100*(markup_p-1), 100*frac, '-b', 'LineWidth', 2);
	title({'Markup distribution', ''}, 'Interpreter', 'latex');
	xlabel('Net markup (\%)', 'Interpreter', 'latex');
	ylabel('Probability (\%)', 'Interpreter', 'latex');
	xlim([0,100]);
	l3 = legend([p1, p2], 'Data', 'Model', 'Interpreter', 'latex', 'Location', 'northeast');
	set(l3, 'box', 'off');

	
	saveas(gcf, strcat("Output/Figures_Paper/fig1_", ...
					   model_run, ".eps"), 'epsc');
	saveas(gcf, strcat("Output/Figures_Paper/fig1_", ...
					   model_run, ".png"));



	close all 
	figure;

	
	set(gcf, 'PaperUnits', 'inches');
	x_width=8;
	y_width=3.5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]);

	subplot(1,3,1);
	X = categorical({'Q1','Q2','Q3','Q4', 'Q5'});
	X = reordercats(X,{'Q1','Q2','Q3','Q4', 'Q5'});
	b1 = bar(X, 100*[profvol_targets.sd(1:5), m.std_table.mean_fun1_profitgrowth_w(1:5)], 'FaceColor', 'flat');
	b1(1).CData = [0,0,1];
	b1(2).CData = [1,0,0]; 
	title({'Profit volatility'; ''}, 'Interpreter', 'latex');
	ax = gca;
	ax.XAxis.TickLabelInterpreter= 'latex';
	ax.XAxis.FontSize = 8;
	xlabel('Size quintile, Q5 is smallest quintile', 'Interpreter', 'latex');
	l1 = legend('Data', 'Model', 'Interpreter', 'latex', 'Location', 'northwest');
	set(l1, 'box', 'off');


	subplot(1,3,2)
	X = categorical({'Top 50\% to Bottom 50\%', 'Top 10\% to Bottom 90\%', 'Top 20\% to Top 20\%', 'Bottom 60\% to Top 20\%'});
	X = reordercats(X,{'Top 50\% to Bottom 50\%', 'Top 10\% to Bottom 90\%', 'Top 20\% to Top 20\%', 'Bottom 60\% to Top 20\%'});
	b1 = bar(X, 100*[top50bottom50, top50bottom50_model;
			     	 top10bottom90, top10bottom90_model;
				 	 top20top20, top20top20_model;
				 	 bottom60top20, bottom60top20_model], ...
			'FaceColor', 'flat');
	b1(1).CData = [0,0,1];
	b1(2).CData = [1,0,0]; 
	ax = gca;
	ax.XAxis.TickLabelInterpreter= 'latex';
	ax.XAxis.FontSize = 8;

	title({'Transition rates'; ''}, 'Interpreter', 'latex');

	subplot(1,3,3);
	X = categorical({'WITHIN', 'BETWEEN', 'CROSS'});
	X = reordercats(X, {'WITHIN', 'BETWEEN', 'CROSS'}); 
	b1 = bar(X, [FHK_targ(1), 100*m.WITHIN_out;
				 FHK_targ(2), 100*m.BETWEEN_out;
			 	 FHK_targ(3), 100*m.CROSS_out], 'FaceColor', 'flat');
	b1(1).CData = [0,0,1];
	b1(2).CData = [1,0,0]; 
	title({'FHK decomposition'; ''}, 'Interpreter', 'latex');
	ax = gca;
	ax.XAxis.TickLabelInterpreter= 'latex';
	ax.XAxis.FontSize = 8;

	saveas(gcf, "Output/Figures_Paper/2x2_new.eps", 'epsc');
	saveas(gcf, "Output/Figures_Paper/2x2_new.png");

    rng(2022+08+30, 'Threefry');
    draw_Lerner = betarnd(1.36, 8, 1e6, 1);
    draw_markup = 1./(1-draw_Lerner) - 1; 

    target_filter_vec = ["three_dropm1"]; 

    rdsales_targets = readtable(strcat("Output/Store_Data/rdsales_t_targets_", target_filter_vec(1), ".csv")); 
    transmat_targets = readmatrix(strcat("Output/Store_Data/transmat_targets_", target_filter_vec(1), ".csv")); 
    reg_targets = readtable(strcat("Output/Store_Data/reg_targets_", target_filter_vec(1), ".csv"));
    profvol_targets = readtable(strcat("Output/Store_Data/profitvol_t_", target_filter_vec(1), ".csv")); 
    uncond_IO_targets = readmatrix(strcat("Output/Store_Data/kogan_targets.csv"));
    cdf_targs = readmatrix("Output/Store_Data/cdf_current.csv");
	markup_targets = [mean(draw_markup), prctile(draw_markup, [50, 90])]; 
    growth_targ = 1.033; 
   

    clear draw_Lerner;
    clear draw_markup;

	t1 = strcat("\begin{tabular}{ccc|clcc}  \hline \hline  \multicolumn{2}{c}{Parameter estimates} & & &  ", ...
				" \multicolumn{3}{c}{Moments used in estimation} \\  \cline{1-2} \cline{5-7}  ",  ...
				" \multicolumn{1}{l}{Parameter} &  \multicolumn{1}{c}{Value} & & &  ", ...
				" \multicolumn{1}{l}{Description} &  \multicolumn{1}{c}{Model} & \multicolumn{1}{c}{Data} \\ ", ...
				" \hline $\phi$ & ", num2str(m.phi_wt, '%0.3f'), " &  &  & Productivity growth & ", ...
					num2str(100*(power(1+m.growth_rate,12)-1), '%0.2f'), "\%", " & ", num2str(growth_targ, '%0.2f'), "\%", " \\ ", ...
				" $\lambda$ & ", num2str(m.lambda, '%0.3f'), ...
				" & & &  Mean markup & ", num2str(m.markup*100, '%0.2f'), "\%", " & ", ...
				 	num2str(markup_targets(1)*100, '%0.2f'), "\%", " \\ ", ...
				 " $B$ & ", num2str(12*m.B, '%0.3f'), " & & & Profit volatility & & \\ ", ...
				" & & & & \hspace{0.2in} All firms & ", num2str(100*m.std_table.mean_fun1_profitgrowth_w(end), '%0.2f'), "\%", ...
					" & ", num2str(profvol_targets.sd(end)*100, '%0.2f'), "\%", " \\ ", ...
				 " & & & & \hspace{0.2in} Top profit quintile & ", num2str(100*m.std_table.mean_fun1_profitgrowth_w(1), '%0.2f'), "\%", ...
				 " & ", num2str(profvol_targets.sd(1)*100, '%0.2f'), "\%", " \\ ", ...
				 " & & & & R\&D to sales & & \\ & & & & \hspace{0.2in} All firms & ", ...
				 	num2str(100*m.std_table.mean_fun2_rdsales_w(end), '%0.2f'), "\%", " & ", ...
					num2str(100*rdsales_targets.p50(end), '%0.2f'), "\%", " \\ ", ...
				 " & & & & \hspace{0.2in} Top profit quintile & ", num2str(100*m.std_table.mean_fun2_rdsales_w(1), '%0.2f'), "\%", ...
				 	" & ",  num2str(100*rdsales_targets.p50(1), '%0.2f'), "\%", " \\ ", ...
				 " & & & & Innovation output & & \\ & & & & \hspace{0.2in} Mean & ", ...
				 	num2str(100*m.uncond_IO(1), '%0.2f'), "\%", " & ", ...
					num2str(uncond_IO_targets(1), '%0.2f'), "\%", " \\ ", ...
				 " & & & & \hspace{0.2in} 90th percentile & ", num2str(100*m.uncond_IO(9), '%0.2f'), "\%", ...
				 	" & ",  num2str(uncond_IO_targets(9), '%0.2f'), "\%", " \\ \hline \end{tabular} ");



	writematrix(t1, "Output/LaTeX_Output/t1.txt", 'QuoteStrings', false); 



	
	close all
	figure;

	set(gcf, 'PaperUnits', 'inches');
	x_width=6.5;
	y_width=2.5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


	subplot(1,2,1)
	p3 = plot(-m.s_max:m.s_max, 12*m_tmp_vec{rho_idx(1)}.investment, 'LineStyle', ':', 'color', [0,.5,0], 'LineWidth', 2);
	hold on ;
	p1 = plot(-m.s_max:m.s_max, 12*m_tmp_vec{rho_idx(2)}.investment, '-b', 'LineWidth', 1.75);
	p2 = plot(-m.s_max:m.s_max, 12*m_tmp_vec{rho_idx(3)}.investment, '--r', 'LineWidth', 2);
	title({'Cross section of firm innovation rates, $x_\sigma$', ''}, 'Interpreter', 'latex');
	xline(0, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	xlabel('Firm technology position, $\sigma$', 'Interpreter', 'latex');	
	xlim([-10,30]);
	l1 = legend([p3, p1,p2], strcat("$\rho = $ ", num2str(100*(power(1+m_tmp_vec{rho_idx(1)}.rho,12)-1)), "\%"), ...
			   strcat("$\rho = $ ", num2str(100*(power(1+m_tmp_vec{rho_idx(2)}.rho,12)-1)), "\%"), ...
			   strcat("$\rho = $ ", num2str(100*(power(1+m_tmp_vec{rho_idx(3)}.rho,12)-1)), "\%"), ...
			   'Interpreter','latex','Location', 'northeast');
	set(l1, 'box', 'off');

	s3 = "However, laggards, no matter how far behind, maintain an annual innovation rate of " + ...
			"at least " + num2str(round(1200*min([m_tmp_vec{rho_idx(2)}.investment(1:m.s_max)], [], 'all'), -1)) + ...
			"\% per year.";  
	s4 = "Growth therefore increases from about " + num2str(round(1200*m_tmp_vec{rho_idx(2)}.growth_rate, 1), '%0.1f') + "\%" + ... 
			" to " + num2str(round(1200*m_tmp_vec{rho_idx(3)}.growth_rate, 1), '%0.1f') + "\%, reflecting an intensive margin " + ... 
		 	"effect (higher innovation rates conditional on a firm's position) and an extensive " + ...
			"margin effect (a greater share of industries in high-R\&D, competitive industries).";



	subplot(1,2,2)
	plot(0:m.s_max, m_tmp_vec{rho_idx(1)}.firm_distributions, 'LineStyle', ':', 'color', [0,.5,0], 'LineWidth', 2);
	hold on ;
	plot(0:m.s_max, m_tmp_vec{rho_idx(2)}.firm_distributions, '-b', 'LineWidth', 1.75);
	plot(0:m.s_max, m_tmp_vec{rho_idx(3)}.firm_distributions, '--r', 'LineWidth', 2);
	title({'Distribution of technology gaps, $\mu_s$', ''}, 'Interpreter', 'latex');
	xlabel('Industry technology gap, $s$', 'Interpreter', 'latex');
	xlim([0,30]);
	ylim([0,.12]);
	yticks([0:.02:.12]); 

	saveas(gcf, strcat("Output/Figures_Paper/bgp_obj_", ...
					   model_run, ".eps"), 'epsc')
	

	g_vec_whiteboard = zeros(1, length(m_tmp_vec)); 
	r_whiteboard = zeros(1, length(m_tmp_vec)); 
	markup_whiteboard = zeros(size(g_vec_whiteboard));
	mu_store = zeros(size(g_vec_whiteboard, 2), m.s_max+1);
	for (ii = 1:length(m_tmp_vec))
		g_vec_whiteboard(ii) = m_tmp_vec{ii}.growth_rate;
		r_whiteboard(ii) = m_tmp_vec{ii}.real_interest_rate;
		markup_whiteboard(ii) = m_tmp_vec{ii}.markup; 
		mu_store(ii, :) = m_tmp_vec{ii}.firm_distributions; 
	end

	g_vec = linspace(power(1.007,1/12)-1, power(1.013,1/12)-1, 15);
	r_rho2 = g_vec + m.rho;
	r_rho1 = g_vec + (power(1.045,1/12)-1);
	r_rho0 = g_vec + (power(1+2.5e-3,1/12)-1);

	growth_rate_whiteboard = g_vec_whiteboard;
	real_int_whiteboard = r_whiteboard;
	r_rho0_whiteboard = r_rho0;
	r_rho2_whiteboard = r_rho2; 
	r_rho1_whiteboard = r_rho1; 
	g_vec_whiteboard = g_vec;




	close all
	figure;

	set(gcf, 'PaperUnits', 'inches');
	x_width=6.5;
	y_width=2.5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


	subplot(1,2,1)
	p1 = plot(100*(power(1+real_int_whiteboard,12)-1), 100*(power(1+growth_rate_whiteboard,12)-1), '-b', 'LineWidth', 2);
	hold on ;
	p2 = plot(100*(power(1+r_rho0_whiteboard,12)-1), 100*(power(1+g_vec_whiteboard,12)-1), 'color', [.5,.5,.5], 'LineStyle', '--', 'LineWidth', 1);
	p3 = plot(100*(power(1+r_rho2_whiteboard,12)-1), 100*(power(1+g_vec_whiteboard,12)-1), 'color', [.33,.33,.33], 'LineStyle', '--', 'LineWidth', 1);
	ylabel('Growth', 'Interpreter', 'latex');
	xlabel('Interest rate', 'Interpreter', 'latex');
	xlim([0.6,6]);
	ylim([0.75,1.25]);
	title({'Growth and the interest rate', 'in general equilibrium'}, 'Interpreter', 'latex');
	l1 = legend([p1, p2, p3], sprintf('Innovation\nschedule'), ...
						  '$\rho = 0.25\%$', ...
						  '$\rho = 2.00\%$', ...
						   'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 7);
	set(l1, 'box', 'off');

	subplot(1,2,2)
	plot(100*(power(1+real_int_whiteboard,12)-1), 100*markup_whiteboard, '-b', 'LineWidth', 2);
	title({'Mean net markup', ''}, 'Interpreter', 'latex');
	xlabel('Interest rate', 'Interpreter', 'latex');
	ylabel('Percent', 'Interpreter', 'latex');
	xlim([0.6,6]);




	saveas(gcf, strcat("Output/Figures_Paper/Figure5_remake_", ...
					   model_run, ".eps"), 'epsc')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	close all
	figure;


	set(gcf, 'PaperUnits', 'inches');
	x_width=7.5;
	y_width=2.5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

	subplot(1,3,1)
	p1 = plot(-m.s_max:m.s_max, -m.parx_parrho/100, '-r', 'LineWidth', 1.75);
	title({'$\mathbf{\partial x}$'; ''}, 'Interpreter', 'latex');
	xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
	xline(0, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	xlim([-30,30]);

	subplot(1,3,2)
	p2 = plot(-m.s_max:m.s_max, m.M(1,:), '-k', 'LineWidth', 1.75);
	title({'Growth multiplier $M_{g}$'; ''}, 'Interpreter', 'latex');
	xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
	xline(0, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	p5 = yline(0, '-k', 'LineWidth', 1.25);
	xlim([-30,30]);
	ylim([-.003,.008]);
	ax = gca;
	ax.YRuler.Exponent = 0;
	
	e_kappa = zeros(m.s_max+1, 2*m.s_max+1);
	for (ii = 1:(size(e_kappa,1)))
		e_kappa(ii,m.s_max+1+ii-1) = 1;
	end

	indicator = zeros(m.s_max+1, 2*m.s_max+1);
	for (ii = 1:(size(indicator,1)))
		indicator(ii,m.s_max+1+ii-1) = 1 + (ii-1 == 0);
	end

	direct_effect = zeros(m.s_max+1, 2*m.s_max+1);
	strategic_channel = zeros(size(direct_effect));
	composition_channel = zeros(size(direct_effect));
	for (ii = 1:size(direct_effect,1))
		direct_effect(ii,:) = log(m.lambda)*indicator(ii,m.s_max+ii).*(m.firm_distributions(ii)*e_kappa(ii,:));
		strategic_channel(ii,:) = log(m.lambda)*indicator(ii,m.s_max+ii).*(m.firm_distributions(ii)*(m.M(m.s_max+1+ii,:)-e_kappa(ii,:)));
		composition_channel(ii,:) = log(m.lambda)*indicator(ii,m.s_max+ii).*(m.investment(m.s_max+ii)*m.M(2*m.s_max+2+ii,:));
	end
	direct_effect = sum(direct_effect,1);
	strategic_channel = sum(strategic_channel,1);
	composition_channel = sum(composition_channel,1);


	subplot(1,3,3)
	p2 = plot(-m.s_max:m.s_max, direct_effect, 'LineStyle', '-', 'LineWidth', 2, 'color', [0,.5,0]);
	hold on 
	p3 = plot(-m.s_max:m.s_max, strategic_channel, ':b', 'LineWidth', 2);
	p4 = plot(-m.s_max:m.s_max, composition_channel, '--r', 'LineWidth', 2);
	p6 = yline(0, '-k', 'LineWidth', 1.25);
	xline(0, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	title({'Growth multiplier'; 'components'}, 'Interpreter', 'latex');
	xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
	xlim([-30,30]);
	ylim([-.003,.008]);
	ax = gca;
	ax.YRuler.Exponent = 0;
	l1 = legend([p2,p3,p4], 'Direct', 'Strategic', 'Composition', ...
		   'Interpreter', 'latex', ...
		   'Location', 'northwest', 'FontSize', 5);
	set(l1, 'box', 'off');


	saveas(gcf, "Output/Figures_Paper/mult_decomp_v3.png")
	saveas(gcf, "Output/Figures_Paper/mult_decomp_v3.eps", 'epsc')


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	close all
	figure;
	if (m.elastic_labor)
		add_it = 2;
	else
		add_it = 3;
	end


	set(gcf, 'PaperUnits', 'inches');
	x_width=5;
	y_width=4;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

	subplot(2,1,1)
	p1 = plot(-m.s_max:m.s_max, m.M(add_it+m.s_max+4, :), '-k', 'LineWidth', 1.75);
	hold on
	xline(-4, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	xline(4, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
	ylabel('Multiplier', 'Interpreter', 'latex');
	title({'$M_{x_4}$'; ''}, 'Interpreter', 'latex');
	xlim([-10,10]);


	subplot(2,1,2)
	p1 = plot(-m.s_max:m.s_max, m.parx_parxsig(:, m.s_max+1+4), 'color', [0, 0, 1], 'LineWidth', 2, ...
			  'LineStyle', ':');
	hold on 
	p2 = plot(-m.s_max:m.s_max, [zeros(1,m.s_max+1), 0, 0, 0, 1, zeros(1, m.s_max-4)], 'LineStyle', '-', ...
			  'color', [0,.5,0], 'LineWidth', 1.75);
	xline(-4, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	xline(4, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	xlabel('Technology position, $\iota$', 'Interpreter', 'latex');
	title({''; '$\frac{\partial x_4}{\partial x^c_\iota}$'; ''}, 'Interpreter', 'latex');
	l1 = legend([p2,p1], 'Direct effect', 'Strategic channel', 'Interpreter', 'latex', 'Location', 'northeast');
	set(l1, 'box', 'off');
	xlim([-10,10]);

	saveas(gcf, strcat("Output/Figures_Paper/strat_M4_", ...
					   model_run, ".png"))
	saveas(gcf, strcat("Output/Figures_Paper/strat_M4_", ...
					   model_run, ".eps"), 'epsc')


	close all
	figure; 

	
	set(gcf, 'PaperUnits', 'inches');
	x_width=4.3;
	y_width=3;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

	profitshare = zeros(1, size(mu_store, 1)); 
	for (ii = 1:size(mu_store, 1))
		profitshare(ii) = sum((flip(m.profit(1:(m.s_max+1))) + m.profit((m.s_max+1):end)).*mu_store(ii, :)/(1-m.tax_rate(2)));
	end
	p1 = plot(100*(power(1 + real_int_whiteboard, 12) - 1), profitshare*100, '-b', 'LineWidth', 2);
	xlim([1, 7]);
	ylabel('Profit share (\%)', 'Interpreter', 'latex');
	xlabel('Interest rate, $r$', 'Interpreter', 'latex');

	saveas(gcf, "Output/Figures_Paper/profitshare.eps", 'epsc');
	saveas(gcf, "Output/Figures_Paper/profitshare.png");



	s1 = "In the data, about " + num2str(100*data_transmat_halves(1,2), '%0.1f') + "% " + ...
				"of firms in the top half each year are ranked in the bottom" + ...
				"half the following year, close to the value in the model, " + ...
				num2str(100*m.transmat_out_halves(1,2), '%0.1f') + "%.";

	s2 =  "In the data, only " + num2str(data_transmat_top10(1,2)*100, '%0.1f') +  ...
		"% of firms in the top 10% each year are ranked in the bottom 90% the following"+ ...
		" year, also close to the value in the model, " + num2str(m.transmat_out_top10(1, 2)*100, '%0.1f') + "%.";



	s5 = "Figure \ref{fig:nontargeted}, middle panel, shows transition rates." + ...
		  " Following \citet{Acemoglu_et_al_AER_2018}, we start by focusing on the transition rate from the larger " + ...
		  "half of firms to the smaller half.  In the data, about " + num2str(100*top50bottom50, '%0.1f') +  ...
		  "\% of firms in the top half each year are ranked in the bottom half the following year, close to " + ...
		  "the value in the model " + num2str(top50bottom50_model*100, '%0.1f') + "\%.  We also take " + ...
		  "a ``best versus the rest'' perspective and examine exit from the top 10\% of firms.  " + ...
		  "In the data, only " + num2str(top10bottom90*100, '%0.1f') + "\% of firms in the top 10\% "+ ...
		  "each year are ranked in the bottom 90\% the following year, also close to the value in the model.  " + ...
		  "In the model and the data, almost all firms from the top 20\% remain in the top quintile the " + ...
		  "following year, while almost no firms from the bottom 60\% reach the top quintile.  " + ...
		  "Thus, our model matches well salient facts about the persistence and attainment of leadership.";  


	writematrix([s1;s2;s3;s4;s5], "Output/Figures_Paper/number_sentences.csv");





end
