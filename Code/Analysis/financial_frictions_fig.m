%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: financial_frictions_fig.m
% Author: Craig A. Chikis
% Date: 03/09/2021
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = financial_frictions_fig()



	rho_vec = sort([power(1+linspace(1e-3, .05, 30),1/12)-1, power(1.02,1/12)-1]);
	alpha_vec = [.0001/12, .001/12, .01/12, .1/12, 1/12,  Inf];

	rho_star = find(abs(rho_vec - (power(1.02,1/12)-1)) == min(abs(rho_vec - (power(1.02,1/12)-1))),1);

	 
	store_results = cell(3, length(alpha_vec));

	for (ii = 1:size(store_results,1))
		if (ii == 1)
			m = benchmark_current();
		elseif (ii == 2)
			m = LMS_param_paramcorrect();
		elseif (ii == 3)
			m = benchmark_current_fric();
		end

		for (jj = 1:length(alpha_vec))
			m.alpha_fric = alpha_vec(jj); 
			m_iter = cell(1,length(rho_vec));
			for (kk = 1:length(rho_vec))
				m_tmp = m;
				m_tmp.rho = rho_vec(kk); 
				m_iter{kk} = m_tmp;
			end
			parfor (kk = 1:length(rho_vec))

				m_iter{kk} = workhorse_nested_robust(m_iter{kk});
			

			end

			store_results{ii, jj} = m_iter; 

		end


	end


	close all
	figure;

	s_max=  m.s_max;

	set(gcf, 'PaperUnits', 'inches');
	x_width=7.25;
	y_width=2;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

	subplot(1,4,1)
	p1 = plot(-s_max:s_max, 12*store_results{3,5}{rho_star}.investment, '-r', 'LineWidth', 1.75);
	hold on 
	p2 = plot(-s_max:s_max, 12*store_results{3,4}{rho_star}.investment, '--b', 'LineWidth', 1.75);
	p3 = plot(-s_max:s_max, 12*store_results{3,3}{rho_star}.investment, 'LineStyle', ':', ...
	     	  'color', [0,.5,0], 'LineWidth', 1.75);
	xline(0, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	ylabel('$x_\sigma$', 'Interpreter', 'latex');
	xlim([-10,30]);
	title({'Innovation rate'; ''},'Interpreter', 'latex');

	subplot(1,4,2)
	p1 = plot(0:s_max, store_results{3,5}{rho_star}.firm_distributions, '-r', 'LineWidth', 1.75);
	hold on 
	p2 = plot(0:s_max, store_results{3,4}{rho_star}.firm_distributions, '--b', 'LineWidth', 1.75);
	p3= plot(0:s_max, store_results{3,3}{rho_star}.firm_distributions, 'LineStyle', ':', ...
	     'color', [0,.5,0], 'LineWidth', 1.75);
	ylabel('$\mu_s$', 'Interpreter', 'latex');
	xlim([0,30]);
	title({'Gap distribution'; ''}, 'Interpreter', 'latex');
	l1 = legend([p1,p2,p3], strcat("$\beta = $ ", num2str(alpha_vec(5)*12)), ...
		   strcat("$\beta = $ ", num2str(alpha_vec(4)*12)), ...
		   strcat("$\beta = $ ", num2str(alpha_vec(3)*12)), 'Interpreter', 'latex', ...
		   'Location', 'northeast', 'FontSize', 5);
	set(l1, 'box', 'off');

	gplot = zeros(3, length(rho_vec));
	markupplot = zeros(size(gplot));
	for (ii = 1:length(rho_vec))
		gplot(1, ii) = 100*(power(1+store_results{3,5}{ii}.growth_rate,12)-1);
		gplot(2, ii) = 100*(power(1+store_results{3,4}{ii}.growth_rate,12)-1);
		gplot(3, ii) = 100*(power(1+store_results{3,3}{ii}.growth_rate,12)-1);
		
		markupplot(1, ii) = 100*store_results{3,5}{ii}.markup;
		markupplot(2, ii) = 100*store_results{3,4}{ii}.markup;
		markupplot(3, ii) = 100*store_results{3,3}{ii}.markup;

	end
	subplot(1,4,3)
	p1 = plot(100*(power(1+rho_vec,12)-1), gplot(1,:), '-r', 'LineWidth', 1.75);
	hold on 
	p2 = plot(100*(power(1+rho_vec,12)-1), gplot(2,:), '--b', 'LineWidth', 1.75);
	p3 = plot(100*(power(1+rho_vec,12)-1), gplot(3,:), 'LineStyle', ':', ...
	     'color', [0,.5,0], 'LineWidth', 1.75);
	title({'Growth'; ''}, 'Interpreter', 'latex');
	ylabel('Productivity growth, $g$', 'Interpreter', 'latex');
	xlim([0,4]);


	subplot(1,4,4)
	plot(100*(power(1+rho_vec,12)-1), markupplot(1,:), '-r', 'LineWidth', 1.75);
	hold on 
	plot(100*(power(1+rho_vec,12)-1), markupplot(2,:), '--b', 'LineWidth', 1.75);
	plot(100*(power(1+rho_vec,12)-1), markupplot(3,:), 'LineStyle', ':', ...
	     'color', [0,.5,0], 'LineWidth', 1.75);
	title({'Mean net markup'; ''}, 'Interpreter', 'latex');
	ylabel('Mean net markup', 'Interpreter', 'latex');
	xlim([0,4]);

	sgtitle('Estimated model under alternative pledgeability assumptions', 'Interpreter', 'latex');

	saveas(gcf, "Output/Figures_Paper/finfric_update1.png")
	saveas(gcf, "Output/Figures_Paper/finfric_update1.eps", 'epsc')

	close all
	figure;

	set(gcf, 'PaperUnits', 'inches');
	x_width=7.25;
	y_width=2;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

	subplot(1,4,1)
	plot(-s_max:s_max, 12*store_results{2,3}{rho_star}.investment, '-r', 'LineWidth', 1.75);
	hold on 
	plot(-s_max:s_max, 12*store_results{2,2}{rho_star}.investment, '--b', 'LineWidth', 1.75);
	plot(-s_max:s_max, 12*store_results{2,1}{rho_star}.investment, 'LineStyle', ':', ...
	     'color', [0,.5,0], 'LineWidth', 1.75);
	xline(0, 'LineStyle', ':', 'color', [.7,.7,.7], 'LineWidth', 1);
	ylabel('$x_\sigma$', 'Interpreter', 'latex');
	xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
	xlim([-10,30]);
	title({'Innovation rate'; ''}, 'Interpreter', 'latex');

	subplot(1,4,2)
	plot(0:s_max, store_results{2,3}{rho_star}.firm_distributions, '-r', 'LineWidth', 1.75);
	hold on 
	plot(0:s_max, store_results{2,2}{rho_star}.firm_distributions, '--b', 'LineWidth', 1.75);
	plot(0:s_max, store_results{2,1}{rho_star}.firm_distributions, 'LineStyle', ':', ...
	     'color', [0,.5,0], 'LineWidth', 1.75);
	ylabel('$\mu_s$', 'Interpreter', 'latex');
	xlabel('Technology gap, $s$', 'Interpreter', 'latex');
	xlim([0,30]);
	title({'Gap distribution'; ''}, 'Interpreter', 'latex');
	l1 = legend(strcat("$\beta = $ ", num2str(alpha_vec(3)*12)), ...
		   strcat("$\beta = $ ", num2str(alpha_vec(2)*12)), ...
		   strcat("$\beta = $ ", num2str(alpha_vec(1)*12)), 'Interpreter', 'latex', ...
		   'Location', 'northeast', 'FontSize', 5);
	set(l1, 'box', 'off');

	gplot2 = zeros(3, length(rho_vec));
	markupplot2 = zeros(size(gplot));
	for (ii = 1:length(rho_vec))
		gplot2(1, ii) = 100*(power(1+store_results{2,3}{ii}.growth_rate,12)-1);
		gplot2(2, ii) = 100*(power(1+store_results{2,2}{ii}.growth_rate,12)-1);
		gplot2(3, ii) = 100*(power(1+store_results{2,1}{ii}.growth_rate,12)-1);
		
		markupplot2(1, ii) = 100*(store_results{2,3}{ii}.markup_wt-1);
		markupplot2(2, ii) = 100*(store_results{2,2}{ii}.markup_wt-1);
		markupplot2(3, ii) = 100*(store_results{2,1}{ii}.markup_wt-1);

	end

	subplot(1,4,3)
	plot(100*(power(1+rho_vec,12)-1), gplot2(1,:), '-r', 'LineWidth', 1.75);
	hold on 
	plot(100*(power(1+rho_vec,12)-1), gplot2(2,:), '--b', 'LineWidth', 1.75);
	plot(100*(power(1+rho_vec,12)-1), gplot2(3,:), 'LineStyle', ':', ...
	     'color', [0,.5,0], 'LineWidth', 1.75);
	title({'Growth'; ''}, 'Interpreter', 'latex');
	ylabel('Productivity growth, $g$', 'Interpreter', 'latex');
	xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
	xlim([0,4]);

	subplot(1,4,4)
	plot(100*(power(1+rho_vec,12)-1), markupplot2(1,:), '-r', 'LineWidth', 1.75);
	hold on 
	plot(100*(power(1+rho_vec,12)-1), markupplot2(2,:), '--b', 'LineWidth', 1.75);
	plot(100*(power(1+rho_vec,12)-1), markupplot2(3,:), 'LineStyle', ':', ...
	     'color', [0,.5,0], 'LineWidth', 1.75);
	title({'Mean net markup'; ''}, 'Interpreter', 'latex');
	ylabel('Mean net markup', 'Interpreter', 'latex');
	xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
	xlim([0,4]);

	sgtitle('LMS (2022) under alternative pledgeability assumptions', 'Interpreter', 'latex');


	saveas(gcf, "Output/Figures_Paper/finfric_update2.png");
	saveas(gcf, "Output/Figures_Paper/finfric_update2.eps", 'epsc');







end
