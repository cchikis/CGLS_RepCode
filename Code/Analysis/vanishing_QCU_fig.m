%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: vanishing_QCU_fig.m
% Author: Craig A. Chikis
% Date: 12/17/2021
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = vanishing_QCU_fig()

    target_filter_vec = "three_dropm1"; 
    kogan_str = "Output/Store_Data/kogan_targets.csv"; 

    rdsales_targets = readtable(strcat("Output/Store_Data/rdsales_t_targets_", target_filter_vec(1), ".csv")); 
    profvol_targets = readtable(strcat("Output/Store_Data/profitvol_t_", target_filter_vec(1), ".csv")); 
    uncond_IO_targets = readmatrix(kogan_str(1));


    m = eta3p_param();

    write_log = false;

    rho_vec = sort([power(1 + linspace(1e-3, 0.07, 25), 1/12) - 1, power(1.02, 1/12) - 1]);
    bar_j_vec = [round(log(2)/log(m.lambda)), round(log(2.5)/log(m.lambda)), 68, 78, m.s_max];

    alpha_vec = [log(0.1/m.phi_wt)/(1-log(2)/log(m.lambda)), log(0.1/m.phi_wt_exog)/(1-log(2)/log(m.lambda))];
    nu_vec = [max(alpha_vec, 0);
            0, 0;
            0.01, 0;
            0.005, 0;
            0.001, 0];
    exer = {alpha_vec, nu_vec}; 

    rng(2022+08+30, 'Threefry')
    draw_Lerner = betarnd(1.36, 8, 1e6, 1);
    draw_markup = 1./(1-draw_Lerner) - 1; 

    markup_targ = 100*[mean(draw_markup), prctile(draw_markup, [50, 90])]; 
    clear draw_Lerner;
    clear draw_markup; 
    smm_err = zeros(length(exer),size(nu_vec, 1),length(rho_vec));
    weights = [5, 5, ones(1, 6)]; 
    targets = [1.033, markup_targ(1), 100*profvol_targets.sd(end), 100*profvol_targets.sd(1), 100*rdsales_targets.p50(end), ...
                100*rdsales_targets.p50(1), uncond_IO_targets(1), uncond_IO_targets(9)]; 
    m_collect = cell(length(exer), size(nu_vec, 1), length(rho_vec)); 
    for (ii = 1:length(exer))
        for (jj = 1:size(nu_vec, 1))

            m = eta3p_param(); 


            if (ii == 1)
                jump_max = bar_j_vec(jj);
                jump_max_p = bar_j_vec(jj);
                jump_max_E = bar_j_vec(jj);

                nu_phi = 0;
                nu_zeta = 0;
                nu_phi_E = 0;
            elseif (ii == 2)
                nu_phi = nu_vec(jj, 1);
                nu_zeta = nu_vec(jj, 2);
                nu_phi_E = nu_vec(jj, 1);

                jump_max = m.s_max;
                jump_max_p = m.s_max;
                jump_max_E = m.s_max;
            end
            if (ii == 1)
                phi_bar = m.phi_wt;
                zeta_bar = m.phi_wt_exog;
                phi_E_bar = m.phi_wt_tilde;
            elseif (ii == 2)
                phi_bar = m.phi_wt;
                zeta_bar =0;
                phi_E_bar = m.phi_wt_tilde;
            end                    

            [prob, probE, prob_exog] = kernel_exp_s(phi_bar, zeta_bar, phi_E_bar,...
                                                    nu_phi, nu_zeta, nu_phi_E, ...
                                                    m.s_max, jump_max, jump_max_p, jump_max_E);

            
            m.prob = prob;
            m.probE = probE; 
            m.prob_exog = prob_exog; 


            m = workhorse_nested_robust(m); 
            m = sim_wrap(m);
            m = CSTAT_sim(m); 
            m = innovation_output(m); 

            m_vec_iter = cell(1, length(rho_vec)); 
            for (kk = 1:length(rho_vec))
                m_tmp = m; 
                m_tmp.rho = rho_vec(kk); 
                m_vec_iter{kk} = m_tmp; 
            end

            smm_iter = nan(1, length(rho_vec)); 
            parfor (kk = 1:length(rho_vec))
                try
                    m_vec_iter{kk} = workhorse_nested_robust(m_vec_iter{kk}); 
                    m_vec_iter{kk} = sim_wrap(m_vec_iter{kk});
                    m_vec_iter{kk} = CSTAT_sim(m_vec_iter{kk});
                    m_vec_iter{kk} = innovation_output(m_vec_iter{kk});

                    output = [100*(power(1+m_vec_iter{kk}.growth_rate,12)-1), ...
                            100*m_vec_iter{kk}.markup, ...
                            100*m_vec_iter{kk}.std_table.mean_fun1_profitgrowth_w(end), ...
                            100*m_vec_iter{kk}.std_table.mean_fun1_profitgrowth_w(1), ...
                            100*m_vec_iter{kk}.std_table.mean_fun2_rdsales_w(end), ...
                            100*m_vec_iter{kk}.std_table.mean_fun2_rdsales_w(1), ...
                            100*m_vec_iter{kk}.uncond_IO(1), ...
                            100*m_vec_iter{kk}.uncond_IO(9)]; 

                    smm_iter(kk) = sum(weights.*abs(output - targets)./(0.5*(abs(output) + abs(targets))));
                catch
                end

            end
            smm_err(ii, jj, :) = smm_iter; 
            m_collect(ii, jj, :) = m_vec_iter; 
          

        end

    end


    rho_idx = find(abs(rho_vec - (power(1.02,1/12)-1)) == min(abs(rho_vec - (power(1.02,1/12)-1))), 1);


    vary_smm = [round(linspace(5, 35, 10)), 36:40, round(linspace(41, m.s_max, 10));
                linspace(0, 0.01, 15), linspace(0.02, 0.1, 10)];
    smm_val = nan(size(vary_smm));
    rho = rho_vec(rho_idx);
    for (ii = 1:size(vary_smm, 1))
        for (jj = 1:size(vary_smm, 2))
            try
                m = eta3p_param(); 
                
                if (ii == 1)
                    jump_max = vary_smm(ii, jj);
                    jump_max_p = vary_smm(ii, jj);
                    jump_max_E = vary_smm(ii, jj);
                
                    nu_phi = 0;
                    nu_zeta = 0;
                    nu_phi_E = 0;
                elseif (ii == 2)
                    nu_phi = vary_smm(ii, jj);
                    nu_zeta = vary_smm(ii, jj);
                    nu_phi_E = vary_smm(ii, jj);
                
                    jump_max = m.s_max;
                    jump_max_p = m.s_max;
                    jump_max_E = m.s_max;
                end
                if (ii == 1)
                    phi_bar = m.phi_wt;
                    zeta_bar = m.phi_wt_exog;
                    phi_E_bar = m.phi_wt_tilde;
                elseif (ii == 2)
                    phi_bar = m.phi_wt;
                    zeta_bar = 0;
                    phi_E_bar = m.phi_wt_tilde;
                end

                [prob, probE, prob_exog] = kernel_exp_s(phi_bar, zeta_bar, phi_E_bar,...
                                                        nu_phi, nu_zeta, nu_phi_E, ...
                                                        m.s_max, jump_max, jump_max_p, jump_max_E);

                
                m.prob = prob;
                m.probE = probE; 
                m.prob_exog = prob_exog; 

                m = workhorse_nested_robust(m);
                m = sim_wrap(m);
                m = CSTAT_sim(m); 
                m = innovation_output(m); 

                output = [100*(power(1+m.growth_rate,12)-1), ...
                            100*m.markup, ...
                            100*m.std_table.mean_fun1_profitgrowth_w(end), ...
                            100*m.std_table.mean_fun1_profitgrowth_w(1), ...
                            100*m.std_table.mean_fun2_rdsales_w(end), ...
                            100*m.std_table.mean_fun2_rdsales_w(1), ...
                            100*m.uncond_IO(1), ...
                            100*m.uncond_IO(9)]; 

                smm_val(ii,jj) = sum(weights.*abs(output - targets)./(0.5*(abs(output) + abs(targets))));
            catch
            end



        end
    end




    d2 = '%0.2f';
    d3 = '%0.3f';

    
    s1 = strcat("\begin{tabular}{lccccccc} \hline \hline  & ", ...
                " \multicolumn{2}{c}{\shortstack{Capped maximum \\ productivity increase}} & & ", ...
                " \multicolumn{2}{c}{\shortstack{Vanishing probability of \\ quick catch-up}} & & Data  ", ...
                " \\  \cline{2-3} \cline{5-6} & & & &  ", " &  & ", "\\ ",...
                " Moments & $\bar{\lambda}=2.0$ & $\bar{\lambda} = 2.5$ & & $\alpha = $ ", num2str(nu_vec(3, 1), d3),...
                "& $\alpha = $ ", num2str(nu_vec(2,1)), ...
                " & &  \\ ", ...
                " \hline  Productivity growth  & ", num2str(100*(power(1+m_collect{1,1,rho_idx}.growth_rate,12)-1), d2), "\%", ...
                " & ", num2str(100*(power(1+m_collect{1,2,rho_idx}.growth_rate,12)-1), d2), "\%", ...
                " & & ", num2str(100*(power(1+m_collect{2,3,rho_idx}.growth_rate,12)-1), d2), "\%", ...
                " & ", num2str(100*(power(1+m_collect{2,2,rho_idx}.growth_rate,12)-1), d2), "\%", ...
                " & & ", num2str(targets(1), d2), "\%", " \\ ", ...
                "  Mean markup  & ", num2str(100*m_collect{1,1,rho_idx}.markup, d2), "\%", ...
                " & ", num2str(100*m_collect{1,2,rho_idx}.markup, d2), "\%", ...
                " & & ", num2str(100*m_collect{2,3,rho_idx}.markup, d2), "\%", ...
                " & ", num2str(100*m_collect{2,2,rho_idx}.markup, d2), "\%", ...
                " & & ", num2str(targets(2), d2), "\%", " \\ ", ...
                " Profit volatility \\ ", ...
                " \hspace{0.2in} All firms  & ", num2str(100*m_collect{1,1,rho_idx}.std_table.mean_fun1_profitgrowth_w(end), d2), "\%", ...
                " & ", num2str(100*m_collect{1,2,rho_idx}.std_table.mean_fun1_profitgrowth_w(end), d2), "\%", ...
                " & & ", num2str(100*m_collect{2,3,rho_idx}.std_table.mean_fun1_profitgrowth_w(end), d2), "\%", ...
                " & ", num2str(100*m_collect{2,2,rho_idx}.std_table.mean_fun1_profitgrowth_w(end), d2), "\%", ...
                " & & ", num2str(targets(3), d2), "\%", " \\ ", ...
                " \hspace{0.2in} Top profit quintile  & ", num2str(100*m_collect{1,1,rho_idx}.std_table.mean_fun1_profitgrowth_w(1), d2), "\%", ...
                " & ", num2str(100*m_collect{1,2,rho_idx}.std_table.mean_fun1_profitgrowth_w(1), d2), "\%", ...
                " & & ", num2str(100*m_collect{2,3,rho_idx}.std_table.mean_fun1_profitgrowth_w(1), d2), "\%", ...
                " & ", num2str(100*m_collect{2,2,rho_idx}.std_table.mean_fun1_profitgrowth_w(1), d2), "\%", ...
                " & & ", num2str(targets(4), d2), "\%", " \\ ", ...
                " R\&D to sales  \\ ", ...
                " \hspace{0.2in} All firms  & ", num2str(100*m_collect{1,1,rho_idx}.std_table.mean_fun2_rdsales_w(end), d2), "\%", ...
                " & ", num2str(100*m_collect{1,2,rho_idx}.std_table.mean_fun2_rdsales_w(end), d2), "\%", ...
                " & & ", num2str(100*m_collect{2,3,rho_idx}.std_table.mean_fun2_rdsales_w(end), d2), "\%", ...
                " & ", num2str(100*m_collect{2,2,rho_idx}.std_table.mean_fun2_rdsales_w(end), d2), "\%", ...
                " & & ", num2str(targets(5), d2), "\%", " \\ ", ...
                " \hspace{0.2in} Top profit quintile  & ", num2str(100*m_collect{1,1,rho_idx}.std_table.mean_fun2_rdsales_w(1), d2), "\%", ...
                " & ", num2str(100*m_collect{1,2,rho_idx}.std_table.mean_fun2_rdsales_w(1), d2), "\%", ...
                " & & ", num2str(100*m_collect{2,3,rho_idx}.std_table.mean_fun2_rdsales_w(1), d2), "\%", ...
                " & ", num2str(100*m_collect{2,2,rho_idx}.std_table.mean_fun2_rdsales_w(1), d2), "\%", ...
                " & & ", num2str(targets(6), d2), "\%", " \\ ", ...
                " Innovation output  \\ ", ...
                " \hspace{0.2in} Mean  & ", num2str(100*m_collect{1,1,rho_idx}.uncond_IO(1), d2), "\%", ...
                " & ", num2str(100*m_collect{1,2,rho_idx}.uncond_IO(1), d2), "\%", ...
                " & & ", num2str(100*m_collect{2,3,rho_idx}.uncond_IO(1), d2), "\%", ...
                " & ", num2str(100*m_collect{2,2,rho_idx}.uncond_IO(1), d2), "\%", ...
                " & & ", num2str(targets(7), d2), "\%", " \\ ", ...
                " \hspace{0.2in} 90th percentile  & ", num2str(100*m_collect{1,1,rho_idx}.uncond_IO(9), d2), "\%", ...
                " & ", num2str(100*m_collect{1,2,rho_idx}.uncond_IO(9), d2), "\%", ...
                " & & ", num2str(100*m_collect{2,3,rho_idx}.uncond_IO(9), d2), "\%", ...
                " & ", num2str(100*m_collect{2,2,rho_idx}.uncond_IO(9), d2), "\%", ...
                " & & ", num2str(targets(8), d2), "\%", " \\ ", "\hline \end{tabular}");
              


    writematrix(s1, "Output/LaTeX_Output/vanishing_QCU_table_" + ...
                        string(false) + ".txt", 'QuoteStrings', false);

    


    g_plot = zeros(size(m_collect, 1), size(m_collect, 2), length(rho_vec)); 
    for (ii = 1:size(m_collect, 1))
        for (jj = 1:size(m_collect, 2))
            for (kk = 1:length(rho_vec))
                g_plot(ii,jj,kk) = 100*(power(1+m_collect{ii,jj,kk}.growth_rate,12)-1); 
            end
        end
    end


    close all
    figure;

    set(gcf, 'PaperUnits', 'inches');
    x_width=6.75;
    y_width=5.05;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

    rho_plot = 100*(power(1+rho_vec,12)-1); 
    subplot(2,2,1);
    p1 = plot(rho_plot, squeeze(g_plot(1, 1, :)), '--r', 'LineWidth', 2);
    hold on;
    p2 = plot(rho_plot, squeeze(g_plot(1, 2, :)), ':b', 'LineWidth', 2);
    xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
    ylabel('Growth rate, $g$', 'Interpreter', 'latex');
    title({'Capped maximum productivity increase'; ''}, 'Interpreter', 'latex');
    l1 = legend('$\bar{\lambda} = 2.0$', '$\bar{\lambda} = 2.5$', 'Interpreter', 'latex');
    xlim([0.1, 7]);
    xticks(1:7);
    set(l1, 'box', 'off');

    subplot(2,2,3);
    xplot = vary_smm(1, (1:size(smm_val, 2)) >= 14 & ~isnan(smm_val(1, :)));
    yplot = smm_val(1, (1:size(smm_val, 2)) >= 14 & ~isnan(smm_val(1, :)));
    p1 = plot(m.lambda.^xplot, yplot, '-k', 'LineWidth', 1.75);
    xlabel('Maximum productivity increase, $\bar{\lambda}$', 'Interpreter', 'latex');
    ylabel('SMM criterion', 'Interpreter', 'latex');
    xlim([m.lambda^vary_smm(1,14), m.lambda^m.s_max]);
    axes('Position',[0.32, 0.33, 0.1, 0.1]);
    box on;
    plot(m.lambda.^vary_smm(1, 1:13), smm_val(1, 1:13), 'LineStyle', '--', 'color', [0.65, 0.65, 0.65], 'LineWidth', ...
        2);

    subplot(2,2,2);
    p1 = plot(rho_plot, squeeze(g_plot(2, 1, :)), '--r', 'LineWidth', 2);
    hold on;
    p2 = plot(rho_plot, squeeze(g_plot(2, 3, :)), ':b', 'LineWidth', 2);
    p3 = plot(rho_plot, squeeze(g_plot(2, 2, :)), 'LineStyle', '-.', 'color', [0, 0.5, 0], 'LineWidth', 2);
    xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
    ylabel('Growth rate, $g$', 'Interpreter', 'latex');
    title({'Vanishing probability of'; 'quick catch-up'}, 'Interpreter', 'latex');
    l1 = legend(strcat("$\alpha = $ ", num2str(nu_vec(1,1), d2)), ...
                strcat("$\alpha = $ ", num2str(nu_vec(3,1), d2)), ...
                strcat("$\alpha = $ ", num2str(nu_vec(2,1), d2)), ...
                'Interpreter', 'latex');
    set(l1, 'box', 'off');
    xlim([0.1, 7]);
    xticks(1:7);

    subplot(2,2,4);
    p1 = plot(vary_smm(2, vary_smm(2, :) <= 0.01), smm_val(2, vary_smm(2, :) <= 0.01), '-k', 'LineWidth', 1.75);
    xlabel('Vanishing parameter, $\alpha$', 'Interpreter', 'latex');
    ylabel('SMM criterion', 'Interpreter', 'latex');
    xlim([0, 0.01]);
    axes('Position',[0.64, 0.32, 0.1, 0.1]);
    box on;
    plot(vary_smm(2, vary_smm(2, :) >= 0.01), smm_val(2, vary_smm(2, :) >= 0.01), 'LineStyle', '--', 'color', [0.65, 0.65, 0.65], 'LineWidth', ...
        2);
    xlim([0.01, 0.06]);

    saveas(gcf, "Output/Figures_Paper/vanishing_QCU_fig_" + string(false) + ".eps", 'epsc');
    saveas(gcf, "Output/Figures_Paper/vanishing_QCU_fig_" + string(false) + ".png");

end