%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: run_multistep.m
% Author: Craig A. Chikis
% Date: 10/16/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = run_multistep() 
    

	
    m = multistep_param(); 
    

    rho_vec = power(1 + linspace(1e-3, 0.05, 25), 1/12) - 1;

    m = workhorse_nested_robust(m);


    g_vec = zeros(1, length(rho_vec));
    markup_vec = zeros(size(g_vec));

    m_vec_iter = cell(size(rho_vec));
    for (ii = 1:length(rho_vec))
        m_iter = m;
        m_iter.rho = rho_vec(ii); 
        m_vec_iter{ii} = m_iter; 
    end

    parfor (jj = 1:length(rho_vec))

        m_vec_iter{jj} = workhorse_nested_robust(m_vec_iter{jj});

        g_vec(jj) = m_vec_iter{jj}.growth_rate;
        markup_vec(jj) = m_vec_iter{jj}.markup_wt-1; 

    end

    try 
        [m.p50_approx, ~, ~] = markup_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.5);
        [m.p90_approx, m.markup_array, ~] = markup_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.9);
    catch
        p50_approx = nan;
        p90_approx = nan;
        markup_array = nan(3, s_max+1); 
    end

    med_large = find(m.markup_array(3, :) - 0.5 >= 0, 1);
    if (~(m.s_max+1+med_large < m.l))
        jump_med = sum(m.prob(m.s_max+1-med_large+1, (m.s_max+1-med_large+1+2):end));
    else
        jump_med = 0.5*(sum(m.prob(m.s_max+1-med_large+1, (m.s_max+1-med_large+1+2):end)) + ...
                        sum(m.prob(m.s_max+1+med_large+1, (m.s_max+1+med_large+1+2):end)));
    end


    iter_this = find(-m.s_max:m.s_max < m.l); 
    collect_radprob = zeros(size(iter_this));
    for (ii = 1:length(iter_this))
        collect_radprob(ii) = sum(m.prob(iter_this(ii), (iter_this(ii)+2):end));
    end
    mu_sigma = [flip(m.firm_distributions(2:end))*0.5, m.firm_distributions(1), m.firm_distributions(2:end)*0.5];
    mu_sigma(~ismember(1:length(mu_sigma), iter_this)) = 0; 
    mu_sigma = mu_sigma/sum(mu_sigma);
    avg_prob = sum(mu_sigma(iter_this).*collect_radprob);

    savg_jump = "The estimated parameter values imply that, in an industry with the median gap, an innovating laggard has a probability of a " + ...
                "more-than-incremental advance equal to " + num2str(100*jump_med, '%0.0f') + "\%." + ...
                " For firms with $\sigma < l$, the probability of a more-than-incremental advance is, on average, about " + ...
                num2str(100*avg_prob, '%0.0f') + "\%." ; 
    writematrix(savg_jump, "Output/Figures_Paper/incremental_multistep.txt", 'QuoteStrings', false); 


    mean_collect_radprob = sum(collect_radprob); 

    prctile_inp_vec = nan; 
    m_vec_prctile = cell(1, length(prctile_inp_vec)); 
    for (ii = 1:length(prctile_inp_vec))
        m_tmp = m; 
        m_tmp.prctile_inp = prctile_inp_vec(ii);



        m_tmp = sim_wrap(m_tmp); 
        m_tmp = CSTAT_sim(m_tmp);
        m_tmp = transmat_sim(m_tmp);
        m_tmp = reg_rdsales_sim(m_tmp); 
        

        m_vec_prctile{ii} = m_tmp; 
      
    end


    try
        growth_out = 100*(power(1+m.growth_rate,12)-1);
    catch
        growth_out = nan;
    end

    try
        if (m.kappa >= 9999)
            markup_out = m.markup*100;
        else
            markup_out = (m.markup_wt-1)*100;
        end
        markup_p50_out = (m.p50_approx-1)*100;
        markup_p90_out = (m.p90_approx-1)*100;
    catch
        markup_out = nan;
        markup_p50_out = nan;
        markup_p90_out = nan;
    end


    try 
        coef_rdsales_bucket_profit_tm1 = 100*m_vec_prctile{1}.reg_rdsales_profit_tm1_bucket.Coefficients.Estimate(...
                ismember(m_vec_prctile{1}.reg_rdsales_profit_tm1_bucket.Coefficients.Row, ...
                         {'bucket_tm1_2', 'bucket_tm1_3', 'bucket_tm1_4', 'bucket_tm1_5'}) ...
                    )' ;
    catch
        coef_rdsales_bucket_profit_tm1 = nan(1,4); 
    end


    try
        profitgrowth_vol_sim = 100*m_vec_prctile{1}.std_table.mean_fun1_profitgrowth_w'; 
        rdsales_median_sim = 100*m_vec_prctile{1}.std_table.mean_fun2_rdsales_w'; 

    catch
        profitgrowth_vol_sim = nan(1,6);
        rdsales_median_sim = nan(1,6); 
    end 

 
    rng(2022+08+30, 'Threefry')
    draw_Lerner = betarnd(1.36, 8, 1e6, 1);
    draw_markup = 1./(1-draw_Lerner) - 1; 
   
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
                    
    rdsales_targets = readtable(strcat("Output/Store_Data/rdsales_t_targets_", target_filter_vec(1), ".csv")); 
    transmat_targets = readmatrix(strcat("Output/Store_Data/transmat_targets_", target_filter_vec(1), ".csv")); 
    reg_targets = readtable(strcat("Output/Store_Data/reg_targets_", target_filter_vec(1), ".csv"));
    profvol_targets = readtable(strcat("Output/Store_Data/profitvol_t_", target_filter_vec(1), ".csv")); 
    markup_targets = [mean(draw_markup), prctile(draw_markup, [50, 90])]; 
    growth_targ = 1.033; 
   
    clear draw_Lerner;
    clear draw_markup;   

    params = [m.phi_wt, m.B, m.lambda, m.zeta, m.kappa]';
    names = ["phi", "B", "lambda", "zeta", "kappa"]';
    parameterization = array2table(params, 'VariableNames', ["values"]);
    parameterization.parameter = names; 


    s1 = strcat("\begin{tabular}{lcc} \hline \hline  Moments  & Model & Data  \\  \hline  ", ...
                " Productivity growth  & ", num2str(growth_out, '%0.2f'), "\%", " & ", ...
                    num2str(growth_targ, '%0.2f'), "\%", " \\ ", ...
                " Net markup  & ", num2str(markup_out, '%0.2f'), "\%", " & ", ...
                    num2str(markup_targets(1)*100, '%0.2f'), "\%", " \\ ", ...
                "Profit volatility \\ ", ...
                " \hspace{0.2in} All firms & ", num2str(profitgrowth_vol_sim(end), '%0.2f'), "\%", " & ", ...
                    num2str(profvol_targets.sd(end)*100, '%0.2f'), "\%", " \\ ", ...
                " \hspace{0.2in} Top profit quintile & ", num2str(profitgrowth_vol_sim(1), '%0.2f'), "\%", " & ", ...
                    num2str(profvol_targets.sd(1)*100, '%0.2f'), "\%", " \\ ", ...
                " R\&D to sales \\ ", ...
                " \hspace{0.2in} All firms & ", num2str(rdsales_median_sim(end), '%0.2f'), "\%", " & ", ...
                    num2str(rdsales_targets.p50(end)*100, '%0.2f'), "\%", " \\ ", ...
                " \hspace{0.2in} Top profit quintile & ", num2str(rdsales_median_sim(1), '%0.2f'), "\%", " & ", ...
                    num2str(rdsales_targets.p50(1)*100, '%0.2f'), "\%", " \\ ", ...
                " R\&D to sales, effect relative to top profit quintile \\ ", ...
                " \hspace{0.2in} Second quintile & ", num2str(coef_rdsales_bucket_profit_tm1(1), '%0.2f'), " & ", ...
                    num2str(reg_targets.rdsales_tm1(1)*100, '%0.2f'), " \\ ", ...
                " \hspace{0.2in} Third quintile & ", num2str(coef_rdsales_bucket_profit_tm1(2), '%0.2f'), " & ", ...
                    num2str(reg_targets.rdsales_tm1(2)*100, '%0.2f'), " \\ ", ...
                " \hspace{0.2in} Fourth quintile & ", num2str(coef_rdsales_bucket_profit_tm1(3), '%0.2f'), " & ", ...
                    num2str(reg_targets.rdsales_tm1(3)*100, '%0.2f'), " \\ ", ...
                " \hspace{0.2in} Smallest quintile & ", num2str(coef_rdsales_bucket_profit_tm1(4), '%0.2f'), " & ", ...
                    num2str(reg_targets.rdsales_tm1(4)*100, '%0.2f'), " \\ ", ...
                " \hline \hline Parameters \\ ", ...
                " \hline $\phi_M$ & ", num2str(m.phi_wt, '%0.3f'), " \\ ", ...
                " $\lambda$ & ", num2str(m.lambda, '%0.3f'), " \\ ", ...
                " $B$ & ", num2str(12*(m.B*(1-m.subsidy(2))^(-m.gamma)), '%0.3f'), " \\ ", ...
                " $\kappa$ & ", num2str(m.kappa, '%0.3f'), " \\ ", ...
                " $l_M$ & ", num2str(m.l), " \\ ", ...
                "\hline \end{tabular}");

    writematrix(s1, strcat("Output/LaTeX_Output/multistep_param", num2str(m.numCPC), ".txt"), ...
                'QuoteStrings', false); 

    

    s_text = "As in the benchmark model and data, leadership is persistent in the estimated multi-step model, with " + ...
               num2str(100*m_vec_prctile{1}.transmat_out_top10(1,1), '%0.0f') + "\% of firms in the top decile remaining there the following year.";
               
    writematrix(s_text, "Output/LaTeX_Output/ms_top10" + num2str(m.numCPC) + ".txt", 'QuoteStrings', false);

    close all
    figure;

    set(gcf, 'PaperUnits', 'inches');
    x_width=3.2;
    y_width=7/2;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); 


    subplot(2,1,1);
    p1 = plot(100*(power(1+rho_vec, 12) - 1), 100*(power(1+g_vec, 12)-1), ':b', 'LineWidth', 2);
    xlim([0,5]);
    xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
    ylabel('Growth, $g$', 'Interpreter', 'latex');

    subplot(2,1,2);
    p1 = plot(100*(power(1+rho_vec, 12) - 1), 100*markup_vec, ':b', 'LineWidth', 2);
    xlim([0,5]);
    ylim([17,22]);
    xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
    ylabel('Mean net markup', 'Interpreter', 'latex');


    saveas(gcf, "Output/Figures_Paper/g_r_robustness_multistep.png");
    saveas(gcf, "Output/Figures_Paper/g_r_robustness_multistep.eps", 'epsc');




end
