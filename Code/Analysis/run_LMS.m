%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: run_LMS.m
% Author: Craig A. Chikis
% Date: 11/12/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = run_LMS() 
    
	


    m = LMS_param_paramcorrect();
    m.N = 50000;
    m.T = 12*12; 
    m.winsor_vec = [1, 99]; 

    rho_outer = power(1 + [0.003, 0.036], 1/12) - 1;
    m_collect_main = cell(1, length(rho_outer)); 
    for (oi = 1:length(rho_outer))
        m_tmp = m; 
        m_tmp.rho = rho_outer(oi); 

        m_tmp = workhorse_nested_robust(m_tmp);
        m_tmp = sim_wrap(m_tmp); 
        m_tmp = CSTAT_sim(m_tmp); 
        m_tmp = innovation_output(m_tmp); 
        m_tmp = FHK_sim(m_tmp); 
        m_collect_main{oi} = m_tmp; 

        [m_collect_main{oi}.p50_approx, ~, ~] =  markup_quant(m_tmp.kappa, m_tmp.lambda, m_tmp.s_max, ...
                                                              m_tmp.firm_distributions, m_tmp.nu_s, 0.5);

        [m_collect_main{oi}.p90_approx, m_collect_main{oi}.markup_array, ~] =  markup_quant(m_tmp.kappa, m_tmp.lambda, m_tmp.s_max, ...
                                                                                            m_tmp.firm_distributions, m_tmp.nu_s, 0.9);

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
    

     


    LMS_tab = "\begin{tabular}{lccccccc}  \hline\hline " + ...
              " Moments & $\rho = 0.3\%$  & $\rho = 3.6\%$ & Data \\ " + ...
              " \hline Markup \\  " + ...
              " \hspace{0.2in}  Mean  & " + num2str(100*(m_collect_main{1}.markup_wt-1), '%0.2f') + "\%" + " & "  + ...
                num2str(100*(m_collect_main{2}.markup_wt-1), '%0.2f') + "\%" + " & " + ...
                num2str(markup_targ(1), '%0.2f') + "\%" + " \\ " + ...
              " \hspace{0.2in}  50th percentile  & " + num2str(100*(m_collect_main{1}.p50_approx-1), '%0.2f') + "\%" + " & "  + ...
                num2str(100*(m_collect_main{2}.p50_approx-1), '%0.2f') + "\%" + " & " + ...
                num2str(markup_targ(2), '%0.2f') + "\%" + " \\ " + ...
              " \hspace{0.2in}  90th percentile  & " + num2str(100*(m_collect_main{1}.p90_approx-1), '%0.2f') + "\%" + " & "  + ...
                num2str(100*(m_collect_main{2}.p90_approx-1), '%0.2f') + "\%" + " & " + ...
                num2str(markup_targ(3), '%0.2f') + "\%" + " \\ " + ...
              " Profit volatility \\ " + ...
              " \hspace{0.2in}  All firms  & " + num2str(100*m_collect_main{1}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + "\%" + " & "  + ...
                num2str(100*m_collect_main{2}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + "\%" + " & " + ...
                num2str(profvol_targets.sd(end)*100, '%0.2f') + "\%" + " \\ " + ...
              " \hspace{0.2in}  Top profit quintile  & " + num2str(100*m_collect_main{1}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + "\%" + " & "  + ...
                num2str(100*m_collect_main{2}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + "\%" + " & " + ...
                num2str(profvol_targets.sd(1)*100, '%0.2f') + "\%" + " \\ " + ...
              " R\&D to sales \\ " + ...
              " \hspace{0.2in}  All firms  & " + num2str(100*m_collect_main{1}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + "\%" + " & "  + ...
                num2str(100*m_collect_main{2}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + "\%" + " & " + ...
                num2str(100*rdsales_targets.p50(end), '%0.2f') + "\%" + " \\ " + ...
              " \hspace{0.2in}  Top profit quintile  & " + num2str(100*m_collect_main{1}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + "\%" + " & "  + ...
                num2str(100*m_collect_main{2}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + "\%" + " & " + ...
                num2str(100*rdsales_targets.p50(1), '%0.2f') + "\%" + " \\ " + ...
              " Innovation output \\ " + ...
              " \hspace{0.2in}  Mean  & " + num2str(100*m_collect_main{1}.uncond_IO(1), '%0.2f') + "\%" + " & "  + ...
                num2str(100*m_collect_main{2}.uncond_IO(1), '%0.2f') + "\%" + " & " + ...
                num2str(uncond_IO_targets(1), '%0.2f') + "\%" + " \\ " + ...
              " \hspace{0.2in}  90th percentile  & " + num2str(100*m_collect_main{1}.uncond_IO(9), '%0.2f') + "\%" + " & "  + ...
                num2str(100*m_collect_main{2}.uncond_IO(9), '%0.2f') + "\%" + " & " + ...
                num2str(uncond_IO_targets(9), '%0.2f') + "\%" + " \\ " + ...
                " FHK within  & " + num2str(100*m_collect_main{1}.WITHIN_out, '%0.2f') + "\%" + " & "  + ...
                num2str(100*m_collect_main{2}.WITHIN_out, '%0.2f') + "\%" + " & " + ...
                num2str(FHK_targ(1), '%0.2f') + "\%" + " \\ " + ...
                " \hline \end{tabular}";

    writematrix(LMS_tab, 'Output/LaTeX_Output/LMS_tab.txt', 'QuoteStrings', false); 

                
                

   


end