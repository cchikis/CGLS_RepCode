%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: error_analysis.m
% Author: Craig A. Chikis
% Date: 01/06/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = error_analysis()


    models = {@benchmark_current, @entry_current, @LMS_param_paramcorrect};
    rho_vec = power(1 + linspace(1e-3, 0.05, 25), 1/12)-1;

    res_vec = cell(length(models), length(rho_vec));
    for (ii = 1:length(models))

        m = models{ii}(); 

        res_iter = cell(1, length(rho_vec));
        for (jj = 1:length(rho_vec))
            m_tmp = m;
            m_tmp.rho = rho_vec(jj); 
            res_iter{jj} = m_tmp; 
        end
        parfor (jj = 1:length(rho_vec))

           res_iter{jj} = workhorse_nested_robust(res_iter{jj});

        end
        res_vec(ii, :) = res_iter;
    end

    rho_plot = 100*(power(1+rho_vec, 12)-1);
    vfi_error_vec = zeros(size(res_vec));
    vfi_hjb_error_vec = zeros(size(res_vec));
    end_hjb_error_vec = zeros(size(res_vec));
    vfi_error_rho_vec = zeros(size(res_vec));
    vfi_hjb_error_rho_vec = zeros(size(res_vec));
    end_hjb_error_rho_vec = zeros(size(res_vec));
    for (ii = 1:size(res_vec, 1))
        for (jj = 1:size(res_vec, 2))
            vfi_error_vec(ii,jj) = res_vec{ii,jj}.vfi_error;
            vfi_hjb_error_vec(ii,jj) = res_vec{ii,jj}.vfi_hjb_error;
            end_hjb_error_vec(ii,jj) = res_vec{ii,jj}.end_hjb_error;
            vfi_error_rho_vec(ii,jj) = res_vec{ii,jj}.vfi_error_rho;
            vfi_hjb_error_rho_vec(ii,jj) = res_vec{ii,jj}.vfi_hjb_error_rho;
            end_hjb_error_rho_vec(ii,jj) = res_vec{ii,jj}.end_hjb_error_rho;
        end
    end

    close all
    figure; 


    set(gcf, 'PaperUnits', 'inches');
    x_width=6.75*1.5;
    y_width=5.05/2;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


    subplot(1,3,1);
    plot(rho_plot, log(end_hjb_error_rho_vec(1, :))/log(10), '-k', 'LineWidth', 1.75);
    xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
    ylabel('$\log_{10} RMSE_{HJB}$', 'Interpreter', 'latex');
    ylim([-20, -10]);
    xlim([0, 5]);
    title({'Benchmark'; ''}, 'Interpreter', 'latex');

    subplot(1,3,2);
    plot(rho_plot, log(end_hjb_error_rho_vec(2, :))/log(10), '-k', 'LineWidth', 1.75);
    xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
    ylabel('$\log_{10} RMSE_{HJB}$', 'Interpreter', 'latex');
    ylim([-20, -10]);
    xlim([0, 5]);
    title({'Entry'; ''}, 'Interpreter', 'latex');


    subplot(1,3,3);
    plot(rho_plot, log(end_hjb_error_rho_vec(3, :))/log(10), '-k', 'LineWidth', 1.75);
    xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
    ylabel('$\log_{10} RMSE_{HJB}$', 'Interpreter', 'latex');
    xlim([0, 5]);
    ylim([-20, -10]);
    title({'LMS (2022)'; ''}, 'Interpreter', 'latex');


    saveas(gcf, "Output/Figures_Paper/err6.png");
    saveas(gcf, "Output/Figures_Paper/err6.eps", 'epsc');

end