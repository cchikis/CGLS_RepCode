%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: transition_solve.m
% Author: Craig A. Chikis
% Date: 09/12/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = transition_solve() 



    m = benchmark_current();


    years = 500;
    yearsrhochange_vec = [1/12, 25, 25];
    no_surprise_vec = [true, true, false]; 
    rhoend = 0.01;
    rhostart = 0.03;
    justprodlabor = true;
    res_cell = cell(1, length(yearsrhochange_vec) + 1);
    feedmu = nan;
    for (ii = 1:length(yearsrhochange_vec))
        yearsrhochange = yearsrhochange_vec(ii); 
        no_surprise = no_surprise_vec(ii); 

        [prod, quality, dL_t_dt, ...
        muder, xder, rdlabor, prodlabor, ...
        prodlaborder, markup_vec, ...
        rho_path, tplot, ximpact, ...
        muimpact, ximpactback, muimpactback, tplotback, ...
        quality2, ...
        prod2, prod3] = ...
                    transition_wrapper(m, rhostart, rhoend, years, yearsrhochange, justprodlabor, feedmu, no_surprise);

        res_cell{ii} = {prod, quality, dL_t_dt, ...
                            muder, xder, rdlabor, prodlabor, ...
                            prodlaborder, markup_vec, ...
                            rho_path, tplot, ximpact, ...
                            muimpact, ximpactback, muimpactback, tplotback, ...
                            quality2, ...
                            prod2, ...
                            prod3};
    end



    close all 
    figure;

        
    set(gcf, 'PaperUnits', 'inches');
    x_width=8;
    y_width=0.7*3.5;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

    subplot(1,3,1);
    p1 = plot([-5:-1, res_cell{1}{11}(1:(end-1))], ...
              [res_cell{1}{10}(1)*ones(1,5), res_cell{1}{10}(2:end)],'LineStyle', '-', 'LineWidth', 2, ...
            'color', [0.4, 0.4, 0.4]);
    hold on ;
    p2 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], ...
              [res_cell{2}{10}(1)*ones(1,5), res_cell{2}{10}(2:end)], '--r', 'LineWidth', 2) ;
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Discount rate'; ''}, 'Interpreter', 'latex');
    xlim([-5,50]);

    subplot(1,3,2);
    p1 = plot([-5:-1, res_cell{1}{11}(1:(end-1))], ...
              100*(power(1+[res_cell{1}{1}(1)*ones(1,5), res_cell{1}{1}(2:end)],12)-1), 'LineStyle', ...
            '-', 'LineWidth', 2, 'color', [0.4, 0.4, 0.4]);
    hold on ;
    p2 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], ...
             100*(power(1+[res_cell{2}{1}(1)*ones(1,5), res_cell{2}{1}(2:end)],12)-1), ...
            '--r', 'LineWidth', 2);
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Productivity growth'; ''}, 'Interpreter', 'latex');
    xlim([-5,50]);
    ylim([0.9, 1.15]);

    subplot(1,3,3);
    p1 = plot([-5:-1, res_cell{1}{11}(1:(end-1))], ...
              (power(1+[res_cell{1}{9}(2,1)*ones(1,5), res_cell{1}{9}(2,2:end)],1)-1), 'LineStyle', ...
            '-', 'LineWidth', 2, 'color', [0.4, 0.4, 0.4]);
    hold on ;
    p2 = plot([-5:-1, res_cell{1}{11}(1:(end-1))], ...
              (power(1+[res_cell{1}{9}(3,1)*ones(1,5), res_cell{1}{9}(3,2:end)],1)-1), 'LineStyle', ...
            '--', 'LineWidth', 2, 'color', [0.4, 0.4, 0.4]);
    p3 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], ...
              (power(1+[res_cell{2}{9}(2,1)*ones(1,5), res_cell{2}{9}(2,2:end)],1)-1), '-r', 'LineWidth', 2);
    p4 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], ...
              (power(1+[res_cell{2}{9}(3,1)*ones(1,5), res_cell{2}{9}(3,2:end)],1)-1), '--r', 'LineWidth', 2);
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Net markups'; ''}, 'Interpreter', 'latex');
    xlim([-5,50]);
    text(10,37.5,'90th percentile', 'Interpreter', 'latex');
    text(10,16,'Median', 'Interpreter', 'latex');


    saveas(gcf, "Output/Figures_Paper/transition_1x3_v2.png");
    saveas(gcf, "Output/Figures_Paper/transition_1x3_v2.eps", 'epsc');




    close all 
    figure;

        
    set(gcf, 'PaperUnits', 'inches');
    x_width=8;
    y_width=0.7*3.5;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

    subplot(1,3,1);
    p3 = plot([-5:-1, res_cell{3}{11}(1:(end-1))], [res_cell{3}{10}(1)*ones(1,5), res_cell{3}{10}(2:end)], ...
            'LineStyle', ':', 'color', [0, 0.5, 0], 'LineWidth', 2);
    hold on ;
    p4 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], [res_cell{2}{10}(1)*ones(1,5), res_cell{2}{10}(2:end)], ...
            'LineStyle', '--', 'color', [1, 0, 0], 'LineWidth', 2);
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Discount rate'; ''}, 'Interpreter', 'latex');
    xlim([-5,50]);

    subplot(1,3,2);
    p3 = plot([-5:-1, res_cell{3}{11}(1:(end-1))], 100*(power(1+[res_cell{3}{1}(1)*ones(1,5), res_cell{3}{1}(2:end)],12)-1), ...
            'LineStyle', ':', 'color', ...
            [0, 0.5, 0], 'LineWidth', 2);
    hold on ;
    p4 =  plot([-5:-1, res_cell{2}{11}(1:(end-1))], 100*(power(1+[res_cell{2}{1}(1)*ones(1,5), res_cell{2}{1}(2:end)],12)-1), ...
            'LineStyle', '--', 'color', ...
            [1, 0, 0], 'LineWidth', 2);
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Productivity growth'; ''}, 'Interpreter', 'latex');
    xlim([-5,50]);
    ylim([0.9, 1.15]);
    l1 = legend([p4,p3], {"One-time" + newline + "($t=0$) surprise", "Continuously" + newline + "surprised"}, 'Interpreter', 'latex', 'Location', 'southeast');
    set(l1, 'box', 'off');
    l1.ItemTokenSize = l1.ItemTokenSize*0.5;

    subplot(1,3,3);
    p6 = plot([-5:-1, res_cell{3}{11}(1:(end-1))], (power(1+[res_cell{3}{9}(1,1)*ones(1,5), res_cell{3}{9}(1,2:end)],1)-1), 'LineStyle', ':', 'color', [0,0.5,0], ...
            'LineWidth', 2);
    hold on ;
    p7 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], (power(1+[res_cell{2}{9}(1,1)*ones(1,5), res_cell{2}{9}(1,2:end)],1)-1), 'LineStyle', '--', 'color', [1,0,0], ...
             'LineWidth', 2);
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Net markups'; ''}, 'Interpreter', 'latex');
    xlim([-5,50]);


    saveas(gcf, "Output/Figures_Paper/transition_1x3_v3.png");
    saveas(gcf, "Output/Figures_Paper/transition_1x3_v3.eps", 'epsc');







end


