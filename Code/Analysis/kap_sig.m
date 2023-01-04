%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: kap_sig.m
% Author: Craig A. Chikis
% Date: 01/05/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = kap_sig()


    m = benchmark_current();

    kappa_vec = [12, 9999];

    profit_vec = zeros(length(kappa_vec), 2*m.s_max+1);
    nu_s_vec = zeros(length(kappa_vec), m.s_max+1);
    markup_vec = zeros(length(kappa_vec), 2*m.s_max+1);
    firm_distributions = ones(1, m.s_max+1)/(m.s_max+1);
    quant = 0.5;
    for (ii = 1:length(kappa_vec))

        [profit_vec(ii, :), nu_s_vec(ii, :)] = get_profit(m.s_max, m.lambda, kappa_vec(ii), m.tax_rate, m.LMS_bar_s, m.two_adjust);

        [~, ~, gmarkup_all] =  markup_quant(kappa_vec(ii), m.lambda, m.s_max, firm_distributions, nu_s_vec(ii, :), quant);
        markup_vec(ii, :) = gmarkup_all;
    end




    close all
    figure;


    set(gcf, 'PaperUnits', 'inches');
    x_width=7.5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


    subplot(1,2,1);
    p1 = plot(-m.s_max:m.s_max, markup_vec(1, :), '--r', 'LineWidth', 2);
    hold on;
    p4 = plot(-m.s_max:m.s_max, markup_vec(2, :), '-k', 'LineWidth', 1.75);
    xlim([-25, 25]);
    xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
    ylabel('Gross markup, $\psi_\sigma$', 'Interpreter', 'latex');
    l1 = legend(strcat("$\kappa = $ ", num2str(kappa_vec(1))), ...
                strcat("Perfect substitution"), ...
                'Interpreter', 'latex', 'Location', 'northwest');
    set(l1, 'box', 'off');


    subplot(1,2,2);
    p1 = plot(-m.s_max:m.s_max, profit_vec(1, :), '--r', 'LineWidth', 2);
    hold on;
    p4 = plot(-m.s_max:m.s_max, profit_vec(2, :), '-k', 'LineWidth', 1.75);
    xlim([-25, 25]);
    xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
    ylabel('Opertating profit, $\pi_\sigma$', 'Interpreter', 'latex');

    saveas(gcf, "Output/Figures_Paper/prof_mark.png");
    saveas(gcf, "Output/Figures_Paper/prof_mark.eps", 'epsc');










end