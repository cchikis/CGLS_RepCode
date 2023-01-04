%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: wrapper_lit_review.m
% Author: Craig A. Chikis
% Date: 12/15/2021
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = wrapper_lit_review()

    cd Code/lit_review_code/LMS_replicationkit/
    calibration_EMA_submit();

    cd ../../../Code/lit_review_code/aaabk_aer2018/Baseline/
    aaabk_rep();

    cd ../../akcigitkerr_jpe2018/model/
    ak_rep();

    cd ../../peters_econometrica2019/
    peters_rep();

    cd ../../../
    own_rep_check();

    load('Output/Store_Data/aaabk_res.mat')
    aaabk = 100*save_this;

    load('Output/Store_Data/peters_res.mat')
    peters = 100*save_this;

    load('Output/Store_Data/ak_res.mat')
    ak = 100*save_this;

    load('Output/Store_Data/own_res_rep.mat')
    aa = 100*save_this(1:2, :);
    ak_at = 100*save_this([1, 3], :);

    close all
    figure;
    
	set(gcf, 'PaperUnits', 'inches');
	x_width=6.25;
	y_width=5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

    p1 = plot(aaabk(1, :), aaabk(2, :), '--r', 'LineWidth', 2);
    hold on ;
    p2 = plot(peters(1, :), peters(2, :), ':b', 'LineWidth', 2);
    p3 = plot(ak(1, :), ak(2, :), 'LineStyle', '--', 'color', [0, 0.5, 0], 'LineWidth', 2);
    p4 = plot(aa(1, :), aa(2, :), 'LineStyle', ':', 'color', [0.5, 0, 0.5], 'LineWidth', 2);
    p5 = plot(ak_at(1, :), ak_at(2, :), 'LineStyle', '--', 'color', [0.75, 0.75, 0.75], 'LineWidth', 2);
    xlabel('Discount Rate, $\rho$', 'Interpreter', 'latex');
    ylabel('Growth rate, $g$', 'Interpreter', 'latex');
    l1 = legend('Acemoglu et al. (2018)', 'Peters (2020)', 'Akcigit and Kerr (2018)',  ...
           'Acemoglu and Akcigit (2012)', 'Akcigit and Ates (2022)', 'Interpreter', 'latex');
    set(l1, 'box', 'off');

    saveas(gcf, "Output/Figures_Paper/lit_review_fig.png");
    saveas(gcf, "Output/Figures_Paper/lit_review_fig.eps", 'epsc');


end