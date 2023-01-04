%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: SCUexplore.m
% Author: Craig A. Chikis
% Date: 09/07/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = SCUexplore()


       s_max = 150; 
       eta_cell = {[0, (0.024/12)*ones(1,s_max)], [0, 0.5*0.024/12*ones(1, s_max)]};
       phi_wt_exog_vec =[1, 0.5]; 



       rho_vec = sort([power(1 + linspace(1e-3, 0.05, 25), 1/12) - 1, power(1.015, 1/12)-1]);
       output_cell = cell(length(eta_cell), length(rho_vec));
       for (jj = 1:length(eta_cell))
              eta = eta_cell{jj};
              phi_wt_exog = phi_wt_exog_vec(jj); 



              for (ii = 1:length(rho_vec)) 
                     m = AA2012_scu(s_max, phi_wt_exog, eta);
                     rho = rho_vec(ii);
                     m.rho = rho; 
                     m = workhorse_nested_robust(m);


                     output_cell{jj, ii} = m; 
              end


       end

       G = @(subsidy, wage_share, x, B, gamma) (1-subsidy).*wage_share.*(x./B).^(1./gamma);

       g_vec = zeros(size(output_cell));
       markup_vec = zeros(size(output_cell)); 
       w_vec = zeros(size(g_vec));
       g_vec_parts = zeros(size(output_cell, 1), s_max+1, size(output_cell, 2));  
       mu_vec = zeros(size(g_vec_parts));
       x_vec = zeros(size(output_cell, 1), 2*s_max+1, length(output_cell));
       lab_demand = zeros(size(output_cell, 1), 2, 2*s_max+1, length(output_cell));
       for(jj = 1:size(output_cell, 1))
              for (ii = 1:size(output_cell, 2))
                     g_vec(jj, ii) = 100*(power(1+output_cell{jj, ii}.growth_rate,12)-1);
                     markup_vec(jj, ii) = 100*output_cell{jj, ii}.markup;
                     w_vec(jj, ii) = output_cell{jj, ii}.wage_share;
                     g_vec_parts(jj, :, ii) = 12*log(m.lambda)*output_cell{jj,ii}.firm_distributions.*[2, ones(1,s_max)].*...
                                                        output_cell{jj,ii}.investment((m.s_max+1):end); 
                     mu_vec(jj, :, ii) = output_cell{jj, ii}.firm_distributions;
                     x_vec(jj, :, ii) = output_cell{jj, ii}.investment; 
                     lab_demand(jj, 1, :, ii) =  G(m.subsidy, w_vec(jj, ii), x_vec(jj, :, ii), m.B, m.gamma);
                     lab_demand(jj, 2, :, ii) = [zeros(1,m.s_max), 0.5/(w_vec(jj, ii)*m.lambda^min(0, m.LMS_bar_s)), ...
                                                        1./(w_vec(jj, ii)*m.lambda.^min(m.LMS_bar_s, (1:m.s_max)))];
              end
       end
       rho_plot = 100*(power(1+rho_vec, 12)-1); 
       idx = [find(abs(rho_plot - 0.1) == min(abs(rho_plot - 0.1)), 1), ...
              find(abs(rho_plot - 1.5) == min(abs(rho_plot - 1.5)) ,1)]; 




       close all
       figure;

       set(gcf, 'PaperUnits', 'inches');
       x_width=6.5;
       y_width=3;
       set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


       subplot(1,2,1);
       p1 = plot(rho_plot, g_vec(1, :), '-k', 'LineWidth', 2);
       xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
       ylabel('Growth rate, $g$', 'Interpreter', 'latex');
       title({'Acemoglu and Akcigit (2012) Slow Catch-up'; ''}, 'Interpreter', 'latex');
       xlim([0,5]);

       subplot(1,2,2);
       p1 = plot(rho_plot, g_vec(2, :), '--r', 'LineWidth', 2);
       xlabel('Discount rate, $\rho$', 'Interpreter', 'latex');
       ylabel('Growth rate, $g$', 'Interpreter', 'latex');
       title({'Version with severe constraints'; 'on creative destruction'}, 'Interpreter', 'latex');
       xlim([0,5]);

       saveas(gcf, "Output/Figures_Paper/AASCU1.png");
       saveas(gcf, "Output/Figures_Paper/AASCU1.eps", 'epsc');

       close all

       idx1 = find(abs(rho_plot - 1.5) == min(abs(rho_plot - 1.5)));
       idx2 = find(abs(rho_plot - 0.1) == min(abs(rho_plot - 0.1)));

       set(gcf, 'PaperUnits', 'inches');
       x_width=6.5;
       y_width=3*1.5;
       set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

       subplot(2,2,1);
       p1 = plot(0:m.s_max, x_vec(2, (m.s_max+1):end, idx1)*12, ':b', 'LineWidth', 2);
       hold on ;
       p2 = plot(0:m.s_max, flip(x_vec(2, 1:(m.s_max+1), idx1)*12), 'LineStyle', '-.', 'color', [0,0.5,0], 'LineWidth', 2);
       xlabel('Technology gap, $s$', 'Interpreter', 'latex');
       ylabel('$x_{s}$', 'Interpreter', 'latex');
       title({strcat("Innovation rates (", "$\rho = $", num2str(rho_plot(idx1)), "\%)"); ""}, 'Interpreter', 'latex');
       ylim([0,6]);
       xlim([0,25]);
       l1 = legend('Leader', 'Laggard', 'Interpreter', 'latex');
       set(l1, 'box', 'off');

       subplot(2,2,2);
       p1 = plot(0:m.s_max, x_vec(2, (m.s_max+1):end, idx2)*12, '--r', 'LineWidth', 2);
       hold on ;
       p2 = plot(0:m.s_max, flip(x_vec(2, 1:(m.s_max+1), idx2)*12), 'LineStyle', '-', 'color', [0.5,0.5,0.5], 'LineWidth', 2);
       xlabel('Technology gap, $s$', 'Interpreter', 'latex');
       ylabel('$x_{s}$', 'Interpreter', 'latex');
       title({strcat("Innovation rates (", "$\rho = $", num2str(rho_plot(idx2)), "\%)"); ""}, 'Interpreter', 'latex');
       ylim([0,6]);
       xlim([0, 25]);
       l1 = legend('Leader', 'Laggard', 'Interpreter', 'latex');
       set(l1, 'box', 'off');

       subplot(2,2,3);
       p1 = plot(-m.s_max:m.s_max, 12*(x_vec(2, :, idx2) - x_vec(2, :, idx1)), '-k', 'LineWidth', 2);
       xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
       ylabel('$x_\sigma(\rho = 0.1\%) - x_\sigma(\rho = 1.5\%)$', 'Interpreter', 'latex');
       title({'Change in innovation rate'; ''}, 'Interpreter', 'latex');
       xlim([-25, 25]);

       subplot(2,2,4);
       p1 = plot(0:m.s_max, 100*(mu_vec(2, :, idx2) - mu_vec(2, :, idx1)), '-k', 'LineWidth', 2);
       xlabel('Technology gap, $s$', 'Interpreter', 'latex');
       ylabel('$\mu_s(\rho = 0.1\%) - \mu_s(\rho = 1.5\%)$ (\%)', 'Interpreter', 'latex');
       title({'Change in gap distribution'; ''}, 'Interpreter', 'latex');
       xlim([0, 10]);

       saveas(gcf, "Output/Figures_Paper/AASCU2.png");
       saveas(gcf, "Output/Figures_Paper/AASCU2.eps", 'epsc');



end