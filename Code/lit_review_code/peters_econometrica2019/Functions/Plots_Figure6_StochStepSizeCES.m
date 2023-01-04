%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6: "Misallocation and Markups: Stochastic Step Size and CES Preferences"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [] = Plots_Figure6_StochStepSizeCES(Print)

load('MisallocationStochStepSize.mat')
load('MisallocationCES.mat')

f1 = figure;

gamma_M_plot = plot(MisallocationStochStepSize.kappa, MisallocationStochStepSize.M);

hold on
gamma_Mp_plot = plot(MisallocationStochStepSize.kappa, MisallocationStochStepSize.Mp);
hold off
ylim([0.85 1])
xlim([0 0.85])
yticks([0.85 0.9 0.95 1])
xticks([0.2 0.4 0.6 0.8])
grid on

f2 = figure;

sigma_M_plot = plot(MisallocationCES.sigma, MisallocationCES.M);
hold on
sigma_Mnaive_plot = plot(MisallocationCES.sigma, MisallocationCES.Mnaive);
hold off
ylim([0.95 1])
yticks([.95 .96 .97 .98 .99 1])
grid on



set([gamma_M_plot,sigma_M_plot]                        , ...
  'LineStyle'       , '-'           , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 9             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.7 .7 .7]    );

set([gamma_Mp_plot,sigma_Mnaive_plot]                        , ...
  'LineStyle'       , '--'           , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 2             , ...
  'MarkerSize'      , 9             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.9 .9 .9]    );

figure(f1)
hXLabel = xlabel('$\kappa$','Interpreter','latex');
hYLabelM = ylabel('TFP Misallocation ($\mathcal{M}$)','Interpreter','latex');

lbl_text = text(0.6,0.985,'Calibrated Model','Interpreter','latex');
set(lbl_text,'FontName','Helvetica','FontSize', 14, 'Color', [.3 .3 .3]);

lbl2_text = text(0.57,0.89,'Holding $\vartheta_{I}$ fixed','Interpreter','latex');
set(lbl2_text,'FontName','Helvetica','FontSize', 14, 'Color', [.3 .3 .3]);

figure(f2)
hXLabel = xlabel('$\sigma$','Interpreter','latex');
hYLabelM = ylabel('TFP Misallocation ($\mathcal{M}$)','Interpreter','latex');

lbl_text = text(3.5,0.995,'$\mathcal{M}$','Interpreter','latex');
set(lbl_text,'FontName','Helvetica','FontSize', 14, 'Color', [.3 .3 .3]);

lbl2_text = text(3.5,0.983,'$\mathcal{M}^{N}$','Interpreter','latex');
set(lbl2_text,'FontName','Helvetica','FontSize', 14, 'Color', [.3 .3 .3]);


set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 12          );

set( [hXLabel,hYLabelM]                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 18          );


set(gcf, 'PaperPositionMode', 'auto');
if Print == 1
print(f1,'-depsc','Figures/Figure6_StochasticStepSize.eps');
print(f2,'-depsc','Figures/Figure6_CES.eps');
end


end