%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3: "The life-cycle of Markups in Indonesia"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = Plots_Figure3_MarkupsIndonesia(Print)

mu_lifecylce = readtable("Data_MarkupLifeCycle.csv");  % ReadData

age = mu_lifecylce(:,1); age = age{:,:};
N = mu_lifecylce(:,2); N = N{:,:};
ln_mu = mu_lifecylce(:,3); ln_mu = ln_mu{:,:};
mu_H = mu_lifecylce(:,4); mu_H= mu_H{:,:};
mu_L = mu_lifecylce(:,5); mu_L= mu_L{:,:};

marker_weights = N/max(N)*700;

%%%%%%%%%%%%%% PLOT FIGURE

close all
LC_plot = figure;
plot_lnmu = plot(age, ln_mu);
hold on
plot_sdlow = plot(age, mu_L);
plot_sdhigh = plot(age, mu_H);
scatter_dots = scatter(age,ln_mu,marker_weights);
hold off
grid on

%%%%%%%%%%%%%% SET STYLES

set(plot_lnmu , ...
  'LineStyle'       , '-'           , ...
  'Marker'          , 'none'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 3             );

set([plot_sdlow,plot_sdhigh] , ...
  'LineStyle'       , '--'           , ...
  'Marker'          , 'none'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 1.5             );

set(scatter_dots , ...
  'MarkerEdgeColor' , [.3 .3 .3]    , ...
  'MarkerFaceColor' , [.3 .3 .3]    );

hXLabel = xlabel('Firm age','Interpreter','latex');
hYLabel_mu = ylabel('Average log markup','Interpreter','latex');

set(gca , 'FontName','Helvetica' ,'FontSize', 12);
set([hXLabel,hYLabel_mu],'FontName','Helvetica','FontSize', 18);

set(gcf, 'PaperPositionMode', 'auto');

if Print == 1
print(LC_plot,'-depsc','Figures/Figure3_MarkupsIndonesia.eps');
end

end
