%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: "The life-cycle of Markups"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [] = Plots_Figure1_LifeCycleMarkup(lambda, EqObjects, Print)

AgeGrid = [0:3:40];
AgeGrid(1) = 0.01; 

LimitMU = log(lambda) * (1 + EqObjects.I/EqObjects.tau);    % Asymptotic markup of old firms
LCmarkup = zeros(length(AgeGrid),1);    

for a = 1:length(AgeGrid)
% Expected Product Age
ap = ExpectedProductAge_FirmAge(EqObjects.x,EqObjects.tau,AgeGrid(a));
% Expected log mark-up
LCmarkup(a) = log(lambda)*(1+EqObjects.I*ap);    
end

close all
f1 = figure;

Plot_limit = plot(AgeGrid, LimitMU * ones(length(AgeGrid)));
hold on
Plot_MU_base = plot(AgeGrid, LCmarkup);
hold off
ylim([0 0.15])
grid on

%%%%%%%%%%%%%% SET STYLES

set(Plot_MU_base , ...
  'LineStyle'       , '-'           , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 9             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.7 .7 .7]    );

set(Plot_limit                           , ...
  'LineStyle'       , '--'          , ...
  'Color'           , [.7 .1 .1]    , ...
  'LineWidth'       , 3.5            );


lbl_text = text(30.5,0.116,'$\ln \lambda \times (1 + I / \tau)$','Interpreter','latex');

hXLabel = xlabel('Firm age ($a_{f}$)','Interpreter','latex');
hYLabel_mu = ylabel('Average log markup ($E[\ln \mu | a_{f}]$)','Interpreter','latex');


set(gca , 'FontName','Helvetica' ,'FontSize', 12);
set([hXLabel,hYLabel_mu],'FontName','Helvetica','FontSize', 18);
set(lbl_text,'FontName','Helvetica','FontSize', 14, 'Color', [.7 .1 .1]);
set(gcf, 'PaperPositionMode', 'auto');

if Print == 1
print(f1,'-depsc','Figures/Figure1_LifeCycleMU.eps');
end


end