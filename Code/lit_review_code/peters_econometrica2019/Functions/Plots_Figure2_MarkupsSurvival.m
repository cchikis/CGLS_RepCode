%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: "Life-cycle Dynamics of Markups and Survival"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Plots_Figure2_MarkupsSurvival(lambda, EqObjects, Print)

NCF = 8;
AgeGrid = [0:1:NCF];
AgeGrid(1) = 0.01; 


I = EqObjects.I; tau = EqObjects.tau; x = 0.3*EqObjects.tau;


LCmarkupCF = zeros(NCF,3);
Survivors = zeros(NCF,3);



for a = 1:NCF+1
    
Age = AgeGrid(a);   

%% MARKUPS
% (1) Baseline profile    
ap = ExpectedProductAge_FirmAge(x,tau,Age);
LCmarkupCF(a,1) = log(lambda)*(1+I*ap);  
% (2) High x
x_high = 1.5*EqObjects.tau;
ap = ExpectedProductAge_FirmAge(x_high,tau,Age);
LCmarkupCF(a,2) = log(lambda)*(1+I*ap);  
% (3) High tau
tau_high = 10*tau; 
ap = ExpectedProductAge_FirmAge(x,tau_high,Age);
LCmarkupCF(a,3) = log(lambda)*(1+I*ap);  

%% SURVIVAL

% (1) Baseline profile 
Survivors(a,1) = (tau - x) * exp(-(tau - x)*Age) / (tau - x * exp(-(tau-x)*Age));
% (2) High x
x_high = 1.5*EqObjects.tau;
Survivors(a,2) = (tau - x_high) * exp(-(tau - x_high)*Age) / (tau - x_high * exp(-(tau-x_high)*Age));
% (3) High tau
tau_high = 3*tau; 
Survivors(a,3) = (tau_high - x) * exp(-(tau_high - x)*Age) / (tau_high - x * exp(-(tau_high-x)*Age));

end

%%%%%%%%%%%%%% PLOT LIFECYCLE
f1 = figure;
% Markups (left axis)
Plot_MU_base = plot(AgeGrid, LCmarkupCF(:,1));
hold on
Plot_MU_x = plot(AgeGrid, LCmarkupCF(:,2));
Plot_MU_tau = plot(AgeGrid, LCmarkupCF(:,3));
hold off
grid on

%%%%%%%%%%%%%% PLOT SURVIVORS
f2 = figure;
Plot_Survivors_base = plot(AgeGrid, Survivors(:,1));
hold on
Plot_Survivors_x = plot(AgeGrid, Survivors(:,2));
Plot_Survivors_tau = plot(AgeGrid, Survivors(:,3));
hold off
ylim([0 1.1])
yticks([.2 0.4 0.6 .8 1])
grid on

%%%%%%%%%%%%%% SET STYLES


set([Plot_MU_base,Plot_Survivors_base] , ...
  'LineStyle'       , '-'           , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.7 .7 .7]    );

set([Plot_MU_x,Plot_Survivors_x]                           , ...
  'LineStyle'       , '-'          , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.7 .1 .1]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.7 .1 .1]    );

set([Plot_MU_tau, Plot_Survivors_tau ]                          , ...
  'LineStyle'       , '-'          , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.1 .1 .9]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.1 .1 .9]    );


figure(f1);
hXLabel = xlabel('Firm age ($a_{f}$)','Interpreter','latex');
title_lbl = title('Average log markup ($E[\ln \mu | a_{f}]$)','Interpreter','latex');

lbl_xh = text(6.8,0.105,'Higher $x$','Interpreter','latex');
lbl_tauh = text(6.8,0.14,'Higher $\tau$','Interpreter','latex');
lbl_base = text(6.8,0.17,'Baseline','Interpreter','latex');

set(gca , 'FontName','Helvetica' ,'FontSize', 12);
set([hXLabel, title_lbl],'FontName','Helvetica','FontSize', 18);
set([lbl_xh,lbl_tauh,lbl_base],'FontName','Helvetica','FontSize', 16);
set(lbl_tauh,'Color',[.1 .1 .9]);
set(lbl_xh,'Color',[.7 .1 .1]);
set(lbl_base,'Color',[.3 .3 .3]);

annotation('arrow', [.75 .75], [.75 .665],'Color', [.1 .1 .9],'LineWidth', 2.5 );
annotation('arrow', [.765 .765], [.75 .58],'Color', [.7 .1 .1],'LineWidth', 2.5 );

figure(f2);

hXLabel = xlabel('Firm age ($a_{f}$)','Interpreter','latex');
title_lbl = title('Survival rate ($S(a_{f})$)','Interpreter','latex');


lbl_xh = text(6.1,0.54,'Higher $x$','Interpreter','latex');
lbl_tauh = text(6.1,0.06,'Higher $\tau$','Interpreter','latex');
lbl_base = text(6.1,0.33,'Baseline','Interpreter','latex');


set(gca , 'FontName','Helvetica' ,'FontSize', 12);
set([hXLabel,title_lbl],'FontName','Helvetica','FontSize', 18);

set([lbl_xh,lbl_tauh,lbl_base],'FontName','Helvetica','FontSize', 16);
set(lbl_tauh,'Color',[.1 .1 .9]);
set(lbl_xh,'Color',[.7 .1 .1]);
set(lbl_base,'Color',[.3 .3 .3]);

annotation('arrow', [.67 .67], [.33 .145],'Color', [.1 .1 .9],'LineWidth', 2.5 );
annotation('arrow', [.67 .67], [.36 .469],'Color', [.7 .1 .1],'LineWidth', 2.5 );


set(gcf, 'PaperPositionMode', 'auto');
if Print == 1
print(f1,'-depsc','Figures/Figure2_Markup.eps');
print(f2,'-depsc','Figures/Figure2_Survival.eps');
end



end