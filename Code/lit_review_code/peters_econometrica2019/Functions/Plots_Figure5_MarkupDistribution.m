%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5: "The Stationary Distribution of Markups"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = Plots_Figure5_MarkupDistribution(lambda, EqObjects, Print)

%% Simulate cross-sectional distribution
NumProducts = 15000;

I = EqObjects.I;
x = EqObjects.x;
tau = EqObjects.tau;
z = tau - x;

[Gap, Names, ~] = SimulateCrossSection(I, z, tau, NumProducts);


%% Cross-sectional markup distribution
% Across products
mu_product = lambda.^Gap;
% Across firms
tmp = grpstats(mu_product.^(-1),Names,'mean');
mu_firm = tmp.^(-1);

%% Put the two distribution in bins
mu_vec_plotting_bins = [1.10; 1.2; 1.3; 1.4; 1.5; 10000];
mass_product = zeros(length(mu_vec_plotting_bins), 1);
mass_firm = zeros(length(mu_vec_plotting_bins), 1);

for j = 1:length(mu_vec_plotting_bins)
    if j == 1
    fil_p = mu_product <= mu_vec_plotting_bins(1);
    fil_f = mu_firm <= mu_vec_plotting_bins(1);
    end    
    if j > 1    
    fil_p = mu_vec_plotting_bins(j-1) < mu_product & mu_product <= mu_vec_plotting_bins(j);
    fil_f = mu_vec_plotting_bins(j-1) < mu_firm & mu_firm <= mu_vec_plotting_bins(j);
    end
    mass_product(j) =  sum(fil_p);
    mass_firm(j) =  sum(fil_f);
end

dens_product = mass_product/sum(mass_product);
dens_firm = mass_firm/sum(mass_firm);


%% Plot distributions

fig1 = figure;

barplot = bar([dens_product, dens_firm]);

barplot(1).FaceColor = [.4 .4 .4];
barplot(2).FaceColor = [.7 .7 .7];
barplot(1).LineWidth = 2;
barplot(2).LineWidth = 2;
barplot(1).BarWidth = 1;
grid on


label = {'\leq1.1','[1.1,1.2]','[1.2,1.3]','[1.3,1.4]','[1.4,1.5]','\geq1.5'};
set(gca, 'XTickLabel',label)
set(gca,'XLim',[0.5 6.5])


densylabel = ylabel('Density of markups','Interpreter','latex');
mulabel = xlabel('Markup','Interpreter','latex');


set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 12          );

set( [mulabel, densylabel]                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 18          );

FullLegend = legend( ...
    'Across products ($\mu$)','Across firms ($\mu_{f}$)');

set(FullLegend,'Interpreter','latex','FontSize',18,'Location', 'northeast','Orientation','vertical');


set(gcf, 'PaperPositionMode', 'auto');
if Print == 1
print(fig1,'-depsc','Figures/Figure5_MarkupDistribution.eps');
end


end












