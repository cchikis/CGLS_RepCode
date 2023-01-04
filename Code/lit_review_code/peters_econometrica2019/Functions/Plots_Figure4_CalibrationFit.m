%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4: "Non-targeted Moments: Model vs Data"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [] = Plots_Figure4_CalibrationFit(lambda, EqObjects, Print)

N = 10;
Age_vec = linspace(0,9,N);
lnmuAge = zeros(N,1);
lnEmplAge = zeros(N,1);
Survivors = zeros(N,1);


I = EqObjects.I; tau = EqObjects.tau; x = EqObjects.x;


%% Get data

Data = readtable("Data_IndonesiaPanel.csv");
Data_lnmu = Data(:,2); Data_lnmu = Data_lnmu{:,:};
Data_lnl = Data(:,3); Data_lnl = Data_lnl{:,:};

Data = readtable("Data_IndonesiaSurvival.csv");
Data_year = Data(:,4); Data_year = Data_year{:,:};
Data_Survivors = Data(:,1); Data_Survivors = Data_Survivors{:,:};
Data_Cohort = Data(:,2); Data_Cohort = Data_Cohort{:,:};

Filter = Data_year == 1991 + 9;
Data_SurvivorsCrossSection  = Data_Survivors(Filter);
Data_SurvivorsPanel = Data_Survivors(Data_Cohort == 1991);



for j = 1:N   
   
   %% Lifecycle of markups 
   lnmuAge(j) = log(lambda) * I * (ExpectedProductAge_FirmAge(x,tau,Age_vec(j)+0.5) - ExpectedProductAge_FirmAge(x,tau,0.5)); 
    
   %% Lifecycle of Employment
   Nvec = linspace(1,200,200)';
   % Number of products at age a
   lnn_age = (1-GammaFunction(Age_vec(j)+0.5,x,tau))/GammaFunction(Age_vec(j)+0.5,x,tau) * sum(log(Nvec).*GammaFunction(Age_vec(j)+0.5,x,tau).^(Nvec)); 
   % Number of products at entry
   lnn_entry = (1-GammaFunction(0.5,x,tau))/GammaFunction(0.5,x,tau) * sum(log(Nvec).*GammaFunction(0.5,x,tau).^(Nvec));
   % Employment  = number of products - markups 
   lnEmplAge(j) = lnn_age - lnn_entry - lnmuAge(j); 
   
   %% Survival rates
   Survivors(j) = (tau - x) * exp(-(tau - x)*Age_vec(j)) / (tau - x * exp(-(tau-x)*Age_vec(j)));
   
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MARKUPS
f1 = figure;

Results_Model_lnmu = plot(Age_vec, lnmuAge);
hold all
Results_Data_lnmu = plot(Age_vec, Data_lnmu);
grid on
ylim([-0.01,0.1])
yticks([0 0.02 0.04 0.06 0.08 0.1])


set(Results_Model_lnmu                        , ...
  'LineStyle'       , '-'           , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.7 .7 .7]    );

set(Results_Data_lnmu                            , ...
  'LineStyle'       , '--'          , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.7 .1 .1]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.7 .1 .1]    , ...
  'MarkerFaceColor' , [.7 .1 .1]    );

hXLabel = xlabel('Age' ,'Interpreter','latex');
title_lbl = title('log markup' ,'Interpreter','latex');

label = {'0','1','2','3','4','5','6','7','8','9'};
set(gca, 'XTickLabel',label)

set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 12          );

set( [hXLabel,title_lbl]                      , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 18          );


lbl_data = text(8,0.062,'Data','Interpreter','latex');
lbl_model = text(8,0.091,'Model','Interpreter','latex');

set([lbl_data,lbl_model],'FontName','Helvetica','FontSize', 16);
set(lbl_data,'Color',[.7 .1 .1]);
set(lbl_model,'Color',[.3 .3 .3]);


set(gcf, 'PaperPositionMode', 'auto');
if Print == 1
print(f1,'-depsc','Figures/Figure4_MarkupFit.eps');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EMPLOYMENT
f2 = figure;

Results_Model_lnEmpl = plot(Age_vec, lnEmplAge);
hold all
Results_Data_lnEmpl = plot(Age_vec, Data_lnl);
grid on

set(Results_Model_lnEmpl                        , ...
  'LineStyle'       , '-'           , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.7 .7 .7]    );

set(Results_Data_lnEmpl                            , ...
  'LineStyle'       , '--'          , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.7 .1 .1]   , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.7 .1 .1]    , ...
  'MarkerFaceColor' , [.7 .1 .1]    );

hXLabel = xlabel('Age','Interpreter','latex');
title_lbl = title('log employment','Interpreter','latex');

label = {'0','1','2','3','4','5','6','7','8','9'};
set(gca, 'XTickLabel',label)

set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 12          );

set( [hXLabel,title_lbl]                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 18          );


lbl_data = text(8,0.65,'Data','Interpreter','latex');
lbl_model = text(8,0.5,'Model','Interpreter','latex');

set([lbl_data,lbl_model],'FontName','Helvetica','FontSize', 16);
set(lbl_data,'Color',[.7 .1 .1]);
set(lbl_model,'Color',[.3 .3 .3]);


set(gcf, 'PaperPositionMode', 'auto');

if Print == 1
print(f2,'-depsc','Figures/Figure4_EmploymentFit.eps');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SURVIVAL

f3 = figure;

plotModel = plot(Age_vec, Survivors);
hold on
plotDataCS = plot(Age_vec, Data_SurvivorsCrossSection);
plotDataPanel = plot(Age_vec, Data_SurvivorsPanel);
hold off
grid on

axis([-Inf, Inf, 0,  1.1]);

set(plotModel                        , ...
  'LineStyle'       , '-'           , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.7 .7 .7]    );

set(plotDataCS                            , ...
  'LineStyle'       , '--'          , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.7 .1 .1]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.7 .1 .1]    , ...
  'MarkerFaceColor' , [.7 .1 .1]    );

set(plotDataPanel                            , ...
  'LineStyle'       , ':'          , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.7 .1 .1]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.7 .1 .1]    , ...
  'MarkerFaceColor' , [.7 .1 .1]   );


hXLabel = xlabel('Age' ,'Interpreter','latex');
title_lbl = title('Survival rate ($S(a)$)','Interpreter','latex');

set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 12          );

set( [hXLabel, title_lbl]                      , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 18          );

lbl_model = text(7,0.28,'Model','Interpreter','latex');
lbl_dataPanel = text(4.1,0.71,'Data (Panel)','Interpreter','latex');
lbl_dataCS = text(1.7,0.93,'Data (Cross-section)','Interpreter','latex');

set([lbl_model,lbl_dataPanel,lbl_dataCS],'FontName','Helvetica','FontSize', 16);
set([lbl_dataPanel,lbl_dataCS],'Color',[.7 .1 .1]);
set(lbl_model,'Color',[.3 .3 .3]);


set(gcf, 'PaperPositionMode', 'auto');
if Print == 1
print(f3,'-depsc','Figures/Figure4_SurvivalFit.eps');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SALES DISTRIBUTION

sales_cdf = readtable("Data_SalesDistribution.csv");

FirmShare_data = sales_cdf(:,2); FirmShare_data = FirmShare_data{:,:};
SalesShare_data = sales_cdf(:,1); SalesShare_data = SalesShare_data{:,:};

N = 100;
n_vec = [1:N]';

theta = x/tau;
SalesShare_cdf =  1 - theta.^n_vec;


F  = (1 - theta)/theta * log(1/(1-theta)); %Mass of firms
FractionofFirms = 1./n_vec .* (1 - theta)/theta .* theta.^n_vec ./ F; 
FirmShare_cdf = cumsum(FractionofFirms);

ind = find(FirmShare_cdf > 0.99999,1);
FirmShare_cdf = FirmShare_cdf(1:ind);
SalesShare_cdf = SalesShare_cdf(1:ind);


x_val  = linspace(min(FirmShare_cdf),max(FirmShare_cdf),8);
SalesConcentration_theory = interp1(FirmShare_cdf,SalesShare_cdf,x_val);
SalesConcentration_data = interp1(FirmShare_data,SalesShare_data,x_val);


f4 = figure;

plot_theory = plot(x_val,SalesConcentration_theory);
hold on
plot_data = plot(x_val,SalesConcentration_data);
hold off
grid on


set(plot_data                            , ...
  'LineStyle'       , '--'          , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.7 .1 .1]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.7 .1 .1]    , ...
  'MarkerFaceColor' , [.7 .1 .1]    );

set(plot_theory                        , ...
  'LineStyle'       , '-'           , ...
  'Marker'          , 'o'           , ...
  'Color'           , [.3 .3 .3]    , ...
  'LineWidth'       , 2.5             , ...
  'MarkerSize'      , 6             , ...
  'MarkerEdgeColor' , [.2 .2 .2]    , ...
  'MarkerFaceColor' , [.7 .7 .7]    );


hXLabel = xlabel('Fraction of firms','Interpreter','latex');
title_lbl = title('Fraction of output','Interpreter','latex');

set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 12          );

set( [hXLabel, title_lbl]                      , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 18          );


lbl_model = text(0.51,0.25,'Model','Interpreter','latex');
lbl_data = text(0.8,0.15,'Data','Interpreter','latex');

set([lbl_model,lbl_data],'FontName','Helvetica','FontSize', 16);
set(lbl_data,'Color',[.7 .1 .1]);
set(lbl_model,'Color',[.3 .3 .3]);



set(gcf, 'PaperPositionMode', 'auto');
if Print == 1
print(f4,'-depsc','Figures/Figure4_SalesConcentrationFit.eps');
end


end