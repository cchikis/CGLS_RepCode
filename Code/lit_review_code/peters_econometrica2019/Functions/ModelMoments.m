% This function calculates the moments for given parameters and equilibrium
% policies (x,I,tau)


function m = ModelMoments(EqOutcomes, param)



% ENTRY RATE
m.EntryRate = EqOutcomes.x / (log(EqOutcomes.tau/(EqOutcomes.tau-EqOutcomes.x)));

% GROWHT RATE
m.AggGrowthRate = log(param.lambda) * (EqOutcomes.I + EqOutcomes.tau);

% LIFE CYCLE MARKUPS
lnmu_OLD = Avglnmu_Age(7.5, EqOutcomes.I, EqOutcomes.tau, EqOutcomes.x, param.lambda);
lnmu_Young = Avglnmu_Age(0.5, EqOutcomes.I, EqOutcomes.tau, EqOutcomes.x, param.lambda);

m.LCMU = lnmu_OLD-lnmu_Young;

% LIFE CYCLE EMPLOYMENT
lnlOLD = lnNbyAge(7.5, EqOutcomes.x, EqOutcomes.tau) - Avglnmu_Age(7.5, EqOutcomes.I, EqOutcomes.tau, EqOutcomes.x, param.lambda);
lnlYoung = lnNbyAge(0.5, EqOutcomes.x, EqOutcomes.tau) - Avglnmu_Age(0.5, EqOutcomes.I, EqOutcomes.tau, EqOutcomes.x, param.lambda);

m.LCEmpl = lnlOLD - lnlYoung;


end