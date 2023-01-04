% This function calibrates the entry and expansion effciencies for the US

function [paramUS, fval] = Calibration_US(NewMoments,param)

Minobject = @(y) Calibration_US_Objective(y, NewMoments, param);

y_init = [param.phix, param.phiz]; 
LB = [0, 0];                           % Lower bound on productivites

options = optimoptions('fmincon','Display','off');
[y_sol, fval] = fmincon(Minobject, y_init, [], [], [], [], LB, [], [], options);


paramUS = param;
paramUS.phix = y_sol(1);
paramUS.phiz = y_sol(2);


end


function diff = Calibration_US_Objective(y, NewMoments, param)

param.phix = y(1);      % Expansion productivity
param.phiz = y(2);      % Entry productivity

% New equilibrium 
EqOutcomes = FindEquilibrium(param);

% New Moments
NewEntryRate = EqOutcomes.x / (log(EqOutcomes.tau/(EqOutcomes.tau-EqOutcomes.x)));


lnL_old = lnNbyAge(12.5, EqOutcomes.x, EqOutcomes.tau) - Avglnmu_Age(12.5, EqOutcomes.I, EqOutcomes.tau, EqOutcomes.x, param.lambda);
lnL_young = lnNbyAge(2.5, EqOutcomes.x, EqOutcomes.tau) - Avglnmu_Age(2.5, EqOutcomes.I, EqOutcomes.tau, EqOutcomes.x, param.lambda);
NewLClnEmpl = lnL_old - lnL_young;


% Objective function
diff = (NewMoments.EntryRate - NewEntryRate)^(2) + (NewMoments.LClnEmpl - NewLClnEmpl)^(2);


end
