% This function calibrates the innovation rate to the life-cycle of log
% mark-ups given an aggregate growth rate



function [I_sol, fval] = CalibrateI(Moments,x,tau)
Minobject = @(I) FindIObjective(I,Moments,x,tau);
[I_sol, fval] = fminbnd(Minobject,0,10);
end



function diff = FindIObjective(I,Moments,x,tau)

% Step 1: Get lambda to match growth rate
lambda = exp(Moments.AggGrowthRate/(I + tau));

% Step 2: Calculate log markups
lnmu_OLD = Avglnmu_Age(7.5, I, tau, x, lambda);
lnmu_Young = Avglnmu_Age(0.5, I, tau, x, lambda);

% Step 3: Calculate difference
diff = ((lnmu_OLD-lnmu_Young) - Moments.LCMU)^2;

end