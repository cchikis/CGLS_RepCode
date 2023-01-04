% This function calculates the sensitivity matrix given a common change in
% the structual parameters by "delta" (i.e. the new parameter is given by
% "(1+delta)*initial parameter". 


function [ChangesEqOutcomes,ChangesMoments] = CalculateSensitivityMatrix(param, delta)


param_orig = param;

%% Calculate calibrated model to save results in one place

% Equilibrium
tmp = FindEquilibrium(param_orig);
    
EquilibriumOutcomesCalibration = zeros(4,1);
EquilibriumOutcomesCalibration(1) = tmp.I;
EquilibriumOutcomesCalibration(2) = tmp.x;
EquilibriumOutcomesCalibration(3) = tmp.z;
EquilibriumOutcomesCalibration(4) = tmp.tau;

% Moments
m = ModelMoments(tmp, param_orig);
   
MomentsCalibration = zeros(4,1);
MomentsCalibration(1) = m.LCMU;
MomentsCalibration(2) = m.LCEmpl;
MomentsCalibration(3) = m.EntryRate;
MomentsCalibration(4) = m.AggGrowthRate;

% Structural Parameter
struct_param = ones(4,1);
struct_param(1) = param_orig.phiI;
struct_param(2) = param_orig.phix;
struct_param(3) = param_orig.phiz;
struct_param(4) = param_orig.lambda;


%% Calculate matrix of structual parameters with changes
param_variation = repmat(struct_param,1,4);

for j = 1:3
    param_variation(j,j) = struct_param(j)*(1 + delta);    
end

param_variation(4,4) = 1 + (struct_param(4)-1)*(1 + delta);

%% Calculate equilibrium for all variations

EnodgenousOutcomes = zeros(4,4);

for j = 1:4
    % Set structural parameters in changes 
    param.phiI = param_variation(1,j);
    param.phix = param_variation(2,j);
    param.phiz = param_variation(3,j);
    param.lambda = param_variation(4,j);
    
    % Calculate equilibrium
    EqOutcomes = FindEquilibrium(param);
    
    % Save results
    EnodgenousOutcomes(1,j) = EqOutcomes.I;
    EnodgenousOutcomes(2,j) = EqOutcomes.x;
    EnodgenousOutcomes(3,j) = EqOutcomes.z;
    EnodgenousOutcomes(4,j) = EqOutcomes.tau;
end

%% Given endogenous outcomes calculate moments
Moments = zeros(4,4);
for j = 1:4
    EqOutcomesSensitivity.I = EnodgenousOutcomes(1,j);
    EqOutcomesSensitivity.x = EnodgenousOutcomes(2,j);
    EqOutcomesSensitivity.z = EnodgenousOutcomes(3,j);
    EqOutcomesSensitivity.tau = EnodgenousOutcomes(4,j);
    
    param.lambda = param_variation(4,j); %Lambda is the only parameter we need
    
    m = ModelMoments(EqOutcomesSensitivity, param);
   
    Moments(1,j) = m.LCMU;
    Moments(2,j) = m.LCEmpl;
    Moments(3,j) = m.EntryRate;
    Moments(4,j) = m.AggGrowthRate;
    
end


%% Export results 

% Levels
MomentSensitivity = [MomentsCalibration,Moments];
EqOutcomesSensitivity = [EquilibriumOutcomesCalibration,EnodgenousOutcomes];

% Changes 

ChangesMoments = MomentSensitivity(:,2:5)./repmat(MomentsCalibration,1,4) - 1;
ChangesEqOutcomes = EqOutcomesSensitivity(:,2:5)./repmat(EquilibriumOutcomesCalibration,1,4) - 1;

end