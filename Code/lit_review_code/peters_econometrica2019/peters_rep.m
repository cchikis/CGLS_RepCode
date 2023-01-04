%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file calibrates the model in 
% "Heterogenous Markups,Growth and Endogenous Misallocation"
% and genertes all tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = peters_rep()
    %clc
    clear
    close all
    
    addpath('Functions/')
    addpath('DataInputs/')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Moments for calibration
    
    Data = readtable("Data_IndonesiaPanel.csv");
    Data_lnmu = Data(:,2); Data_lnmu = Data_lnmu{:,:};
    Data_lnl = Data(:,3); Data_lnl = Data_lnl{:,:};
    
    Moments.EntryRate = 0.104;          % Entry rate
    Moments.AggGrowthRate = 0.03;       % Rate of aggregate productivity growth
    Moments.LCMU = Data_lnmu(8);        % Markup lifecycle (8 years versus entrants) - for model: 7.5 versus 0.5
    Moments.LCEmpl = Data_lnl(8);       % Employment lifecyel (8 years versus entrants) - for model: 7.5 versus 0.5
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Structural parameters
    
    param.zeta = 2;                     % Convexity of innovation costs
    param.rho = 0.05;                   % Discount rate
    param.L = 1;                        % Labor force
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calibrate model (Table 3)
    
    % Step 1: Calibrate tau and x from entry rate and employment lifecycle
    [y_sol, fvalxtau] = CalibrateTauX(Moments);
    
    EqObjects.x = y_sol(1);
    EqObjects.tau = y_sol(2);
    EqObjects.z = EqObjects.tau - EqObjects.x;
    
    % Step 2: Calibrate I and lambda from markup licecycle and growth rate
    
    [I_sol, fvalI] = CalibrateI(Moments,EqObjects.x,EqObjects.tau);
    EqObjects.I = I_sol;
    param.lambda = exp(Moments.AggGrowthRate/(EqObjects.I + EqObjects.tau));
    
    % Step 3: Solve for accompanying structural parameters
    param = RecoverCosts(EqObjects,param);
    
    % Step 4: Check that model returns targeted moments
    EqOutcomes = FindEquilibrium(param);
    m = ModelMoments(EqOutcomes, param);
    
    rho_vec = linspace(5e-4, 0.065, 25);
    g_vec = zeros(1, size(rho_vec, 2));
    for (ii = 1:size(rho_vec, 2))
        param.rho = rho_vec(ii);
        EqOutcomes = FindEquilibrium(param);
        m = ModelMoments(EqOutcomes, param);
        
        g_vec(ii) = m.AggGrowthRate;
    end
    save_this = [rho_vec; g_vec];
    save('../../../Output/Store_Data/peters_res.mat', 'save_this')
    end
    
    % % Save results reported in Table 3
    % writetable(struct2table(structfun(@(x)round(x,3),param,'UniformOutput',false)),'Tables/Table3_Parameters.txt')
    % writetable(struct2table(structfun(@(x)round(x,3),EqObjects,'UniformOutput',false)),'Tables/Table3_EquilibriumObjects.txt')
    % writetable(struct2table(structfun(@(x)round(x,3),m,'UniformOutput',false)),'Tables/Table3_MomentsModel.txt')
    % writetable(struct2table(structfun(@(x)round(x,3),Moments,'UniformOutput',false)),'Tables/Table3_MomentsData.txt')
    
    % display('Baseline model calibrated (Table 3)')
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% Sensitivity Matrix (Table 4)
    
    % delta = 0.05;
    % [ChangesEqOutcomes,ChangesMoments] = CalculateSensitivityMatrix(param, delta);
    
    % T4_1 = array2table(round(ChangesEqOutcomes*100,2), 'VariableNames', {'phi_I', 'phi_x', 'phi_z', 'lambda'});
    % writetable(T4_1,'Tables/Table4_ChangeOutcomes.txt');
    % T4_2 = array2table(round(ChangesMoments*100,2), 'VariableNames', {'phi_I', 'phi_x', 'phi_z', 'lambda'});
    % writetable(T4_2,'Tables/Table4_ChangeMoments.txt');
    
    % display('Sensitivity matrix calculated (Table 4)')
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% Markups and Misallocation (Table 5)
    
    % m = EquilibriumMisallocation(EqObjects,param);
    
    % T5 = array2table(round(m,3), 'VariableNames', {'theta', 'Emu', 'stddevmu', 'stddevmu_f', 'M', 'Lambda'});
    % writetable(T5,'Tables/Table5_Misallocation.txt');
    
    % display('Equilibrium misallocation calculated (Table 5)')
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% Counterfactual (Table 6)
    
    % [EI_Indonesia,EI_US] = Counterfactual(param);
    
    % writetable(struct2table(structfun(@(x)round(x,3),EI_Indonesia,'UniformOutput',false)),'Tables/Table6_OutcomesIndonesia.txt')
    % writetable(struct2table(structfun(@(x)round(x,3),EI_US,'UniformOutput',false)),'Tables/Table6_OutcomesUS.txt')
    
    % display('Counterfactual calculated (Table 6)')
    
    
    
    
    
    
    