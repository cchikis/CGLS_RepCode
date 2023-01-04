%% Replication file for "Growth through Heterogeneous Innovations", 
%  by Ufuk Akcigit and William R. Kerr
%  May 2018
% clearvars -global
% clear all;
% close all;
% clc;
function [] = ak_rep()

    % disp('------------------');
    % disp(' BASELINE RESULTS ')
    % disp('------------------');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global alg;
    format short;
    initAlg('Baseline');
    pS = readParam(alg.paramFile);  % load parameters
    p  = putParams(pS);
    paramBaselineTable = table(cell2mat(struct2cell(pS)),'RowNames',alg.paramNames,'VariableNames',{'Parameters'}); % Table 3 and 5
    %fprintf('-------------------------\n')
    rho_vec = linspace(5e-4, 0.065, 25);
    g_vec = zeros(size(rho_vec));
    for (ii = 1:length(rho_vec))
        p.rho = rho_vec(ii);
        [~,~,eq ] = solveEq(p);         % solve model
        g_vec(ii) = eq.growth;
        %eq.growth
        %[~,m,qualityFirmA, nProdA] = callMoment(p,eq);   % firm simulation and moment calculations
    end
    save_this = [rho_vec; g_vec];
    save('../../../../Output/Store_Data/ak_res.mat', 'save_this')
    end
    
    % momentBaselineTable = table(m.momentDataVSModel(m.momentWgt>0,1),m.momentDataVSModel(m.momentWgt>0,2),...
    %                             'RowNames',m.momentName(m.momentWgt>0),'VariableNames',{'Data','Model'})    % Table 4
    % %fprintf('Estimated sigma = %1.3f, Implied sigma + psi = %1.3f\n',p.sigma, p.sigma+ 1/p.psiTilde);                          
    % %fprintf('-------------------------\n')
    % growthDecompTable   = table([eq.growthInt, eq.growthExt,eq.growthEntry]',...
    %                             'RowNames',{'Internal','External','New Entry'},'VariableNames',{'Growth_Decomposition'})  % Table 6
    % %fprintf('-------------------------\n')
    
    % save(['Mat files' filesep 'resultBaseline.mat']);   % results are saved 
    
    % varyingSigma();         % Varying sigma under baseline parametrization (for Figure 8 and 9)
    % figures;                % figures as appeared in the paper. They are saved under subfolder "Figures".
    % %fprintf('\n\n\n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % disp('-------------');
    % disp(' EXTENSTIONS ')
    % disp('-------------');
    
    % disp('Adding Fraction of Top Innovations as a Target')
    
    % clear all;
    % global alg;
    % initAlg('ExtensionFact12');
    % pS = readParam(alg.paramFile);  % load parameters
    % p  = putParams(pS);
    
    % [~,~,eq] = solveEq(p);          % solve model
    % [~,m]    = callMoment(p,eq);    % firm simulation and moment calculations
    
    % momentExtensionFact12Table = table(m.momentDataVSModel(m.momentWgt>0,1),m.momentDataVSModel(m.momentWgt>0,2),...
    %                             'RowNames',m.momentName(m.momentWgt>0),'VariableNames',{'Data','Model'}) % Table 10
    % %fprintf('Estimated sigma = %1.3f, Implied sigma + psi = %1.3f\n',p.sigma, p.sigma+ 1/p.psiTilde);                        
    % %fprintf('-------------------------\n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % disp('Adding Patent per Employment as a Target')
    
    % clear all;
    % global alg;
    % initAlg('ExtensionFact123');
    % pS = readParam(alg.paramFile);  % load parameters
    % p  = putParams(pS);
    
    % [~,~,eq] = solveEq(p);          % solve model
    % [~,m]    = callMoment(p,eq);    % firm simulation and moment calculations
    
    % momentExtensionFact123Table = table(m.momentDataVSModel(m.momentWgt>0,1),m.momentDataVSModel(m.momentWgt>0,2),...
    %                             'RowNames',m.momentName(m.momentWgt>0),'VariableNames',{'Data','Model'}) % Table 11
    % %fprintf('Estimated sigma = %1.3f, Implied sigma + psi = %1.3f\n',p.sigma, p.sigma+ 1/p.psiTilde);                          
    % %fprintf('-------------------------\n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % disp('Alternative Growth Cap')
    
    % clear all;
    % global alg;
    % initAlg('ExtensionGrowthCap');
    % pS = readParam(alg.paramFile);  % load parameters
    % p  = putParams(pS);
    
    % [~,~,eq] = solveEq(p);          % solve model
    % [~,m]    = callMoment(p,eq);    % firm simulation and moment calculations
    
    % momentExtensionGrowthCap30Table = table(m.momentDataVSModel(m.momentWgt>0,1),m.momentDataVSModel(m.momentWgt>0,2),...
    %                             'RowNames',m.momentName(m.momentWgt>0),'VariableNames',{'Data','Model'})  % Table 12
    % %fprintf('Estimated sigma = %1.3f, Implied sigma + psi = %1.3f\n',p.sigma, p.sigma+ 1/p.psiTilde);                          
    % %fprintf('-------------------------\n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % disp('Alternative R&D Elasticities: psi = 0.4')
    
    % clear all;
    % global alg;
    % initAlg('ExtensionPsi04');
    % pS = readParam(alg.paramFile);  % load parameters
    % p  = putParams(pS);
    
    % [~,~,eq] = solveEq(p);          % solve model
    % [~,m]    = callMoment(p,eq);    % firm simulation and moment calculations
    
    % momentExtensionPsi04Table = table(m.momentDataVSModel(m.momentWgt>0,1),m.momentDataVSModel(m.momentWgt>0,2),...
    %                             'RowNames',m.momentName(m.momentWgt>0),'VariableNames',{'Data','Model'})  % Table 13
    % %fprintf('Estimated sigma = %1.3f, Implied sigma + psi = %1.3f\n',p.sigma, p.sigma+ 1/p.psiTilde);                          
    % %fprintf('-------------------------\n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % disp('Alternative R&D Elasticities: psi = 0.6')
    
    % clear all;
    % global alg;
    % initAlg('ExtensionPsi06');
    % pS = readParam(alg.paramFile);  % load parameters
    % p  = putParams(pS);
    
    % [~,~,eq] = solveEq(p);          % solve model
    % [~,m]    = callMoment(p,eq);    % firm simulation and moment calculations
    
    % momentExtensionPsi06Table = table(m.momentDataVSModel(m.momentWgt>0,1),m.momentDataVSModel(m.momentWgt>0,2),...
    %                             'RowNames',m.momentName(m.momentWgt>0),'VariableNames',{'Data','Model'}) % Table 13
    % %fprintf('Estimated sigma = %1.3f, Implied sigma + psi = %1.3f\n',p.sigma, p.sigma+ 1/p.psiTilde);                          
    % %fprintf('-------------------------\n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % % end
    
    
    
    