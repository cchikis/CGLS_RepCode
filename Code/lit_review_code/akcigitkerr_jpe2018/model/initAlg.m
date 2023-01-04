function initAlg(pid)
  global alg
  
  %% Add path
  
  path('Utility Functions',path);

  
  %% Parameters to use
  alg.ptag = 'Baseline';

  % paramet set override
  if (nargin == 1)
    alg.ptag = pid;
  end
  

  %% Parameters for solver
  alg.NN        = 30;                      % Max. number of product lines
  alg.NNExt     = 35;                      
  alg.maxIter   = 30000;    
  alg.maxIter2  = 1000;  
  alg.crit1     = 1e-7;
  alg.crit2     = 1e-5;
  alg.final     = 0;
  alg.options   = optimset(optimset('fsolve'), 'TolFun', 1.0e-12, 'TolX',1.0e-12, 'Display','off');
  alg.lastBn    = zeros(alg.NN+1,1);


  %% Parameters for simulation
  alg.delt      = 0.02;                  % for converting continuous time to discrete time: implies that a year has 1/delt subperiod
  alg.nFirm     = 2^12;                 
  alg.nProd     = 2^10;                  % Number of product lines for product simulation    
  alg.nBurn     = 500;                   % Number of years for burning
  alg.nReal     = 5;                     % Number of years for calculating statistics 
  alg.stepMax   = 50;                    % Max number of step for follow up innovation   
  alg.growthCap = 10;                    % growth cap
  
  alg.parallel  = 1;                     % 1 for parallelization, 0 otherwise 
  alg.nThread   = 2;                     % Number of threads for parallel computing
  alg.seed      = [4324,6565];          
  alg.nCores    = 2;
  
 


  %% Filenames
  alg.paramFile       = ['Params/params'  alg.ptag '.txt'];
  alg.momentDataFile  = ['Moments/moment' alg.ptag '.txt'];
  alg.estResults      = ['Output/result'  alg.ptag '.txt'];
  alg.summaryBestStep = ['Logs/summaryBestStep' alg.ptag '.txt'];

  %% Parameter Names
  alg.paramNames{1}  = 'sigma';
  alg.paramNames{2}  = 'psiTilde';
  alg.paramNames{3}  = 'chiTilde';
  alg.paramNames{4}  = 'psiHat';
  alg.paramNames{5}  = 'chiHat';
  alg.paramNames{6}  = 'rho';
  alg.paramNames{7}  = 'theta';
  alg.paramNames{8}  = 'eta';
  alg.paramNames{9}  = 'lambda';
  alg.paramNames{10} = 'alpha';
  alg.paramNames{11} = 'LProd';
  alg.paramNames{12} = 'beta';
  alg.paramNames{13} = 'zeta';
  alg.paramNames{14} = 'nu';         
  alg.paramNames{15} = 'gammaeta';
  
  alg.nParam         = length(alg.paramNames);
  
  %% Moment Weights
  alg.wgtrdIntense        = 1.0; 
  alg.wgtfractInternal    = 1.0;
  alg.wgtgrowth           = 3.0;
  alg.wgtprofitabilty     = 1.0; 
  alg.wgtentryRate        = 1.0;
  alg.wgtintCite_extCite  = 1.0;
  alg.wgtbetaGrowth       = 1.0;

  alg.wgtbetaPat_emp      = 0.0;
  alg.wgtbetaFracRad      = 0.0;

  if strcmp(pid,'ExtensionFact12')
    alg.wgtbetaFracRad      = 1.0;
  elseif strcmp(pid,'ExtensionFact123')
    alg.wgtbetaFracRad      = 1.0;
    alg.wgtbetaPat_emp      = 1.0;
  elseif strcmp(pid,'ExtensionGrowthCap')  
      alg.growthCap = 30;
  end

  
end

