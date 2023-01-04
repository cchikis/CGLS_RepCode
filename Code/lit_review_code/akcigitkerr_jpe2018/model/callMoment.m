function [score,m,qualityFirmA, nProdA, qual_StepA, stepA,aux ,aux2] = callMoment(params,eq)
%% Calculates moments
global alg;
% fprintf(1,'\n');
% fprintf(1,'---------------------\n');
% fprintf(1,'-- FIRM SIMULATION --\n');
% fprintf(1,'---------------------\n');
% 

pS = params;

% Preallocation
nProd           = cell(1,alg.nThread);
exit            = cell(1,alg.nThread);
qualityFirm     = cell(1,alg.nThread);
numRadical      = cell(1,alg.nThread);
numInternal     = cell(1,alg.nThread);
numInnovation   = cell(1,alg.nThread);
step            = cell(1,alg.nThread);
qual_Step       = cell(1,alg.nThread);
numInnoq1       = cell(1,alg.nThread);
numInnoq2       = cell(1,alg.nThread);
numInnoq3       = cell(1,alg.nThread);
numInnoq4       = cell(1,alg.nThread);

m     = struct();
mm    = struct();

% moment storage
mm.momentVec  = [];
mm.momentName = {};
mm.momentWgt  = [];
mPos          = 1;

seed = alg.seed;

algSim = cell(1,alg.nThread);
for ii = 1:alg.nThread
    algSim{ii} = alg; 
end


% parallelize firm simulation
if alg.parallel == 1
  if exist('gcp')
      pool = gcp('nocreate');
      if isempty(pool)
          parpool('local',alg.nThread);
      end
  else
      try matlabpool local alg.nThread, catch ME, end
  end
end

%% Firm Simulation
fprintf(1,'Simulation...')       
parfor (i = 1:alg.nThread,alg.parallel*alg.nCores)
   [nProd{i},~,exit{i},qualityFirm{i},step{i},numRadical{i},numInternal{i},numInnovation{i},qual_Step{i},numInnoq1{i},numInnoq2{i},numInnoq3{i},numInnoq4{i}] = firmSim(pS,eq,seed(i),algSim{i});
end
fprintf(1,' Done!\n')        
    

% Concatenate threads
nProdA          = vertcat(nProd{:});
qualityFirmA    = vertcat(qualityFirm{:});
exitA           = vertcat(exit{:});
numRadicalA     = vertcat(numRadical{:});
numInternalA    = vertcat(numInternal{:});
numInnovationA  = vertcat(numInnovation{:});
stepA           = vertcat(step{:});       
qual_StepA      = vertcat(qual_Step{:});
numInnoq1A      = vertcat(numInnoq1{:}); 
numInnoq2A      = vertcat(numInnoq2{:});
numInnoq3A      = vertcat(numInnoq3{:});
numInnoq4A      = vertcat(numInnoq4{:});



% Another selection of sample
%[numInnovation,numCreateDest] = firmSimforSelection(p,eq,nProdD,seed,algS)
per1 = 1;
per2 = 2;

firmSize     = sum(qualityFirmA,2);
firmSize     = reshape(firmSize,size(firmSize,1),size(firmSize,3));
firmSize(exitA == 1)   = 0;
nProdA(exitA == 1  )   = 0;


% extend xn a little bit in case!
eq.xn = [eq.xn;eq.xn(end)*ones(alg.NNExt - alg.NN,1)];

%% MOMENTS

% 1 - R&D expenditure by firm
rdExp = pS.chiTilde*(nProdA(:,per1).^pS.sigmaTilde).*(eq.xn(nProdA(:,per1)).^pS.psiTilde) + pS.chiHat*(eq.z^pS.psiHat)*firmSize(:,per1);
% Sales
sales = (((1 - pS.beta)*pS.zeta/pS.betaTilde)^((1 - pS.beta)/pS.beta))*pS.LProdFin*firmSize(:,per1);
% R&D intensity
m.rdIntense = sum(rdExp)/sum(sales);

% 2 - the fraction of patents that are internally developed
m.fractInternal = eq.z/(eq.z + eq.tau);

% 3 - growth rate
m.growth = eq.growth;

% 4 - Profitibility
m.profitabilty = sum(pS.pie*firmSize(:,per1))/sum(sales);

% 5 - entry rate
m.entryRate = eq.xEnt/eq.M;

% 6 - o	(Internal patent expected # of cites) / (external patent expected # of cites)
m.intCite_extCite = pS.lambda/pS.sbar;

% 7 - Regression growth
% Annual
growthFirmAn                    = (firmSize(:,per2) - firmSize(:,per1))./firmSize(:,per1);
growthFirmAn(growthFirmAn>10)   = alg.growthCap;   % Growth cap

X1                              = [ones(length(growthFirmAn),1) log(firmSize(:,per1))];
m.betaGrowth                    = regress(growthFirmAn,X1);

id = (1:1:length(growthFirmAn))';
aux2 = [growthFirmAn firmSize(:,per1) id];

numInnovationA(exitA==1) = 0;
numInternalA(exitA==1)   = 0;
numRadicalA(exitA==1)    = 0;
numInnoq1A(exitA==1)     = 0;
numInnoq2A(exitA==1)     = 0;
numInnoq3A(exitA==1)     = 0;
numInnoq4A(exitA==1)     = 0;

firmSize(logical(exitA(:,2:end))) = 0;      % if firm exits, we set employment equal to zero at the begining of the period. 
firmSize(:,end) = 0;                        % we will not use the last employment

sumInnovation                  = sum(numInnovationA,2);
sumnumInternal                 = sum(numInternalA,2);
sumnumRadical                  = sum(numRadicalA,2);
sumFirmSize                    = sum(firmSize,2);
sumInnoq1                      = sum(numInnoq1A,2);
sumInnoq2                      = sum(numInnoq2A,2);
sumInnoq3                      = sum(numInnoq3A,2);
sumInnoq4                      = sum(numInnoq4A,2);

aux = [sumInnovation,sumnumInternal,sumnumRadical,sumFirmSize,firmSize(:,1) id sumInnoq1 sumInnoq2 sumInnoq3 sumInnoq4];
aux(aux(:,4)==0,:) = [];    % get rid of zero employment firms

count = sum(logical(firmSize>0),2);
count(count==0)=[];

aux   = [aux count]; 

X2                              = [ones(length(aux),1) log(aux(:,4)./count)]; 
XX2                             = [ones(length(aux),1) log(aux(:,5))];  

pat_Emp                         = zscore(aux(:,1)./aux(:,4));
fracRad                         = aux(:,3)./aux(:,1);

m.betaPat_emp                   = regress(pat_Emp,X2);
m.betaFracRad                   = regress(fracRad,XX2);

%% ---------------------------------------------------------------------------------
addMoment(m.rdIntense,      'R&D intensity',alg.wgtrdIntense);                     % moment-1
addMoment(m.fractInternal,  'fraction of internal patents',alg.wgtfractInternal);  % moment-2
addMoment(m.growth,         'growth',alg.wgtgrowth)                                % moment-3
addMoment(m.profitabilty,   'profitabilty',alg.wgtprofitabilty)                    % moment-4
addMoment(m.entryRate,      'entry rate',alg.wgtentryRate)                         % moment-5
addMoment(m.intCite_extCite,'internal/external cite',alg.wgtintCite_extCite)       % moment-6
addMoment(m.betaGrowth(2),  'growth vs size (fact 1)' ,alg.wgtbetaGrowth);         % moment-7
addMoment(m.betaPat_emp(2), 'patent per emp vs size (fact 3)',alg.wgtbetaPat_emp); % moment-8
addMoment(m.betaFracRad(2), 'top innov. vs size (fact 2)',alg.wgtbetaFracRad);     % moment-9
%% ---------------------------------------------------------------------------------


%%%%%%%%%%%%%%%
%% SET SCORE %%
%%%%%%%%%%%%%%%
mm.nMoment = length(mm.momentVec);

% load data moments
mm.momentData = load(alg.momentDataFile);
m.momentData  = mm.momentData;
m.momentName  = mm.momentName;
m.momentWgt   = mm.momentWgt';
mm.mErr       = mm.momentData' - mm.momentVec;

score = (sum(mm.momentWgt'.*abs(mm.mErr')./(0.5*abs(mm.momentVec') + 0.5*abs(mm.momentData))))/mm.nMoment;

m.momentDataVSModel = [mm.momentData mm.momentVec'];

%------------------------------------%
function addMoment(val,str,wgt)
  mm.momentVec(mPos)  = val;
  mm.momentName{mPos} = str;
  mm.momentWgt(mPos)  = wgt;
  mPos                = mPos + 1;
end

end

