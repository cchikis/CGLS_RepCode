function [nProd,age,exit,quality,step,numRadical,numInternal,numInnovation,qual_Step,numInnoq1,numInnoq2,numInnoq3,numInnoq4] = firmSim(p,eq,seed,algS)
%-------------------------%
% simulates firm dynamics
%-------------------------%

rng(seed,'combRecursive');

nSubPer = 1/algS.delt;
nBurn   = algS.nBurn/algS.delt;
nTKeep  = algS.nReal/algS.delt;                      % number of simulation period kept for statistics
nTSim   = nBurn + nTKeep;                            % number of simulation period (total)

                  
cdfPhi       = cumsum(p.theta*(1 - p.theta).^(0:400)');   % Analytical step size dist 
critStep     = find(cdfPhi>.1,1,'first')  - 1;            % critical step size lower than which counted as radical innovation
critStepq(1) = find(cdfPhi>.25,1,'first') - 1;
critStepq(2) = find(cdfPhi>.5,1,'first')  - 1;
critStepq(3) = find(cdfPhi>.75,1,'first') - 1;
critStepq(4) = find(cdfPhi>=1,1,'first')  - 1;


% determine the quality location of internal innovation 
qualExternal = p.eta.*p.alpha.^(0:100);
index1       = find(qualExternal<p.lambda,1,'first') - 1;
index2       = find(critStepq>index1,1,'first');  %  this gives which quartile the internal innnovation belongs.


% Preallocation for reporting
nProd    = ones (algS.nFirm,algS.nReal + 1);
age      = zeros(algS.nFirm,algS.nReal + 1);
exit     = zeros(algS.nFirm,algS.nReal + 1);
quality  = zeros(algS.nFirm,algS.NNExt,algS.nReal + 1);
step     = zeros(algS.nFirm,algS.NNExt,algS.nReal + 1);
numRadical    = zeros(algS.nFirm,algS.nReal + 1);   % number of major innovation in a year
numInternal   = zeros(algS.nFirm,algS.nReal + 1);   % number of internal innovation in a year
numInnovation = zeros(algS.nFirm,algS.nReal + 1);   % number of all innovations in a year

numInnoq1    = zeros(algS.nFirm,algS.nReal + 1);   
numInnoq2    = zeros(algS.nFirm,algS.nReal + 1);   
numInnoq3    = zeros(algS.nFirm,algS.nReal + 1);   
numInnoq4    = zeros(algS.nFirm,algS.nReal + 1);   

%nProdD    = ones(algS.nFirm,1);
nProdD    = zeros(algS.nFirm,1);
in = 1;
for s = 1:algS.NN
    nProdD(in:end,:) = s;
    in = in + ceil(eq.fDist(s)/eq.M*algS.nFirm);
end

qualityD  = zeros(algS.nFirm,algS.NNExt);
qualityD(:,1)  = 1;                             % All firms start with relative quality equal to 1 
stepD     = -1*ones(algS.nFirm,algS.NNExt);     % -1 means the quality entry is unoccupied.
stepD(:,1)= 0;    

for s = 1:algS.nFirm
    qualityD(s,1:nProdD(s)) = 1; 
    stepD(s,1:nProdD(s))= 0;   
end

 % Initial Quality Dist
qual_Step = [qualityD(:) stepD(:)];
qual_Step = qual_Step((qual_Step(:,2)>-1),:);

ageD      = zeros(algS.nFirm,1);
exitD     = zeros(algS.nFirm,1);

                              % start from zero step;  

firmIndex = (1:algS.nFirm)';
numRadicalD    = zeros(algS.nFirm,1);   % number of major innovation in a year
numInternalD   = zeros(algS.nFirm,1);   % number of internal innovation in a year
numInnovationD = zeros(algS.nFirm,1);   % number of all innovations in a year

numInnoqD    = zeros(algS.nFirm,4);


xDelt     = [eq.xn;eq.xn(end)*ones(algS.NNExt - algS.NN,1)]*algS.delt;   % just extending xn a little bit with constant innovation intensity 
zDelt     = eq.z*algS.delt;
tauDelt   = eq.tau*algS.delt;

subPer    = 1; 
period    = 1;    

for j = 1:nTSim
    %% State determination
    draw            = rand(algS.nFirm,1); 
    stateUpMajor    = p.theta*xDelt(nProdD).*nProdD;   
    stateUpFollow   = stateUpMajor + (1 - p.theta)*xDelt(nProdD).*nProdD;   
    stateDown       = stateUpFollow + tauDelt.*nProdD;                  % creative destruction rate is same for all firms, products.
    stateInternal   = stateDown + zDelt.*nProdD;
    if max(stateUpFollow)>=1
        warning('Probability is greater than 1')
    end
    upMajor         = (draw < stateUpMajor);            
    upFollow        = logical((1 - upMajor).*(draw < stateUpFollow));
    down            = logical((1 - upMajor - upFollow).*(draw < stateDown));
    internal        = logical((1 - upMajor - upFollow - down).*(draw < stateInternal));
    
    nProdD          = min(nProdD + upMajor + upFollow - down,algS.NNExt - 1);                         % update number of product line    
    ageD            = ageD + 1;                                   % update age 
%     [xx,yy]=max(nProdD)
 
    % Quality matrix modification due to down
    drawDown  = rand(sum(down),1);       % we randomly drop the quality
    nProdDOld = nProdD(down) + 1; 
    dropQuals = ceil(drawDown.*nProdDOld);
    firmIndDown = firmIndex(down);

    for s = 1:sum(down)
        qualityD(firmIndDown(s),dropQuals(s):end-1) = qualityD(firmIndDown(s),dropQuals(s)+1:end);  % shifting such that one randomly quality drops (this is kind of problematic if all quality matrix fills with product lines)
        stepD(firmIndDown(s),dropQuals(s):end-1) = stepD(firmIndDown(s),dropQuals(s)+1:end);  % shifting such that one randomly quality drops (this is kind of problematic if all quality matrix fills with product lines)
    end
    
    % Check if exit
    ifExit           = (nProdD<1);
    nProdD(ifExit)   = 1;                                       % enters again with one product
    firmIndExit      = firmIndex(ifExit);
    drawMajorOrFollow= rand(sum(ifExit),1);
    majorOrFollow    = (drawMajorOrFollow<p.theta);
    
    newQual_Step     = datasample(qual_Step,sum(ifExit));
    qualityD(firmIndExit,1) = newQual_Step(:,1) + p.eta*(majorOrFollow + (1 - majorOrFollow).*p.alpha.^(newQual_Step(:,2) + 1));
    stepD(firmIndExit,1)    = (1 - majorOrFollow).*(newQual_Step(:,2) + 1);     % update step                        
    ageD(ifExit)     = 0;                                       % reset age
    
   
    %% Update qualities
    
    % Quality matrix for upMajor
    newQual_Step       = datasample(qual_Step,sum(upMajor));

    IndMajor           = (nProdD(upMajor) - 1)*algS.nFirm + firmIndex(upMajor);
    qualityD(IndMajor) = newQual_Step(:,1) + p.eta;
    stepD(IndMajor)    = 0;    % since it is a major innovation
    
    
    % Quality matrix for upFollow
    newQual_Step         = datasample(qual_Step,sum(upFollow));
    IndFollow            = (nProdD(upFollow) - 1)*algS.nFirm + firmIndex(upFollow);
    qualityD(IndFollow)  = newQual_Step(:,1) + p.eta*p.alpha.^(newQual_Step(:,2) + 1);
    stepD(IndFollow)     = newQual_Step(:,2) + 1;
    
    % Quality matrix for internal innovation
    drawInternal   = rand(sum(internal),1); 
    nProdDInternal = nProdD(internal); 
    internalQuals  = ceil(drawInternal.*nProdDInternal);
    
    IndInternal            = (internalQuals - 1)*algS.nFirm + firmIndex(internal);
    qualityD(IndInternal)  = qualityD(IndInternal)*(1 + p.lambda);      
    % no change in step
    
    
    % Finalize new quality matrix
    qualityD          = qualityD /(1+eq.growth*algS.delt);         %adjust for the fact that qbar is growing with g
    
    % New quality distributioon
    qual_Step = [qualityD(:) stepD(:)];
    qual_Step = qual_Step((qual_Step(:,2)>-1),:);
    %qual_Step = sortrows(qual_Step,2);
%     [nelements,centers] = hist(qual_Step(:,1));
%     refreshdata(h,'caller')
%     drawnow
    
    
    if j == nBurn   % this is initial distribution
        nProd(:,period)     = nProdD; 
        age(:,period)       = ageD;
        quality(:,:,period) = qualityD;
        step(:,:,period)    = stepD;
        period              = period + 1;
    end
    
    if j>nBurn
        exitD = exitD + ifExit;    % track whether the firms exits within period after burning period.
        numRadicalD(upMajor)  = numRadicalD(upMajor) + 1;
        numRadicalD(upFollow) = numRadicalD(upFollow) + (stepD(IndFollow)<=critStep);  
        numInternalD          = numInternalD + internal;
        numInnovationD        = numInnovationD + upMajor + internal + upFollow;
        
        numInnoqD(upMajor ,1)  = numInnoqD(upMajor,1)  + 1;
        numInnoqD(upFollow,1)  = numInnoqD(upFollow,1) + (stepD(IndFollow)<=critStepq(1));
        numInnoqD(upFollow,2)  = numInnoqD(upFollow,2) + (stepD(IndFollow)<=critStepq(2)).*(stepD(IndFollow)>critStepq(1));
        numInnoqD(upFollow,3)  = numInnoqD(upFollow,3) + (stepD(IndFollow)<=critStepq(3)).*(stepD(IndFollow)>critStepq(2));
        numInnoqD(upFollow,4)  = numInnoqD(upFollow,4) + (stepD(IndFollow)<=critStepq(4)).*(stepD(IndFollow)>critStepq(3));
        
        numInnoqD(:,index2) = numInnoqD(:,index2) + (internal);  % putting internal innovation to the appropriate quartile

        if subPer == nSubPer
            
            % Reporting the end-year values
            nProd(:,period)         = nProdD; 
            age(:,period)           = ageD;
            quality(:,:,period)     = qualityD;
            step(:,:,period)        = stepD;
            exit(:,period)          = (exitD>0);
            numRadical(:,period)    = numRadicalD;
            numInternal(:,period)   = numInternalD;
            numInnovation(:,period) = numInnovationD;
            
            numInnoq1(:,period)    = numInnoqD(:,1);
            numInnoq2(:,period)    = numInnoqD(:,2);
            numInnoq3(:,period)    = numInnoqD(:,3);
            numInnoq4(:,period)    = numInnoqD(:,4);
            
            % reset
            numRadicalD    = zeros(algS.nFirm,1);   
            numInternalD   = zeros(algS.nFirm,1);   
            numInnovationD = zeros(algS.nFirm,1);  
            
            numInnoqD     = zeros(algS.nFirm,4);
            
            period = period + 1;
            subPer = 1;
        else
            subPer = subPer + 1;
        end
    end
%     
       
end

end


