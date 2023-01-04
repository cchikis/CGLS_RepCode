function [p] = putParams(params)

global alg

% First get deep parameters
if isstruct(params)
   p = params;    % this part used for estimation
   if length(fieldnames(params))~=alg.nParam
        error('Number of arguments does not match')
   end
else 
   p = struct();
   if length(params)~=alg.nParam
       error('Number of arguments does not match')
   end
   [~,names]      = readParam(alg.paramFile);
   for i = 1:length(params)
       p.((char(names(i))))=params(i);
   end
end


% Derived parameters
p.sigmaTilde = (1 - p.sigma)*p.psiTilde;
p.eS         = p.alpha*p.theta*p.eta/(1 - (1 - p.theta)*p.alpha);   % expected step size for follow up improvement
p.sbar       = p.theta*p.eta/(1 - (1 - p.theta)*p.alpha);           % expected step size exploration R&D
p.LProdFin   = p.LProd/(((1 - p.beta)^2)/p.beta + 1);
p.LProdInt   = max(p.LProd - p.LProdFin,0) ;

p.betaTilde = (p.beta^p.beta)*((1 - p.beta)^(1 - 2*p.beta))*(p.zeta^(1 - p.beta));
p.pie       = p.LProdFin*(1 - p.beta)*p.betaTilde;
p.gamma     = p.gammaeta/p.eta;


end
    
    



