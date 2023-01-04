function [Bn,xn,flagBn,flagsigmaTilde_1] = solveBn(p,tau)
%% Calculates Bn by uniformization

global alg


%Derived parameter
phi = (p.chiTilde^(1/p.psiTilde))*p.psiTilde; 
nEnd    = alg.NN;          % maximum number of products


% Preallocation
BnOld  = zeros(nEnd+1,1);
BnNew  = zeros(nEnd+1,1);
xn     = zeros(nEnd,1);
Ptilde = zeros(nEnd,3);

state  = (1:nEnd)';

% First solve Bn for sigmaTilde = 1 without using the loop
funcsigmaTilde_1 = @(x)(((p.rho+tau)*x/(p.psiTilde-1))^((p.psiTilde-1)/p.psiTilde))*phi - p.Omega - x;
init =  0.01;

[increment,fval,flag]   = fsolve(funcsigmaTilde_1,init,alg.options);   % solves for increment = Bn-Bn-1 
%[a,b,c]   = fzero(funcsigmaTilde_1,a);
flagsigmaTilde_1 = 1;
flagBn           = 1;
if flag ~= 1 || abs(fval)>1.0e-5
    flagsigmaTilde_1 = 0;
    flagBn = 0;
end
BnsigmaTilde_1  = increment*(1:nEnd+1)';       % value function when sigmaTilde = 1.
                   
B0      = 0;                                    % value at n = 0 
%BnOld   = alg.lastBn;
BnOld    = BnsigmaTilde_1;                      % we initialize the value function at its value for sigmaTilde=1.
%BnOld     = 0*ones(nEnd+1,1);                  % we can also initialize it at any other value.

% Now the value function iteration
crit = 1;
iter = 0;

if p.sigma==1/p.psiTilde
    Bn = BnsigmaTilde_1;
    xn = max(zeros(nEnd,1),((p.Omega + Bn(2:end,1) - Bn(1:end-1,1))./(p.psiTilde*(state.^(p.sigmaTilde - 1))*p.chiTilde)).^(1/(p.psiTilde-1)));
else
%% Uniformization starts here
xmax  = ((p.Omega + increment)/(p.psiTilde*p.chiTilde))^(1/(p.psiTilde-1)) + .05; 
nuBar = nEnd*(tau + xmax);          % maximal transition rate  
delta = nuBar/(p.rho+nuBar);        % new discount factor

while crit > alg.crit1 && iter<=alg.maxIter  
    
    xn        = max(zeros(nEnd,1),((p.Omega + BnOld(2:end,1) - BnOld(1:end-1,1))./(p.psiTilde*(state.^(p.sigmaTilde - 1))*p.chiTilde)).^(1/(p.psiTilde-1))); % xn is calculated given BnOld
    returnFun = (state.*xn*p.Omega - p.chiTilde*(state.^p.sigmaTilde).*(xn.^p.psiTilde))/(p.rho + nuBar);                                             % return function is calculated given BnOld  
    
    % Calculate 'expected Bnprime'
    Ptilde(:,1:2) = [state.*xn state*tau]/nuBar;
    Ptilde(:,3)   =  1 - sum(Ptilde(:,1:2),2);    % new transition probabilities. The format is [P(n,n+1)  P(n,n-1) P(n,n)]
    if sum(sum(Ptilde<0))
        Ptilde = max(0,Ptilde);
        alg.countZeros = sum(sum((Ptilde==0)));
    else 
        alg.countZeros = sum(sum((Ptilde==0)));
    end
    
    BnMat = [BnOld(2:end) [B0;BnOld(1:end-2)] BnOld(1:end-1)];
    BnExp = sum(Ptilde.*BnMat,2);
    
    % Update value function
    BnNew(1:nEnd) = returnFun(1:nEnd) + delta*BnExp(1:end);
    BnNew(end)    = 2*BnNew(end-1) - BnNew(end-2);              
    
    crit  = sum(abs(BnNew - BnOld));
    iter = iter + 1;
    BnOld = BnNew;
 
end



if iter < alg.maxIter && ~isnan(crit)
     Bn   = BnNew;
     flagBn = 1;
     
else
    if isnan(crit)
        flagBn = 0;
    else
        flagBn  = 2;
    end
    Bn    = Inf;
    xn    = Inf;
end

end

end












