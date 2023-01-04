
function [eq] = firmSizeDist(eq,p)

%% Solves stationary distribution given the value of xn, tau and structural parameters

% f(n)  : share of firms with n product (NOT normalized to 1)
% mu(n) : share of product lines which belongs to a firm with n-product (sum is normalized to one for aggregate)
% M gives total mass of firms, given that total number of product line is normalized to one.

global alg

state  = [1:1:alg.NN+1]';
xn     = eq.xn;
tau    = eq.tau; 

%% Solving stationary distribution      

xEntUP   = 1.8;
xEntDOWN = 0;
critxEnt = 1;
iterAll  = 0;

while critxEnt>alg.crit1 && iterAll<=alg.maxIter
    
    xEnt    = (xEntUP + xEntDOWN)/2;
    fDist   = zeros(alg.NN,1);
    
    uplim   = 2;
    downlim = 0.0;
    crit1   = 1;
    iter    = 0;
    
    %% High type dist
    while crit1>alg.crit1 && iter<=alg.maxIter
        fDist(1) = (uplim+downlim)/2;
        M        =  xEnt/(fDist(1)*tau);
        fDist(2) = ((xn(1) + tau)*fDist(1) - xEnt/M)/(2*tau);
        fDist(2) = min(1,max(0,fDist(2)));
        fDist(3) = ((xn(2) + tau)*fDist(2)*2 - fDist(1)*xn(1))/(3*tau);
        fDist(3) = min(1,max(0,fDist(3)));
        
        for i = 3:alg.NN-1
            fDist(i+1) = (fDist(i)*i*(xn(i)+ tau) - fDist(i-1)*(i-1)*xn(i-1))/(tau*(i + 1));
            fDist(i+1) = min(1,max(0,fDist(i+1)));
        end
        
        summ = sum(fDist) ;
        
        if summ>1
            uplim = uplim - (uplim - fDist(1))/3;
        else
            downlim = downlim + (fDist(1) - downlim)/3;
        end;
        crit1 = abs(summ-1);
        iter = iter + 1;
    end
    
    % Calculate new value for creative destruction
    xEnt_new = max(tau - M*sum(state(1:alg.NN).*xn(1:alg.NN).*fDist),0);
    
    critxEnt  = abs(xEnt_new - xEnt);
    
    if xEnt_new > xEnt; 
         xEntDOWN = xEnt ; 
    else
         xEntUP   = xEnt;
    end
    iterAll = iterAll + 1;
end

if iterAll>alg.maxIter
        disp('Max iter reached')
        disp(critau);
end

eq.fDist = M*fDist;
eq.M     = M;


% Derived statistics

eq.mu = eq.fDist.*(1:1:alg.NN)';
eq.xEnt = xEnt;   % here, xEnt gives innovation rate*mass of entrants

%% Growth decomposition
eq.growthEntry = xEnt*p.sbar;
eq.growthExt   = sum(state(1:alg.NN).*xn(1:alg.NN).*eq.fDist)*p.sbar;
eq.growthInt   = eq.z*p.lambda;

%% R&D Decomposition
eq.RDCostEntry = xEnt*p.nu;
eq.RDCostExt   = sum((state(1:alg.NN).^((1 - p.sigma)*p.psiTilde)).*(xn(1:alg.NN).^p.psiTilde).*eq.fDist)*p.chiTilde;
eq.RDCostInt   = (eq.z^p.psiHat)*p.chiHat;
eq.RDTotal     = eq.RDCostEntry + eq.RDCostExt + eq.RDCostInt;
eq.xEnt        = xEnt;

end
% end of function






       
       



