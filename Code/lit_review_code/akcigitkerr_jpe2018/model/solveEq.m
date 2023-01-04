function [A,Bn,eq,flag] = solveEq(p)
% Solves equilibrium

global alg
% fprintf(1,'-------------------------\n');
% fprintf(1,'-- Solving Equilibrium --\n');
% fprintf(1,'-------------------------\n');

%{ 
Here is the algorithm:
    (i)   Guess tau
    (ii)  Given tau, guess A to find the fix point for A such that r and z are compatible
    (iii) Given A, solve Bn
    (iv)  By using E{V(1,qj + exp increment} = nu*qbar, update tau

%}   

tauUp    = 1.0;
tauDown  = 0.0;
crittau  = 1.0;
iterAll  = 0.0;

while crittau>alg.crit2 && iterAll<=alg.maxIter2
    
    tau    = (tauUp + tauDown)/2;
    
    zUp      = 1.0;
    zDown    = 0.0;
    critz    = 1.0;
    iterz    = 0.0;
    
    while critz>alg.crit2 && iterz<=alg.maxIter2
        
        z       = (zUp + zDown)/2;
        
        growth  = tau*p.sbar + z*p.lambda;
        intRate = p.rho + growth;
        [A,~]   = solveA(p,tau,intRate);
        
        zNew    = ((A*p.lambda)/(p.chiHat*p.psiHat))^(1/(p.psiHat - 1));
        
        if z > zNew
            zUp = z;
        else
            zDown = z;
        end;
        critz = abs(z - zNew);
        iterz = iterz+ 1;
        
    end
    
    if iterz>alg.maxIter2
         warning('Max iter reached in z loop')
    end
    
    % Now calculate Bn 
    p.Omega  = A*(1 + p.theta*p.eta + (1 - p.theta)*p.eS);
    [Bn,xn,flagBn,flagsigmaTilde_1] = solveBn(p,tau);
   
    if flagBn ==1;
        alg.lastBn = Bn;
        EV_1 = A*(1 + p.sbar) + Bn(1);
        if p.nu > EV_1; 
          tauUp = tau; 
        else
          tauDown = tau;
        end
        crittau = abs(p.nu - EV_1);
        
    elseif flagsigmaTilde_1 == 0;  %
        tauDown = tau;
    end
   
   iterAll = iterAll + 1;
end


if iterAll<alg.maxIter2  && ~isnan(crittau)
        flag = 1;
        % disp('Equilibrium solved!')
        eq.xn     = xn;
        eq.z      = z;
        eq.tau    = tau;
        eq.growth = growth;

        eq = firmSizeDist(eq,p);
else
    if isnan(crittau)
        warning('crittau is Nan')
    else
        warning('Max iter reached in tau loop')
    end
    
    flag = 0;
    A    = 0;
    Bn   = 0;
    tau  = 0;
    xn   = 0;
    z    = 0;
    eq   = 0; 
end
    

    
end





