function [a,flag] = solveA(p,tau,r)
% Solves A

global alg

func = @(A)(r + tau)*A - (p.pie + (A^(p.psiHat/(p.psiHat-1)))*((p.lambda/p.psiHat)^(p.psiHat/(p.psiHat-1)))*(p.psiHat - 1)*(p.chiHat^(1/(1 - p.psiHat))));

init = 1;

[a,fval,flag]   = fsolve(func,init,alg.options);


if flag ~= 1 || abs(fval)>1.0e-8
    warning('Could not solve A')
end

end