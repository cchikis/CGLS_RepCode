function [xvec, muvec, gvec, g, res,flag] = gen_compute_eqm_correct(lambda,pivec,param,kappa,r,xinit,c)
% See a description of this function in the comment beginning on Line 30 below.    
    n = (length(pivec)-1)/2;
    if isempty(xinit)
        xinit = zeros(1,2*n); xinit(n+1) = 1;
    end
    tol=1e-12;
    %opt = optimoptions('fsolve','Display','off','FunctionTolerance',tol,'OptimalityTolerance',tol,'MaxFunctionEvaluations',50000,'MaxIterations',1000);
    opt = optimoptions('fsolve','Display','off','FunctionTolerance',tol,'OptimalityTolerance',tol);
    [xvec,~,flag] = fsolve(@(xvec)gen_eqm_eqns_correct(xvec,pivec,param,kappa,r,c),xinit,opt);xvec=abs(xvec);
    res = sum(abs(gen_eqm_eqns_correct(xvec,pivec,param,kappa,r,c)));
    %xvec=xvec/r;
    muvec = zeros(1,n+1);
    muvec(1) = 1; muvec(2) = 2*xvec(n+1) / (xvec(n)+kappa);
    xtmp = [xvec,0];
    for i=2:n
        muvec(i+1) = muvec(i)* xtmp(n+i) / (xtmp(n+1-i)+kappa);
    end
    muvec = muvec./sum(muvec);
    
    [g,gvec]=gen_compute_g(muvec,xvec,lambda,kappa);
    if 1==0
    gvec(1) = 2*muvec(1)*xvec(n+1);
    gvec(2:n+1) = muvec(2:n+1).*xtmp(n+2:(2*n+1));
    gvec = gvec.*log(lambda);
    g = sum(gvec);
    end
end

% Note that all results in our main text take the investment success rates in LMS's 
% paper and replication code as given.  
%
% This function (gen_compute_eqm_correct.m) is used to obtain the middle row 
% of Figure IA.5. 
%
% This function is identical to gen_compute_eqm.m in LMS's replication
% code, except that:
%
% (i) Line 12 is commented out.  Commenting out Line 12 does
% not change the equations solved or the interpretation of the function's
% output, because Line 3 of gen_eqm_eqns_correct is also commented out.
% See our Online Appendix C.1.2 (``Note on ``dividing by r''...) for additional discussion.  
%
% (ii) Line 10 calls gen_eqm_eqns_correct.m, not LMS's gen_eqm_eqns_correct.m
% gen_eqm_eqns_correct.m solves the model using the investment cost assumption
% provided in LMS's paper: (c^2 eta_s)^2.  See our Online Appendix C for details.
% LMS state (p. 213) that, in their quantitative model, the cost of achieving
% investment success rate eta is c(eta)=(c * \eta)^2.  This functional form 
% assumption is used here.  

