function [xvec, muvec, gvec, g, res,flag] = gen_compute_eqm_state_dep_kappa(lambda,pivec,param,kappa,r,xinit)
    
    n = (length(pivec)-1)/2;
    if isempty(xinit)
        xinit = zeros(1,2*n); xinit(n+1) = 1;
    end
    tol=1e-12;
    %opt = optimoptions('fsolve','Display','off','FunctionTolerance',tol,'OptimalityTolerance',tol,'MaxFunctionEvaluations',50000,'MaxIterations',1000);
    opt = optimoptions('fsolve','Display','off','FunctionTolerance',tol,'OptimalityTolerance',tol);
    [xvec,~,flag] = fsolve(@(xvec)gen_eqm_eqns_state_dep_kappa(xvec,pivec,param,kappa,r),xinit,opt);xvec=abs(xvec);
    res = sum(abs(gen_eqm_eqns_state_dep_kappa(xvec,pivec,param,kappa,r)));
    xvec=xvec/r;
    muvec = zeros(1,n+1);
    muvec(1) = 1; muvec(2) = 2*xvec(n+1) / (xvec(n)+kappa(1));
    xtmp = [xvec,0];
    for i=2:n
        muvec(i+1) = muvec(i)* xtmp(n+i) / (xtmp(n+1-i)+kappa(i));
    end
    muvec = muvec./sum(muvec);
    
    n = length(muvec)-1; xtmp = [xvec,0];
    gvec = zeros(1,n+1);
    gvec(1) = 2*muvec(1)*xvec(n+1);
    gvec(2:n+1) = muvec(2:n+1).*xtmp(n+2:(2*n+1));
    gvec = gvec.*log(lambda);
    g = sum(gvec);
end