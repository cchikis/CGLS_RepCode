function [g,gvec]=gen_compute_g(muvec,xvec,lambda,kap)
    n = length(muvec)-1; xtmp = [xvec,0];
    gvec = zeros(1,n+1);
    gvec(1) = 2*muvec(1)*xvec(n+1);
    gvec(2:n+1) = muvec(2:n+1).*xtmp(n+2:(2*n+1));
    gvec = gvec.*log(lambda);
    g = sum(gvec);
    if ~isempty(kap)
        g = muvec(2:n+1)*(xtmp(n:-1:1)+kap)'*log(lambda);
    end
end
% This function is identical to LMS's function of the same name.
