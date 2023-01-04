function [transg, LI,vLovF,mu,g1,vLovF1,vL1,vL2,vF1,vF2] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec)
n = (length(pi)-1)/2;
nperiods = T;
r1 = r  ; r2 = r1-dr;
xinit = zeros(1,2*n); xinit(n+1) = 1;
[xvec1, muvec1, ~, g1] = gen_compute_eqm(lamb,pi,1,kap,r1,xinit); v = zeros(1,2*n); for j=1:2*n; v(j+1) = v(j)+xvec1(j); end; vL1 = v(n+1:2*n+1); vF1 = v(n+1:-1:1);
[xvec2, muvec2, ~, ~] = gen_compute_eqm(lamb,pi,1,kap,r2,xinit); v = zeros(1,2*n); for j=1:2*n; v(j+1) = v(j)+xvec2(j); end; vL2 = v(n+1:2*n+1); vF2 = v(n+1:-1:1);
vLovF1 = vL1*muvec1' / (vF1*muvec1');

etaL2 = xvec2(n+2:end); etaF2 = xvec2(n:-1:1); eta02 = xvec2(n+1);
transmu = zeros(length(muvec1),nperiods);
[transg,mu,vLovF]= deal(zeros(1,nperiods));

mu(1) = muvec1*(0:n)';
for t=1:nperiods
    if t>1
        mudown12 = transmu(2:end,t-1) .* (etaF2+kap)' * eps; mudown12 = [mudown12;0]; 
        muup12 = transmu(1:end-1,t-1) .* [2*eta02,etaL2]' * eps; muup12 = [0;muup12];
        mustay12 = transmu(1:end,t-1) .* (1-[2*eta02,etaF2+kap+[etaL2,0]]*eps)';
        transmu(:,t)= mustay12 + mudown12 + muup12;
    else
        transmu(:,t) = muvec1;
    end
    vLovF(t) = vL2*transmu(:,t) / (vF2*transmu(:,t));
    transg(t) = gen_compute_g(transmu(:,t)',xvec2,lamb,kap);
    LI(t) = pishrvec*transmu(:,t);
    mu(t) = transmu(:,t)'*(0:n)';
end


end