% A description of how this function differs from LMS's function of the
% same name is provided at the end of this function.  

function [transg, transg_comp, LI,vLovF,mu,g1,vLovF1,vL1,vL2,vF1,vF2] = sc18_transition(T,r,dr,eps,lamb, pi, kap,pishrvec, percent_adjust, lL, lF, l0,sig, PDA_leader)
n = (length(pi)-1)/2;
nperiods = T;
r1 = r  ; r2 = r1-dr;
xinit = zeros(1,2*n); xinit(n+1) = 1;
[xvec1, muvec1, ~, g1] = gen_compute_eqm(lamb,pi,1,kap,r1,xinit); v = zeros(1,2*n); for j=1:2*n; v(j+1) = v(j)+xvec1(j); end; vL1 = v(n+1:2*n+1); vF1 = v(n+1:-1:1);
if r-dr < 0.5
    xinit = xvec1;
end
[xvec2, muvec2, ~, ~] = gen_compute_eqm(lamb,pi,1,kap,r2,xinit); v = zeros(1,2*n); for j=1:2*n; v(j+1) = v(j)+xvec2(j); end; vL2 = v(n+1:2*n+1); vF2 = v(n+1:-1:1);
vLovF1 = vL1*muvec1' / (vF1*muvec1');

etaL2 = xvec2(n+2:end); etaF2 = xvec2(n:-1:1); eta02 = xvec2(n+1);
transmu = zeros(length(muvec1),nperiods);
[transg,mu,vLovF]= deal(zeros(1,nperiods));

mu(1) = muvec1*(0:n)';
for t=1:nperiods
    if t>1
        mudown12 = transmu(2:end,t-1) .* (etaF2+kap)' * eps * percent_adjust; mudown12 = [mudown12;0]; 
        muup12 = transmu(1:end-1,t-1) .* [2*eta02,etaL2]' * eps * percent_adjust; muup12 = [0;muup12];
        mustay12 = transmu(1:end,t-1) .* (1-[2*eta02,etaF2+kap+[etaL2,0]]*eps * percent_adjust)';
        transmu(:,t)= mustay12 + mudown12 + muup12;
    else
        transmu(:,t) = muvec1;
    end
    vLovF(t) = vL2*transmu(:,t) / (vF2*transmu(:,t));
    if PDA_leader == 0
        transg(t) = gen_compute_g(transmu(:,t)',xvec2,lamb,kap); % \dot{Q}_F(t) under the Spillover from Followers PDA with k = 1
    else
        transg(t) = (transmu(2:n+1,t)'*[xvec2(n+2:2*n) 0]'+transmu(2,t)*(kap+xvec2(n)))*log(lamb); % \dot{Q}_F(t) under the Spillover from Leaders PDA with k = 1
    end
                    
    LI(t) = pishrvec*transmu(:,t);
    mu(t) = transmu(:,t)'*(0:n)';
end

for t=1:nperiods

if t < nperiods
        dotmu0(t) = (transmu(1,t+1)-transmu(1,t))/eps;
    else
        dotmu0(t) = dotmu0(t-1);
    end
    comp(t) = dotmu0(t)*(log(2^(sig/(sig-1))*l0) - (sig/(sig-1))*log((lamb*lL)^((sig-1)/sig) + ...
                                                                                    lF^((sig-1)/sig))) - ...
                        dotmu0(t)*(2*l0 - (lL + lF))/(2*transmu(1,t)*l0 + (1-transmu(1,t))*(lL + lF)); 
    comp(t) = comp(t) * 100; % So that the composition effects are encoded in percent terms, like eta, kappa, and transg
    transg_comp(t) = transg(t) + comp(t);

end

end
% This function is identical to LMS's function of the same name, except
% that:
%
% (i) Lines 10-12 tell MATLAB to alter the guess for the investment
% success rates whenever the new discount rate is less than 50 basis
% points; when our code is executed, this condition will only be satisfied when
% calculating the transition dynamics following a decline in the discount
% rate from 3.59% to 0.33% (the decline in the discount rate in LMS's main
% calibration exercise).  Lines 10-12 imply that
% when solving for the BGP for a discount rate r = 0.33% , the guess
% for investment success rates is
% equal to the solution for investment success rates when solving for the r =
% 3.59% BGP.  This approach is very similar to how LMS make guesses when
% solving for the BGP at low discount rates (see Line 14 of their
% calibration_EMA_submit.m script).  This approach produces a BGP growth
% rate of 0.76% when r = 0.33%, the same as reported in their paper.  
% 
% (ii) Line 34 calculates \dot{Q}_F(t) under the Spillover from Leaders
% PDA (see Definition 3 of our paper).  Line 34 runs if PDA_leader == 1.
% This result is used only for our Online Appendix B.
%
% (iii) Lines 41-54 calculate the productivity growth composition effects
% in Equation 11 of our paper.  transg(t) corresponds to the change in the
% follower productivity index, \dot{Q}_F(t), in period t \in \{1,...,T\}.  
% dotmu0(t) is an approximation to  \dot{\mu}_0(t), the derivative of the share of tied industries with respect to time.   
% transg_comp(t) is productivity growth at time t, taking into account
% composition effects.  
