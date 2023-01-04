function [pivec,salesL,salesF,pishr,mkp,lL,lF,l0,rhovec] = compute_pi_fast(sigma, lambda, n)
% See the comment on Line 29 below for a description of this function.  
    vec = 0:n;
    for s = vec
       eps = 1;
       rho_init = 1;
       while eps > 0.0001
           rho = ((lambda^(-s))*((sigma + rho_init^(1-sigma)))/(sigma*(rho_init^(1-sigma)) + 1)).^(1/sigma);
           eps = abs((rho - rho_init)/rho_init);
           rho_init = rho;
       end
       rhovec(s+1) = rho_init;
    end
    rhovec_neg = fliplr(rhovec(2:n+1));
    pivec = [(1./(sigma*(rhovec_neg.^(1-sigma)) + 1)), ((rhovec.^(1-sigma))./(sigma + (rhovec.^(1-sigma))))];
    salesL = rhovec.^(1-sigma) ./ (rhovec.^(1-sigma)+1);
    salesF = 1 ./ (rhovec.^(1-sigma)+1);
    mL = (sigma+rhovec.^(1-sigma))/(sigma-1);
    mF = (sigma*rhovec.^(1-sigma)+1)/(sigma-1)./rhovec.^(1-sigma);
    tmp = [mF(end:-1:2),mL];
    pishr = (tmp-1)./tmp;
    mkp=(mL+mF)/2;
  
    lL = pivec(n+1+1).*(mL(2)-1).^-1; % Calculate production cost of leader in a market with a 1 step productivity gap
    lF = pivec(n).*(mF(2)-1).^-1;     % Calculate production cost of follower in a market with a 1 step productivity gap 
    l0 = pivec(n+1).*(mF(1)-1).^-1;   % Calculate production cost of a firm in a tied market
end

% This function is identical to LMS's compute_pi_fast, except Lines 24-26.
% Lines 24-26 follow from Equation IA.4 of our Online Appendix, and provide
% additional output from this function.  
