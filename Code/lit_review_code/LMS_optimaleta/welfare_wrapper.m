function welfare = welfare_wrapper(kap, lamb, pivec, rvec, c, rhovec, sig, n, flatpi, xinit)

    [xvec_mod, muvec_mod, ~, gvec_mod,~,flag_mod] = gen_compute_eqm_correct(lamb,pivec,1,kap/100,rvec/100,xinit,c);

    invcost_mod=[c^2 * xvec_mod.^2,0];
    L = sum(muvec_mod.*(flip(invcost_mod(1:(n+1))) + ...
                        invcost_mod((n+1):end) + ...
                        (rhovec.^(1-sig))*(sig - 1)./((1+rhovec.^(1-sig)).*(sig+rhovec.^(1-sig))) + ...
                        (rhovec.^(1-sig))*(sig - 1)./((1+rhovec.^(1-sig)).*(sig*rhovec.^(1-sig) + 1))));


    markup_s = (sig + rhovec.^(1-sig))./(sig - 1);
    markup_ms = (sig*rhovec.^(1-sig) + 1)./((sig-1)*rhovec.^(1-sig));

    term1 = log((rhovec.^(1-sig))./(1 + rhovec.^(1-sig)));
    term2 = sig/(sig-1);
    term3 = markup_s.^((1-sig)/sig);
    term4 = ((lamb.^min(0:n, flatpi)).*(markup_ms)./(rhovec.^(sig-1))).^((1-sig)/sig);

    lnY0_alt = sum(muvec_mod.*(term1 + term2*log(term3 + term4)));

    welfare = -((lnY0_alt - L)/(rvec/100) + gvec_mod/((rvec/100)^2));

end