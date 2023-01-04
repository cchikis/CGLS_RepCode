function out = gen_eqm_eqns_state_dep_kappa(xvec,pivec,costparam,kappa,r)
    n = (length(pivec)-1)/2;
    xvec=xvec/r;
    xvec = abs([xvec,0]);
    out = zeros(1,2*n);
    flag = 2;
if flag==1 % plain-vanilla quadratic costs
    for i=2:n
        out(n+i) = pivec(n+1+i)-pivec(n+i) + xvec(n+1+i)^2/2 - xvec(n+i)^2/2 ...
                   - xvec(n+i)*(xvec(n+1-i) + kappa(i) + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa(i-1) );
        if i<n
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa(i)+r) + kappa(i-1)*xvec(n-i+2) + xvec(n+i+1)*xvec(n-i);
        else
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa(i)+r) + kappa(i-1)*xvec(n-i+2);
        end
    end
    out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*xvec(n-1) - xvec(n)*(xvec(n+1)+kappa(1)+r);
    out(n+1) = pivec(n+2) - pivec(n+1) + xvec(n+2)^2/2 - xvec(n+1)^2/2 - xvec(n+1)*(kappa(1)+r);


elseif flag==2 % quadratic costs, different cost function for leader and follower
    a=costparam(1);
    for i=2:n
        if i>2
            out(n+i) = pivec(n+1+i)-pivec(n+i) + a*(xvec(n+1+i)^2/2 - xvec(n+i)^2/2 ...
                       - xvec(n+i)*(xvec(n+1-i) + kappa(i) + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa(i-1) ));
        else
            out(n+i) = pivec(n+1+i)-pivec(n+i) + a*(xvec(n+1+i)^2/2 - xvec(n+i)^2/2) ...
                       - a*xvec(n+i)*(xvec(n+1-i) + kappa(i) + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa(i-1) );
        end
        if i<n
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa(i)+r) + kappa(i-1)*xvec(n-i+2) + xvec(n+i+1)*xvec(n-i);
        else
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa(i)+r) + kappa(i-1)*xvec(n-i+2);
        end
    end
    out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*xvec(n-1) - xvec(n)*(xvec(n+1)+kappa(1)+r);
    out(n+1) = pivec(n+2) - pivec(n+1) + a*xvec(n+2)^2/2 - xvec(n+1)^2/2 - xvec(n+1)*(kappa(1)+r);
end