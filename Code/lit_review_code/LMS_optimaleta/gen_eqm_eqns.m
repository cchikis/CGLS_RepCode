function out = gen_eqm_eqns(xvec,pivec,costparam,kappa,r) % This function is identical to LMS's function of the same name.
    n = (length(pivec)-1)/2;
    xvec=xvec/r;
    xvec = abs([xvec,0]);
    out = zeros(1,2*n);
    flag = 2;
if flag==1 % plain-vanilla quadratic costs
    for i=2:n
        out(n+i) = pivec(n+1+i)-pivec(n+i) + xvec(n+1+i)^2/2 - xvec(n+i)^2/2 ...
                   - xvec(n+i)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa );
        if i<n
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2) + xvec(n+i+1)*xvec(n-i);
        else
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2);
        end
    end
    out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*xvec(n-1) - xvec(n)*(xvec(n+1)+kappa+r);
    out(n+1) = pivec(n+2) - pivec(n+1) + xvec(n+2)^2/2 - xvec(n+1)^2/2 - xvec(n+1)*(kappa+r);
elseif flag==2 % quadratic costs, different cost function for leader and follower
    a=costparam(1);
    for i=2:n
        if i>2
            out(n+i) = pivec(n+1+i)-pivec(n+i) + a*(xvec(n+1+i)^2/2 - xvec(n+i)^2/2 ...
                       - xvec(n+i)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa ));
        else
            out(n+i) = pivec(n+1+i)-pivec(n+i) + a*(xvec(n+1+i)^2/2 - xvec(n+i)^2/2) ...
                       - a*xvec(n+i)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa );
        end
        if i<n
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2) + xvec(n+i+1)*xvec(n-i);
        else
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2);
        end
    end
    out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*xvec(n-1) - xvec(n)*(xvec(n+1)+kappa+r);
    out(n+1) = pivec(n+2) - pivec(n+1) + a*xvec(n+2)^2/2 - xvec(n+1)^2/2 - xvec(n+1)*(kappa+r);
elseif flag==3 % state-dependent quadratic cost
    c1=costparam(1); c2=costparam(2); c3=costparam(3);
    cons = c1*((1:n)/c3).^c2;
    for i=2:n
        if i>2
            out(n+i) = pivec(n+1+i)-pivec(n+i) + cons(i)*xvec(n+1+i)^2/2 - cons(i-1)*xvec(n+i)^2/2 ...
                       - cons(i-1)*xvec(n+i)*(xvec(n+1-i) + kappa + r) + cons(i-2)*xvec(n+i-1)*(xvec(n+2-i) + kappa );
        else
            out(n+i) = pivec(n+1+i)-pivec(n+i) + cons(i)*xvec(n+1+i)^2/2 - cons(i-1)*xvec(n+i)^2/2 ...
                       - cons(i-1)*xvec(n+i)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa );
        end
        if i<n
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2) + xvec(n+i+1)*xvec(n-i);
        else
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2);
        end
    end
    out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*xvec(n-1) - xvec(n)*(xvec(n+1)+kappa+r);
    out(n+1) = pivec(n+2) - pivec(n+1) + cons(i)*xvec(n+2)^2/2 - xvec(n+1)^2/2 - xvec(n+1)*(kappa+r);
elseif flag==4 % different degree of concavity
    a=costparam;
    for i=2:n
        out(n+i) = pivec(n+1+i)-pivec(n+i) + xvec(n+1+i)^a*(a-1)/a - xvec(n+i)^a*(a-1)/a ...
                   - xvec(n+i)^(a-1)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)^(a-1)*(xvec(n+2-i) + kappa );
        if i<n
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + (a-1)/a*(xvec(n-i+2)^a-xvec(n-i+1)^a) ...
                     -xvec(n-i+1)^(a-1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2)^(a-1) + xvec(n+i+1)*xvec(n-i)^(a-1);
        else
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + (a-1)/a*(xvec(n-i+2)^a-xvec(n-i+1)^a) ...
                     -xvec(n-i+1)^(a-1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2)^(a-1);
        end
    end
    out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^a*(a-1)/a -xvec(n)^a*(a-1)/a + xvec(n+2)*xvec(n-1)^(a-1) - xvec(n)^(a-1)*(xvec(n+1)+kappa+r);
    out(n+1) = pivec(n+2) - pivec(n+1) + xvec(n+2)^a*(a-1)/a - xvec(n+1)^a*(a-1)/a - xvec(n+1)^(a-1)*(xvec(n)+kappa+r)+xvec(n+1)*xvec(n)^(a-1);
elseif flag==5 % non-zero marginal costs at zero effort
    a=costparam;
    xvec = sqrt(xvec.^2);
    for i=2:n
        out(n+i) = pivec(n+1+i)-pivec(n+i) + xvec(n+1+i)^2/2 - xvec(n+i)^2/2 ...
                   - (xvec(n+i)+a)*(xvec(n+1-i) + kappa + r) + (xvec(n+i-1)+a)*(xvec(n+2-i) + kappa );
        if i<n
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -(xvec(n-i+1)+a)*(xvec(n+i)+kappa+r) + kappa*(xvec(n-i+2)+a) + xvec(n+i+1)*(xvec(n-i)+a);
        else
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -(xvec(n-i+1)+a)*(xvec(n+i)+kappa+r) + kappa*(xvec(n-i+2)+a);
        end
    end
    out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*(xvec(n-1)+a) - (xvec(n)+a)*(xvec(n+1)+kappa+r);
    out(n+1) = pivec(n+2) - pivec(n+1) + xvec(n+2)^2/2 - xvec(n+1)^2/2 - (xvec(n+1)+a)*(kappa+r) + a*(xvec(n+1)-xvec(n));
elseif flag==6
    for i=2:n
        out(n+i) = pivec(n+1+i)-pivec(n+i) + xvec(n+1+i)^2/2 - xvec(n+i)^2/2 ...
                   - (xvec(n+i)+a)*(xvec(n+1-i) + kappa + r) + (xvec(n+i-1)+a)*(xvec(n+2-i) + kappa );
        if i<n
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                        -xvec(n-i+1)*(xvec(n+i)+kappa+r) - c*(xvec(n+i)+r) + kappa*xvec(n-i+2) + xvec(n+i+1)*(xvec(n-i)+c);
        else
            out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + ( 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) - c*(xvec(n+i)+r) + kappa*xvec(n-i+2));
        end
    end
    out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*(xvec(n-1)+c) - (xvec(n)+c)*(xvec(n+1)+kappa+r);
    out(n+1) = pivec(n+2) - pivec(n+1) + xvec(n+2)^2/2 - xvec(n+1)^2/2 - (xvec(n+1)+c)*(kappa+r) +c*(xvec(n+1)- xvec(n));
end
