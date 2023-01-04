function out = gen_eqm_eqns_correct(xvec,pivec,costparam,kappa,r,c)
    n = (length(pivec)-1)/2;
    %xvec=xvec/r;
    xvec = abs([xvec,0]); 
    out = zeros(1,2*n);
    flag = 2;
% This code is used to produce the middle row of Figure IA.5. 
% Lines 2-6 are exactly as in LMS's code, except that Line 3 is commented
% out.  Commenting out Line 3 of this function does
% not change the equations solved or the interpretation of this function's
% output, because Line 12 of gen_compute_eqm_correct is also commented out.
% See our Online Appendix C.1.2 ("Note on `dividing by r'...") for additional discussion.  

% Because Line 6 of LMS's script gen_eqm_eqns.m (and, correspondingly, this script) sets flag = 2,
% statements below associated with "if flag ==1", 
% "if flag == 3," "if flag = 4,"  "if flag = 5," "if flag = 6" are never executed when
% producing the results in our Comment.  We leave
% those statement unmodified; only the statements associated with "if flag
% == 2" are modified.
        
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
            out(n+i) = (pivec(n+1+i)-pivec(n+i))/2/c^2 + a*(xvec(n+1+i)^2/2 - xvec(n+i)^2/2 ...  
                       - xvec(n+i)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa ));
            % Lines 42-44 correspond to our Equation IA.9 in our Result 2.  Lines 48-49, by comparison, are
            % commented out and are from LMS's function.  Relative to LMS's code, pivec here is divided by 
            % 2c^2.  (Note that in LMS's code, pivec is pre-multiplied (in calibration_EMA_submit.m) by c^2.)  
            % out(n+i) = pivec(n+1+i)-pivec(n+i) + a*(xvec(n+1+i)^2/2 - xvec(n+i)^2/2 ...
            %           - xvec(n+i)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa ));               
        else
            % Lines 52-53 also correspond to our Equation IA.9 in our Result 2.
            out(n+i) = (pivec(n+1+i)-pivec(n+i))/2/c^2 + a*(xvec(n+1+i)^2/2 - xvec(n+i)^2/2) ...
                       - a*xvec(n+i)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa );
            %out(n+i) = pivec(n+1+i)-pivec(n+i) + a*(xvec(n+1+i)^2/2 - xvec(n+i)^2/2) ...
            %           - a*xvec(n+i)*(xvec(n+1-i) + kappa + r) + xvec(n+i-1)*(xvec(n+2-i) + kappa );
                   
        end
        if i<n
            % Lines 60-61 correspond to our Equation IA.10 in our Result 2.
            out(n+1-i) = (pivec(n+2-i) - pivec(n+1-i))/2/c^2 + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2) + xvec(n+i+1)*xvec(n-i);
            %out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
            %         -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2) + xvec(n+i+1)*xvec(n-i);           
        else
            % Lines 66-67 correspond to our Equation IA.11 in our Result 2.  
            out(n+1-i) = (pivec(n+2-i) - pivec(n+1-i))/2/c^2 + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
                     -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2);
            %out(n+1-i) = pivec(n+2-i) - pivec(n+1-i) + 1/2*(xvec(n-i+2)^2-xvec(n-i+1)^2) ...
            %         -xvec(n-i+1)*(xvec(n+i)+kappa+r) + kappa*xvec(n-i+2);
        end
    end
    % Lines 73-74 correspond to our Equations IA.12 and IA.13.
    out(n) = (pivec(n+1) - pivec(n))/2/c^2 + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*xvec(n-1) - xvec(n)*(xvec(n+1)+kappa+r);
    out(n+1) = (pivec(n+2) - pivec(n+1))/2/c^2 + a*xvec(n+2)^2/2 - xvec(n+1)^2/2 - xvec(n+1)*(kappa+r);
    %out(n) = pivec(n+1) - pivec(n) + xvec(n+1)^2/2 -xvec(n)^2/2 + xvec(n+2)*xvec(n-1) - xvec(n)*(xvec(n+1)+kappa+r);
    %out(n+1) = pivec(n+2) - pivec(n+1) + a*xvec(n+2)^2/2 - xvec(n+1)^2/2 - xvec(n+1)*(kappa+r);
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
