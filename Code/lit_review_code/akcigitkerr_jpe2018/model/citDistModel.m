
function [fn] = citDistModel(theta,gammaeta,alpha,begin)
%% citation distribution

%theta         = x(1);
%gamma*eta     = x(2);
%alpha         = x(3);  % determines how dist deviates from exponential

stepMax   = 100;
citMax    = 100;

% Citation distribution
M     = 1/theta;
fkn = zeros(stepMax+1,citMax+1);

for k = 0:stepMax
    for n = 0:citMax
        fkn(k+1,n+1) = ((theta*(1 - theta)^k)/(M*(theta + gammaeta*(alpha^k)*(1 - theta))))...
                      *(gammaeta*(alpha^k)*(1 - theta)/(theta + gammaeta*(alpha^k)*(1 - theta)))^n;
    end
end

fn = sum(fkn)';
fn = fn(begin:end)/sum(fn(begin:end));

end















