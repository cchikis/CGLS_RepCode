% This function calculates the average log mark-up by firm age



function lnmu = Avglnmu_Age(Age, I, tau, x, lambda)
% Expected Product Age
ap = ExpectedProductAge_FirmAge(x,tau,Age);

% Expected log mark-up
lnmu = log(lambda)*(1+I*ap);
end