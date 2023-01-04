% This function calculates the expected product age as a function of the
% age of the firm

function ap = ExpectedProductAge_FirmAge(x,tau,Age)

eapNI = ExpectedProductAgeNotInitial(x,tau,Age);

ap = eapNI + (Age - eapNI) * (1-ProbNoInitialProductAge(x,tau,Age)) * ...
    AvgInverseNAge(x,tau,Age);

end


function y = AvgInverseNAge(x,tau,Age)

gamma_a = x * (1-exp(-(tau-x)*Age)) / (tau - x*exp(-(tau-x)*Age));

y = (1-gamma_a)/(gamma_a) * log(1/(1-gamma_a));

end

function p = ProbNoInitialProductAge(x,tau,Age)

p = 1 - (tau*exp(-x*Age) - x*exp(-tau*Age))/(tau - x);

end

function eap = ExpectedProductAgeNotInitial(x,tau,Age)

Num1 = 1/tau * (1-exp(-tau*Age));
Num2 = 1/x * exp(-tau*Age)*(1-exp(-x *Age));
Den = 1 - exp(-(x+tau)*Age);

eap = (Num1 - Num2) / Den;

end